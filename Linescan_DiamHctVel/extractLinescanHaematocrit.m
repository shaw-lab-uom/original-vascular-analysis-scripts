function extractLinescanHaematocrit(fname, prefs)

% created by Kira, February 2020
% analyses RBC density (as a %) from binarised line scan data

% subfunctions needed in path:
% findFolders, avgDataOverSlidingWindow, cleanLoco, norm_01, findRogueStims,
% RemoveSpikes, smoothKnit

% input: fname
% fname is the experimental directory
% this directory should contain a tif file called 'RBCV_binary.tif' (user creates)
% this is just the velocity trace binarised - make sure it is the same num
% of frames as the loco channel - if need to remove bad frames, delete from
% both in preprocessing or after extracting
% this directory should also contain a .ini (notepad) file, with parameter
% info and channel 3 with the loco info

% output: continuous time series data (hct and locomotion), as well as
% extracted linescan channels (so can look for bad frames if need to
% remove)
% and saves graphs (into exp_dir)

%% PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    
    %%%%%%%% TIMING PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %manually send in the linescan timing info by entering numbers into the
    %prefs below :
    %1 if want to read mspline/fps timing info from ini file
    %0 for manual input
    prefs.loadINI = 0; %1 %0
    prefs.pxInfoXML = 1; %pixel info comes from xml (ThorLabs Scope)
    %if want the code to automatically read in the timing info from an ini
    %file (outputted by sciscan software) leave as empty, i.e. []
    %this function is edited to read in a sciscan ini file - will need to
    %edit
    prefs.mspline = [];
    prefs.fps = [];
    
    %window size must be divisble by 4 (as per Drew paper) 
    %this is because the step size for overlapping time window is 1/4 of
    %window size (aka 40ms time window has an overlapping 10ms step)
    prefs.windowSz = 200; %ms 
    %40 recommended by Drew. 
    % %200ms advised by Rungta paper- esp for the slower scan speed of Thor scope. 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%% CHANNEL PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %name of binarised linescan tif file (for calculating haematocrit from)
    prefs.tifNm = 'RBCV_binary.tif';
    %name of loco and stim channel tif files to search:
    %if have no loco or stim channel info, put open brackets and they won't
    %be loaded
    prefs.locoChNm = '*ch_3.tif'; %[]
    prefs.stimChNm = '*ch_4.tif'; %[]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%% TRACE CLEAN PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % prefs for calling function to remove spikes in trace:
    % RemoveSpikes.m
    %call function?
    prefs.removeSpikes = 1; %1 - run removeSpikes function, %0 dont call func
    %outlier removal method:
    prefs.method = 'mean'; % 'mean', 'median', 'gsed', 'quartile'
    %output plots:
    prefs.plotRemovedSpikes = 1; %1 - will show plots of where spikes removed
    
    % Prefs for smoothing trace (will ignore nans so dnt effect smooth):
    % smoothKnit.m
    %call function?
    prefs.doSmooth = 1; %1 - run smoothKnit function, %0 dont call func
    prefs.smoothFactor = 0.05; %how much to smooth by
    prefs.plotSmoothKnit = 0; %1 to plot the smoothed trace and the original trace.
    %this is set to 0 by default, as the smoothed trace is already overlaid
    %over the original trace in the final continuous trace plot outputted
    %from this function
    %note that both the original unsmoothed and smoothed trace will be
    %saved out as variables
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% find the RBCV binary tif file
% NB it is very important that this is binarised as it relies on that in the
% find pixels >0
find_tif_file = findFolders(fname, prefs.tifNm);

%for plots:
screenSz = get(0,'Screensize');

for a = 1:size(find_tif_file,2) %loop all tif files found
    
    % clear variables used throughout loop to stop interference btwn exp
    % dirs
    clearvars -except find_tif_file a fname screenSz prefs;
    % disp progress to user
    disp([num2str(a),'/',num2str(size(find_tif_file,2))]);
    
    % find individual exp dir for saving vars and plots into
    [expDir,~] = fileparts(find_tif_file{1,a});

    % Get the number of frames, the width and height of the binary linescan
    % note the user may have cropped this tif to only include section of
    % the linescan trajectory which goes down centre of vessel - there
    % could thus be user errors, where the crop does not include the full
    % height of the tif (which would create timing errors)
    % find the height of this tif and compare to unedited loco or stim ch
    % tif to check no user errors
    info = imfinfo(find_tif_file{1,a});
    width = info(1).Width;
    height = info(1).Height;
    num_images = numel(info);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% TIMING INFO FOR CREATING WINDOW SIZE TO MEASURE HCT WITHIN:
    
    % looking for .ini file, which contains parameter info
    if prefs.loadINI == 1
        find_ini_file = findFolders(expDir, '*.ini');
        x = cell2mat(find_ini_file);
        ini_file = ini2struct(x);
    end
    
    % Extracts the mspline (for time vector)/fps:
    if prefs.loadINI == 1 && isempty(prefs.mspline) %read sciscan ini file
        %this is specific to ini file name from sciscan software
        %edit if using different software and mspline info is stored in
        %diff format
        mspline = str2num(ini_file.x_.ms0x2ep0x2eline);
    elseif prefs.pxInfoXML == 1 %ThorLabs microscope
        S = readstruct([expDir,filesep,'Experiment.xml']);
        % fps / number of lines in a frame
        %multiply by 1000 to convert to ms???
        mspline=((1/S.LSM.frameRateAttribute)/S.LSM.pixelXAttribute)*1000; %0;
    else %load manually inputted mspline
        mspline = prefs.mspline;
    end
     lps = 1000/mspline; 
    if prefs.loadINI == 1 && isempty(prefs.fps) %read sciscan ini file
        %this is specific to ini file name from sciscan software
        %edit if using different software and fps info is stored in
        %diff format
        fps = str2num(ini_file.x_.frames0x2ep0x2esec);
    elseif prefs.pxInfoXML == 1 %ThorLabs microscope
        %1/frameRate=" - think this is time in seconds for 1 frame, so to
        %convert to number of frames per second, divide 1 (s) by this value
        fps = 1/S.LSM.frameRateAttribute; %get frames per second info
    else %load manually inputted fps
        fps = prefs.fps;
    end
    
    % use mspline info to create window size
    % input prefs specify size of window preferred (40ms by default)
    windowsize = round((prefs.windowSz/mspline)/4)*4;
    
% Get the number of frames, the width and height of the binary linescan
    % note the user may have cropped this tif to only include section of
    % the linescan trajectory which goes down centre of vessel - there
    % could thus be user errors, where the crop does not include the full
    % height of the tif (which would create timing errors)
    % find the height of this tif and compare to unedited loco or stim ch
    % tif to check no user errors
    info = imfinfo(find_tif_file{1,a});
    width = info(1).Width;
    height = info(1).Height;
    num_images = numel(info);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% LOAD TIF FILES (for RBCV binary tif, locomotion, stimulation) %%%%%%

    otherChFlag = 0;
    
 % get loco channel:
    % e.g. so can remove motion artefacts
    if ~isempty(prefs.locoChNm)
        loco_ch = cell2mat(findFolders(expDir, prefs.locoChNm));
        % loco channel info:
        % NB num of frames and height (256) should be the same between binary
        % RBCV channel and loco channel (as this is important for timing info)
        % width can be different due to crop - so need to extract this
        % individually
        info = imfinfo(loco_ch);
        height_ttt = info(1).Height;
        num_images_ttt = numel(info);
        width_raw = info(1).Width;
        clear info;
        otherChFlag = 1;
    else %no loco channel loaded, fill with nan
        disp('No loco channel loaded');
        loco_ch = NaN;
    end
    
    %get stim channel:
    if ~isempty(prefs.stimChNm)
        stim_ch = cell2mat(findFolders(expDir, prefs.stimChNm));
        %if loading stim but no loco ch, can still get height info to check
        %for user errors when cropping
        if isnan(loco_ch)
            info = imfinfo(stim_ch);
            height_ttt = info(1).Height;
            num_images_ttt = numel(info);
            width_raw = info(1).Width;
            clear info;
            otherChFlag = 1;
        end
    else %no stim channel loaded, fill with nan
        disp('No stim channel loaded');
        stim_ch = NaN;
    end

       %check loco and ls diam channel have same amount of time in them
    %(else there has been a user error when cropping the RBCV binary file)
    %NB can only make this extra check if other channels are being loaded
    if otherChFlag == 1 %other channels have been loaded
        if height ~= height_ttt || num_images ~= num_images_ttt
            %warn user of error and exit function:
            disp(['Exp dir: ', expDir]);
            disp('Check num frames in diam vs loco ch... skipping...');
            continue;
        end
    else
        width_raw = NaN; %no other channels loaded
    end
    clear num_images_ttt height_ttt; %unused vars now passed check
    
    % inform user as process can be quite slow
    disp('Loading images...');
    
    % LOAD DATA TIF FILES INTO MATLAB:
    % Initialize an empty array to store the tiff file in.
    A = zeros(width, height, num_images);
    if otherChFlag == 1
        B = zeros(width_raw, height, num_images);
        C = zeros(width_raw, height, num_images);
    end
    clear otherChFlag;
    for k = 1:num_images %loop frames
        % disp every 500 frames to show progress
        if k == 1 || (mod(k,500) == 0)
            disp(['Frame ', num2str(k), '/', num2str(num_images)])
        end
        % load the linescan diameter
        A(:,:,k) = imread(find_tif_file{1,a}, k)';
        if ~isnan(loco_ch)
            % load locomotion
            B(:,:,k) = imread(loco_ch, k)';
        end
        if ~isnan(stim_ch)
            %load stim
            C(:,:,k)=imread(stim_ch,k)';
        end
    end %end of frames loop
    
    %reorder the data:
    %put all frames on top of each other, so it becomes a 2D image
    %get the time x space line scan info
    binaryLine = reshape(A, [width, height*num_images]);
    if ~isnan(loco_ch)
        locoLine = reshape(B, [width_raw, height*num_images]);
    end
    if ~isnan(stim_ch)
        stimLine = reshape(C, [width_raw, height*num_images]);
    end

    
    %% Extract haematocrit trace: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % calculate haematocrit % for multiple time points along data
    stepsize = round(.25*windowsize);
    % time, i.e. no Lines
    nlines = size(binaryLine,2);
    % space, i.e. Pixels
    npoints = size(binaryLine,1);
    nsteps = floor(nlines/stepsize)-3;
    
    % preallocate to save memory
    thresTP = zeros(nsteps,windowsize,npoints);
    hct = zeros(nsteps,1);
    time_ttt = zeros(nsteps,1)';
    % loop for every step
    for k = 1:nsteps
        %specify time window in terms of number of lines 
        %as it isn't just one data point, take the range of data points,
        %and then go in the middle (will later convert this to ms using mspline,
        %and then seconds by /1000)
        if k < nsteps
            time_ttt(:,k+1) = (k-1)*stepsize+windowsize/2;
        end
        % take info from Continuous data
        % loop through data in steps of step size * 4
        % this will overlap by one time box on each loop
        % 0:T/4 time box = stepsize+windowsize/2 - i.e. first box
        % 4 steps later (i.e. 4 time boxes on, 4*stepsize) - (k-1)*stepsize+windowsize
        % extract thresholded data for every time window
        thresTP(k,:,:) = binaryLine(:,1+(k-1)*stepsize: ...
            (k-1)*stepsize+windowsize)';
        % find RBC density % for every time window
        % in binary, black is set to 255... and white fluorescent bkgrnd to
        % 0, so search for 0s
        hct(k) = (size(find(thresTP(k,:,:)==255),1))./ ...
            (size(thresTP(k,:,:),2)*size(thresTP(k,:,:),3))*100;
    end %end of loop steps
    hct = hct';
    clear thresTP stepsize nlines npoints nsteps; 
    
    %% CLEAN TRACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Subfunctions to clean the hct trace:
    %one to remove outliers and one to smooth (will both be on by
    %default)
    
    % check if remove spikes function turned on - on by default:
    if prefs.removeSpikes == 1
        disp('Removing spikes from trace');
        hct = RemoveSpikes(expDir, hct, prefs);
        if prefs.plotRemovedSpikes
            saveas(gcf, [expDir, filesep,  'hct_SpikesRemovedTrace.fig']);
            close;
        end
    end
    
    % check if smooth function turned on - on by default:
    if prefs.doSmooth == 1
        %call a function to smooth the data (ignoring NaNs)
        [hct_smooth] = smoothKnit(hct, prefs);
        if prefs.plotSmoothKnit
            saveas(gcf, fullfile(expDir, 'hct_smooth.fig'));
            close;
        end
    else
        hct_smooth = NaN;
    end
    
    %% GET TIMING INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % create time vector for plots
    % convert to be in seconds (using mspline info)
    time = (time_ttt*mspline)/1000;
    clear time_ttt;
    
    %calculate seconds per window, useful for timing info
    %used in place of frames per second for line scan data, as scan style
    %is different
    spw = windowsize*mspline/1000;
    
    %% MATCH LOCO / STIM CHANNELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %as the line scan data is averaged over the window size determined and
    %in overlapping steps, the data from the other channels should also be
    %manipulated in the same way
    %i.e. so they are the same size for comparing across channels
    %this function requires time in the 1st dim, so transpose data
    %locomotion channel:
    if ~isnan(loco_ch)
        [rawLoco_tw]=avgDataOverSlidingWindow(locoLine',windowsize);
        locomotion = rawLoco_tw'; clear rawLoco_tw;
        % clean loco so below thresh removed, and data is between 0 and 1
        % so can detect loco events off and on later
        [locomotion] = cleanLoco(locomotion);
        [locomotion] = norm_01(locomotion);
        locomotion(:,find(locomotion<=0.01))=0;
        if size(locomotion,2)~=size(time,2)
            disp(['WARNING: loco trace different size to diam trace']);
        end
    else
        locomotion = NaN;
    end
    %stim channel:
    if ~isnan(stim_ch)
        [rawStim_tw] = avgDataOverSlidingWindow(stimLine',windowsize);
        stim = rawStim_tw'; clear rawStim_tw;
        %find when stim is on, binarise stim channel:
        stim = stim>(mean(stim)+std(stim));
        %call function to check if any spikes in stim signal (which last
        %less than 1s)
        %it will replace the spikes with zeros - this must be a binarised
        %stim trace
        [stim] = findRogueStims(stim,spw);
        if size(stim,2)~=size(time,2)
            disp(['WARNING: stim trace different size to diam trace']);
        end
    else
        stim = NaN;
    end
    
    %% General time series plots of extracted data:
    
    % plot the original line scan continuous data, with thresholded data
    % i.e. check the threshold has worked - as this will be used for % calc
    figure;
    %make fig size of screen
    set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) ...
        screenSz(4)]);
    a1=subplot(3,1,1);
    imagesc(time,[],binaryLine);
    colormap gray;
    title('Original Binary Line Scan Data');
    xlabel('Time (s)');
    a2=subplot(3,1,2);
    plot(time, hct);
    hold on;
    if prefs.doSmooth
        plot(time,hct_smooth,'r');
        legend('norm','smooth');
    end
    title('Haematocrit Plot');
    xlabel('Time (s)');
    ylabel('% RBC density');
    linkaxes([a1,a2],'x');
    if ~isnan(loco_ch)
        a3=subplot(3,1,3);
        plot(time, locomotion);
        title('Locomotion Trace');
        xlabel('Time (s)');
        ylabel('A.U.');
        linkaxes([a1,a2,a3],'x');
    end
    %save continuous plots as matlab fig into exp dir:
    figSave = 'hct_contTrace.fig';
    saveas(gcf, fullfile([expDir,filesep,figSave]));
    close;
    
    
    %% save variables as mat file into individual exp dir
    save([expDir,filesep,'contData_ls_Hct'],'fps','spw',...
        'time','hct','hct_smooth', 'locomotion', 'stim');
    
    
end %end of loop tif files

end %end of func

