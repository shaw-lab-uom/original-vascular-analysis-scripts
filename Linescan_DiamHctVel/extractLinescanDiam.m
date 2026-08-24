function extractLinescanDiam(fname,prefs)

%function to find the diameter of a vessel from a linescan
%this is less complex then for xyFWHM as the vessel is scanned in a
%straight (vertical) line, and there will be no branches
%will step through in 40ms windows and find mean intensity, then FWHM

%this code requires these functions in your path:
%findFolders, pixel4ls, avgByTimeWindow, thresholdVesselIm

%written Dec 2017, by Kira & Orla, updated May 2018 by Kira

%INPUTS-
%fname = directory(s) containing diam.tif - preprocessed vessel in imageJ
%vesselCh = specify which PMT channel the vessel is in, 1 = FITC, 2 = TR,
%if no user input, automatically assumes vessel is in channel 2 (TR)
%borderSz = this should be an integer (in pixels) to crop the edges of each
%frame, may be useful if there is noise around the edge (from image
%registration)
%this is automatically set at 0, unless the user changes this value

%% PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if not enough arguments are inputted
if nargin<2
    
     %%%%%%%% TIMING PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %manually send in the linescan timing info by entering numbers into the
    %prefs below :
    %1 if want to read mspline/fps timing info from ini file
    %0 for manual input
    prefs.loadINI = 1; %1 %0
    prefs.pxInfoXML = 0; %pixel info comes from xml (ThorLabs Scope)
    %if want the code to automatically read in the timing info from an ini
    %file (outputted by sciscan software) leave as empty, i.e. []
    %this function is edited to read in a sciscan ini file - will need to
    %edit
    prefs.mspline = [];
    prefs.fps = [];
    
    %window size must be divisble by 4 (as per Drew paper) 
    %this is because the step size for overlapping time window is 1/4 of
    %window size (aka 40ms time window has an overlapping 10ms step)
    prefs.windowSz = 200; %ms %40 recommended by Drew
    %border for cropping edge of each frame (in pixels)
    prefs.borderSz = 0; %default 0 pixel crop
    
    %%%%%%%%PIXEL SIZE PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %if prefs.loadINI (above) is 1, will try and load the pixel size from
    %the outputted ini file, else if it is 0 and user wants to send in
    %known pizel size manually:
    % in microns
    prefs.pxSz = []; %input number here manually if know pixel size %in um
    
    % EXTRA STEP FOR CALCULATING PIXEL SIZE USING FUNCTION (IN EVENT OF
    % SOFTWARE ERRORS):
    % if user needs to recalculate pixel size (due to error with mismatched
    % zoom/outputted pixel size from sciscan software - call this function:
    % Uses known zoom/pixel size combinations (by default specific to
    % SciScan), and requires a csv file which contains start and end pixels
    % of cropped linescan image (to translate to your linescan trajectory
    % path)
    %1 - means you need to call the 'pixel4ls.m' function to work out the
    %pixel size - as when you set it on sciscan, it errored (open
    %subfunction to see more detailed explanation)
    % 0 - means don't call the function and use the pixel size from ini
    % file or manually inputted above ^
    prefs.pxFlag = 1; %0 %1
    %name of CSV file saved with start and end coords (in pixels) of your
    %linescan path for diameter bounding box (used to crop)
    prefs.csvFileNm = '*_diam.csv';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%% CHANNEL PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %name of diameter tif file
    %should be cropped either side of linescan trajectory section which
    %crosses the vessel horizontally
    prefs.tifNm = '*diam.tif';
    %name of loco and stim channel tif files to search:
    %if have no loco or stim channel info, put open brackets and they won't
    %be loaded
    prefs.locoChNm = '*ch_3.tif'; %[]; %'*ch_3.tif'; %[]
    prefs.stimChNm = '*ch_4.tif'; %[]; %'*ch_4.tif'; %[]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%% TRACE CLEAN PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % prefs for calling func to threshold vessel image from background:
    % thresholdVesselIm.m
    % option if signal not very clear - by default this will be switched
    % off, works better on raw images
    prefs.Thresh = 0; %0 do not threshold image, 1 - threshold vessel image
    prefs.imgThresh = 0.5; %std, used for thresholding image
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% find the diameter tif file:

%check for existence of the diam.tif file
find_tif_file = findFolders(fname, prefs.tifNm);

%loop through all directories with tif file
for a = 1:size(find_tif_file,2)
    
    clearvars -except a find_tif_file prefs fname
    
    %find the local folder containing the tif file
    [expDir,~] = fileparts(find_tif_file{1,a});
    
    disp(['Processing: ', num2str(a), '/', num2str(size(find_tif_file,2))]);
    disp(['ExpDir: ', extractAfter(expDir,fname)]);
    
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
    
    %% TIMING INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %reads the .ini file and extracts the pixel size, fps, and mspline
    %this info is used for creating a time vector
    if prefs.loadINI == 1
        find_ini_file = findFolders(expDir, '*.ini');
        x = cell2mat(find_ini_file);
        ini_file = ini2struct(x);
    end
    if prefs.loadINI == 1 && isempty(prefs.fps) %read sciscan ini file
        %this is specific to ini file name from sciscan software
        %edit if using different software and fps info is stored in
        %diff format
        fps = str2num(ini_file.x_.frames0x2ep0x2esec);
    elseif prefs.pxInfoXML == 1 %ThorLabs microscope
        S = readstruct([expDir,filesep,'Experiment.xml']);
        %1/frameRate=" - think this is time in seconds for 1 frame, so to
        %convert to number of frames per second, divide 1 (s) by this value
        fps = 1/S.LSM.frameRateAttribute; %get frames per second info
    else %load manually inputted fps
        fps = prefs.fps;
    end

        % Extracts the mspline (for time vector)/fps:
    if prefs.loadINI == 1 && isempty(prefs.mspline) %read sciscan ini file
        %this is specific to ini file name from sciscan software
        %edit if using different software and mspline info is stored in
        %diff format
        mspline = str2num(ini_file.x_.ms0x2ep0x2eline);
    elseif prefs.pxInfoXML == 1 %ThorLabs microscope
        % fps / number of lines in a frame 
        %multiply by 1000 to convert to ms??? 
        mspline=((1/S.LSM.frameRateAttribute)/S.LSM.pixelXAttribute)*1000; %0;
    else %load manually inputted mspline
        mspline = prefs.mspline;
    end
    lps = 1000/mspline; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% LOADING PIXEL SIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if prefs.loadINI == 1 && prefs.pxFlag == 0
        %load pixel size from ini file
        %NB/ this may be specific to sciscan ini file - edit the name if
        %needed
        pxsz = str2num(ini_file.x_.x0x2epixel0x2esz) * 1000000;
    elseif prefs.loadINI == 1 && prefs.pxFlag == 1 %software outputted wrong px size
        %call separate function to calculate the pixel size for line scan
        disp('Calling pixel4ls func to calc px sz');
        [pxsz, linePxs]= pixel4ls(expDir, width, prefs.csvFileNm); %this is in microns
    elseif prefs.pxInfoXML == 1 %ThorLabs microscope
       % S = readstruct([expDir,filesep,'Experiment.xml']);
        pxsz = S.LSM.pixelSizeUMAttribute;
        %S.LSM.pixelSizeUMAttribute; %get pixel size info
    else %manually inputted pixel size into prefs
        pxsz = prefs.pxsz;
    end
    
    if ~exist('linePxs','var')
        linePxs = width-1; 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% LOAD TIF FILES (for ls RBCV tif, locomotion, stimulation channels) %
    
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
    rawLine = reshape(A, [width, height*num_images]);
    if ~isnan(loco_ch)
        locoLine = reshape(B, [width_raw, height*num_images]);
    end
    if ~isnan(stim_ch)
        stimLine = reshape(C, [width_raw, height*num_images]);
    end

    
    %explaining the data format:
    %width represents the length of the line you drew
    %multiply height x num_images to recreate the line scan image
    %width (i.e. x) is the number of pixels in the line - i.e. the size of
    %the line you drew
    %height (i.e. y) - each row represents a scan through the line
    %the num_images (i.e. frames) multiplied by length of the y axis is
    %equivalent to how many times the line was scanned - this is just how
    %sciscan stores the info, i.e. can put it back into 2 dimensions
    %normalise the image so the values are between 0 & 1
    %then every line scan (i.e. across multiple sessions) will be on the
    %same scale, and comparible
    
    % Diameter tif file:
    % add border crop to each frame before reshape
    if prefs.borderSz>0
        rawIm_ttt = A(prefs.borderSz:width-prefs.borderSz,...
            prefs.borderSz:height-prefs.borderSz,:);
    else
        rawIm_ttt=A;
    end
    
    %reorder the data:
    %put all frames on top of each other, so it becomes a 2D image
    %get the time x space line scan info
    rawLine = reshape(A, [width, height*num_images]);
    if ~isnan(loco_ch)
        locoLine = reshape(B, [width_raw, height*num_images]);
    end
    if ~isnan(stim_ch)
        stimLine = reshape(C, [width_raw, height*num_images]);
    end
    
    % Normalize the image so the values are between 0 & 1
    % Every line scan (i.e. across multiple sessions) will be on the same
    % scale
    rawLine = rawLine-min(min(rawLine));
    
    % calculate window size in frames from prefs millisecond requested size
    windowsize = round((prefs.windowSz/mspline)/4)*4; %pixels
    

    %% Extract diameter trace:
    
    stepsize = round(.25*windowsize); %divide the no of pixels by 4
    nlines = size(rawLine,2); %time, i.e. no Lines
    npoints = size(rawLine,1); %space, i.e. Pixels
    nsteps = floor(nlines/stepsize)-3;
    
    %loop through frames/rows of the image to find full width half max
    for k = 1:nsteps %step through data
        
        if (mod(k,1000) == 0)
            disp(['Step ', num2str(k), '/', num2str(nsteps)]);
        end
        
        %used to create time vector (sampling frequency)
        timepts(k) = 1+(k-1)*stepsize+windowsize/2;
        
        % take info from Continuous data
        % loop through data in steps of step size * 4
        % this will overlap by one time box on each loop
        % 0:T/4 time box = stepsize+windowsize/2 - i.e. first box
        % 4 steps later (i.e. 4 time boxes on, 4*stepsize) -
        % (k-1)*stepsize+windowsize
        clear data_hold_ttt data_hold;
        data_hold_ttt = ...
            rawLine(:,1+(k-1)*stepsize:(k-1)*stepsize+windowsize)';
        
        %replace any '0' diameter calculations with NaNs
        %as zeros will affect the diameter averages
        for p = 1:size(data_hold_ttt, 1)
            if nanmean(data_hold_ttt(p ,:)) == 0
                data_hold_ttt(p ,:) = nan;
            end
        end
        
        %call a function to threshold the image for getting diameter
        if prefs.Thresh == 1
            if k == 1
                disp('Thresholding image for FWHM.. slow...');
            end
            [data_hold, ~] = thresholdVesselIm(data_hold_ttt, ...
                prefs.imgThresh);
        else
            data_hold = data_hold_ttt;
        end
        
        %for the first window - output a figure of the vessel
        % image for time window specified (used to calc FWHM)
        %if have automated threshold function on, will also display
        %thresholded image
        if k == 1
            figure;
            subplot(211);
            imagesc(data_hold_ttt);
            title('Raw Data');
            if prefs.Thresh == 1
                subplot(212);
                imagesc(data_hold);
                title('Thresholded Data');
            end
            %save figure
            saveas(gcf, fullfile(expDir, 'diam_threshImg.fig'));
            close;
        end
        
        %extract each row of the linescan, for each frame
        data_intensity = nanmean(data_hold,1)';
        
        %  FWHM calculation: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find the half max value.
        halfMax = (min(data_intensity) + max(data_intensity)) / 2;
        %find the nearest pixels to the half way threshold
        [~,I] = sort(abs(halfMax-(data_intensity)),'ascend');
        %find which points of FWHM curve to select (i.e. cnt be next to
        %each other) - or wnt get a diam reading
        for m = 1:size(I,1)
            if abs(I(1) - I(m))>5
                useme = m;
                break
            end
        end
        I = [I(1), I(useme)];
        %take the min val as start of curve, and max as end
        index1 = min(I);
        index2 = max(I);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %check found first and last value after FWHM threshold from intensity curve:
        if ~isempty(index1) && ~isempty(index2)
            % convert the distance between first and last pixels
            % into microns (i.e. for diameter reading) where
            % fluorescence becomes dark background
            if prefs.pxFlag == 1 %used pixel4ls to get pixel size
                diameter(k) = sum(linePxs(index1:index2-1));
            else %have pixel size in um
                diameter(k) = (index2-index1 + 1)*pxsz; %microns
            end
        else %no FWHM calculated from image (due to bad signal)
            diameter(k) = NaN;
        end %end of check if FWHM calculated
        
        %suppressing plot outputs
        %if user is debugging and wants to check that it is working, can
        %use the plots below (and turn off the suppression)
        if false
            %plot the distribution of the data,
            %and the indices for fwhm found
            figure;
            subplot(2,1,1);
            plot(data_hold');
            subplot(2,1,2);
            plot(data_intensity,'b', 'LineWidth',2);
            hold on; 
            plot([index1,index1],[halfMax,halfMax],'go');
            plot([index2,index2],[halfMax,halfMax],'ro');
            title(['window:',num2str(k)]);
        end
        
        
    end %end of stepping through data
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% CLEAN TRACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % check if remove spikes function turned on - on by default:
    if prefs.removeSpikes == 1
        disp('Removing spikes from diam trace');
        diameter = RemoveSpikes(expDir, diameter, prefs);
        if prefs.plotRemovedSpikes
            saveas(gcf, [expDir, filesep,  'diam_SpikesRemovedTrace.fig']);
            close;
        end
    end
    
    % check if smooth function turned on - on by default:
    if prefs.doSmooth == 1
        %call a function to smooth the data (ignoring NaNs)
        [diameter_smooth] = smoothKnit(diameter, prefs);
        if prefs.plotSmoothKnit
            saveas(gcf, fullfile(expDir, 'diam_smoothTrace.fig'));
            close(gcf);
        end
    else
        diameter_smooth = NaN;
    end
    
    % GET TIMING INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate seconds per window, useful for timing info
    %used in place of frames per second for line scan data, as scan style
    %is different
    spw = windowsize*mspline/1000;
    %/1000 to convert from ms to seconds
    time = (timepts*mspline)/1000;
    
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
    
    %continuous data traces - linescan image and diam calculated
    figure;
    screenSz=get(0,'Screensize');
    set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3)/2 ...
        screenSz(4)]);
    
    %line scan diam image plot
    p1=subplot(4,1,1);
    imagesc(time, [], rawLine);
    xlim([time(1),time(end)]);
    xlabel('time(s)'); ylabel('Pixels'); 
    colormap gray
    title('Line Scan Image');
    
    %diameter trace
    p2=subplot(4,1,2);
    plot(time,diameter,'r');
    hold on;
    plot(time,diameter_smooth,'k');
    legend('raw','smooth');
    xlabel('time(s)')
    ylabel('um')
    xlim([time(1),time(end)]);
    
    %line scan image plot
    if ~isnan(loco_ch)
        p3 = subplot(4,1,3);
        plot(time, locomotion);
        xlim([time(1),time(end)]);
        xlabel('time(s)');
        ylabel('A.U.');
        title('Locomotion');
    end
    
    %stim trials plot
    if ~isnan(stim_ch)
        p4=subplot(4,1,4);
        plot(time, stim);
        xlim([time(1),time(end)]);
        xlabel('time(s)');
        title('Stim Trials');
    end
    
    %link x axis so can zoom into raw image and see velocity etc calculated
    if ~isempty(prefs.locoChNm) && ~isempty(prefs.stimChNm)
        linkaxes([p1,p2,p3,p4],'x');
    end
    %save figure
    saveas(gcf, fullfile(expDir, 'diam_contTrace.fig'));
    %close figure
    close;
    
    
    %% save variables
    matfile = fullfile(expDir, 'contData_ls_diam');
    prefs2output = prefs;
    save(matfile,'time','diameter','diameter_smooth', ...
        'pxsz','stim','locomotion','prefs2output','timepts',...
        'fps','mspline', 'spw', 'fps','-v7.3');
    
end %end of looping through tif files

end %end of function

