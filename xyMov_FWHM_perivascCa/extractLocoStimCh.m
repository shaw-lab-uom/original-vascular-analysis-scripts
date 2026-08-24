
function extractLocoStimCh(fname)

%function created by Kira, March 2019

%exp folder needs to contain tif files with ch3 an ch4 in name
%NB these are automatically in the name from sciscan

%loads in ch3 and ch4 tif files from inputted folder
%finds average intensity of each of these channels to give locomotion and
%stim traces
%should put spont/stim/_VR into the exp folder name

%INPUTS-
%fname       : folder name with 'vessel.tif' file(s) - can be a top dir -
%NB this code needs vessel.tif file in dir to work

%OUTPUTS-
%no vars outputted from the function, but it will save some figures and
%mat files into the experimental directories with tif file in

%functions you need in your path for this func to work:
%findFolders, INIinfo, loadTifFileIn2Mat, cleanLoco, findRogueStims

%find all loco channels, for directories to extract loco stim for
all_dirs = findFolders(fname, '*ch_3.tif');

%% load vessel tif images and image info

for a = 1:size(all_dirs,2)
    
    clearvars -except a all_dirs fname
    
    disp(['Processing file ', num2str(a), '/', num2str(size(all_dirs,2))]);
    disp([num2str(all_dirs{1,a})]);
    
    [dir,~] = fileparts(all_dirs{1,a});
    
    if exist([dir,'\','stimLocoTraces.mat']) == 0
        
        clear loco_ch;
        loco_ch = all_dirs{1,a};
        [expDir,~] = fileparts(loco_ch);
        
        %call function to find the fps and pixel size for this recording
        [fps, pxsz] = INIinfo(expDir);
        pxsz_um = pxsz*1000000; % pixel size (um)
        
        
        if ~isempty(loco_ch) %check if the exp dir has correct tif files in
            
            %% locomotion
            disp('Processing loco...');
            %call func to load tif file with loco in
            [movement] = loadTifFileIn2Mat(loco_ch);
            movement = nanmean(nanmean(movement,2),3)';
            
            %process locomotion
            %if there is VR in the exp dir title, then it will load this and update
            %VR flag to 1
            VRflag = 0;
            if regexp(expDir, '_VR')
                VRflag = 1;
                %position shows position within VR
                position = movement;
                %velocity shows speed moved through VR ??
                velocity = diff(movement,1);
                clear movement;
                movement = velocity;
                %filter movement to remove position info
                movement(movement>20) = 0;
                movement(movement<4) = 0;
                %as it has been differentiated, add a 0 to end to make same size as
                %other variables
                movement = [movement, 0];
                velocity = [velocity, 0];
            end
            
            %process raw loco data
            %clean the loco to remove 'flickers' (below a threshold)
            [locomotion] = cleanLoco(movement);
            %also save raw locomotion with 'flickers' still there (no thresh)
            %i.e. incase the threshold is too stringent etc - can check it is ok
            locomotion_raw = abs(diff(movement));
            %as the movement var has been differentialed, add 0 on end to make
            %right size
            locomotion_raw = [locomotion_raw,0];
            
            %% stim
            disp('Processing stim...');
            
            %check if spont or stim in title - if so - if spont, just make it same
            %size as loco but filled with zeros
            if ~isempty(regexpi(expDir,'stim'))
                
                stim_ch = cell2mat(findFolders(expDir, '*ch_4.tif'));
                
                if ~isempty(stim_ch)
                   
                    [stim] = loadTifFileIn2Mat(stim_ch);
                    stim = nanmean(nanmean(stim,2),3)';
                    
                    %process stim:
                    
                    %find when stim is on, binarise stim channel:
                    stim = stim>(mean(stim)+std(stim));
                    
                    %call function to check if any spikes in stim signal (which last
                    %less than 1s)
                    %it will replace the spikes with zeros - this must be a binarised
                    %stim trace
                    [stim] = findRogueStims(stim,fps);
                    
                else
                    
                    %sciscan error, no channel 4 recorded (so dont have
                    %stim info)
                    stim = nan(size(locomotion));
                    %save warning into dir
                    save([expDir,filesep,'warning-noLocoCh.txt']);
                    
                end
                
            else
                
                stim = zeros(size(locomotion));
                
            end %checking is stim exp in title
            matfile = fullfile(expDir,'stimLocoTraces.mat');
            if VRflag == 0 %if VR have extra vars to save
                save(matfile,'movement','locomotion',...
                    'fps','locomotion_raw','stim');
            else
                save(matfile,'movement','locomotion', 'stim',...
                    'locomotion_raw','velocity','position','fps');
            end
            
            
        else
            
            disp('there are no loco and/or stim ch tif in this folder...');
            
        end %end of check if exp dir has correct tif files in
        
    else
        
        disp('Already extracted, skipping...');
    end
    
end

end %end of function
