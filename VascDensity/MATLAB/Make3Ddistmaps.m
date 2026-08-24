% Will find all binarised images in a folder structure, and one by one load
% them and  convert to 3D distance map.

% This uses MIJ to interface Matlab and ImageJ (http://bigwww.epfl.ch/sage/soft/mij/). 
% There are installation instructions here: https://uk.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab
% % Copied here: 
% Installation 
% 1) Put mij.jar into the java directory of Matlab (e.g for Window Machine 'C:\Program Files\MATLAB\R2009b\java\'). 
% 2) Copy also ij.jar (ImageJ) in the java directory of Matlab. Get this file from the ImageJ website: http://rsb.info.nih.gov/ij/ 
% 3) Extend the java classpath to mij.jar, e.g using the Matlab command: javaaddpath 'C:\Program Files\MATLAB\R2009b\java\mij.jar'. 
% 4) Extend the java classpath to ij.jar, e.g using the Matlab command: javaaddpath 'C:\Program Files\MATLAB\R2009b\java\ij.jar'. 
% 5) Copy your entire imagej plugins folder into the matlab java folder
% 6) Start MIJ by running the Matlab command:
% MIJ.start("imagej-plugin-path"); NB this is the folder with your plugins
% inside the matlab java folder

% javaaddpath 'C:\Program Files\MATLAB\R2019b\java\mij.jar'
% javaaddpath 'C:\Program Files\MATLAB\R2019b\java\ij.jar'
% MIJ.start('C:\Program Files\MATLAB\R2019b\java\plugins')

clear all; close all;

%% set the file location to where data is stored:
selpath = 'D:\Dropbox (Brain Energy Lab)\Everything\Manuscripts\Neurophotonics-Methods\ExampleData\vesselStacks\MM1L\Slice1Side1'; %uigetdir;

%% call function to find the data files to be loaded into matlab
%look for all data files to be analysed
findmatfile = findFolders(selpath, 'thresholdedImage.tif');


%% loop through each experimental directory containing the desired file
for a = 1:size(findmatfile, 2)
    
    clearvars -except findmatfile a; 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     STEP ONE:
%     LOAD THE THRESHOLDED IMAGES INTO MATLAB AND GENERATE DISTANCE MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp([num2str(a),'/',num2str(size(findmatfile,2))]);
    
    %save corresponding exp directory into workspace
%     [expDir, ~] = fileparts(findmatfile{a});
    
    filetoopen = ['path=[' findmatfile{a} ']'];
    
    % Open the image in ImageJ
    MIJ.run('Open...', filetoopen); 
    % Invert stack so vessels are black and background white
    MIJ.run('Invert', 'stack');
    % Calculate 3D distance map from inverted stack
    MIJ.run('3D Distance Map', 'map=EDT image=thresholdedImage mask=Same threshold=1');
    % Immport this distance map into Matlab as 3D array    
    bigdistmap = MIJ.getCurrentImage;
    
    % Save this as "bigdistmap.mat" in the same folder as the original
    % image
    [filepathandname,~] = fileparts(findmatfile{a});

    % Do we need to remove any of the stack because masks were removed
    % because too noisy?
    
    maskfilename = fullfile([filepathandname,filesep,'masksRemoved.tif']);
    maskfilenametoopen = ['path=[' maskfilename ']'];
    
    if isfile(maskfilename)
        % File exists so get this stack into Matlab too.
        
        MIJ.run('Open...', maskfilenametoopen);
        masksremoved = MIJ.getCurrentImage;
        %masks to exclude are 255
    else
        % File does not exist so make a zeros file the same size as bigdistmap.
        masksremoved = zeros(size(bigdistmap));
    end
    
    CleanedDistMap = bigdistmap;
    CleanedDistMap(masksremoved==255) = NaN;
    
    save(fullfile([filepathandname, filesep, 'bigdistmap.mat']),...
        'bigdistmap','CleanedDistMap');

    % Close all ImageJ windows: https://forum.image.sc/t/mij-doesnt-exit/1796
    %     MIJ.run('Clear Results');
    %     MIJ.closeAllWindows;
    MIJ.run('Close All');
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     STEP TWO:
    %     GET THE CAPILLARY SPACING INFO FROM THESE MAPS IN UM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %NB should already be calibrated in um (as data comes from an image
    %which has been converted to um x um - may need to check this)
    
    %Remove edges of image in xy (1st 2 dim)
    %remove surface artefacts, by cropping beginning and end of stack...
    %NB this will have already been done somewhat in the preprocessing
%     CleanedDistMap = CleanedDistMap(20:end-20,20:end-20,5:end-5); 
    
    %reshape this temp var to be in 2D - i.e. so get all values
    dist_um = reshape(CleanedDistMap, ...
        [1 size(CleanedDistMap,1)*size(CleanedDistMap,2)*size(CleanedDistMap,3)]);
    %remove zeros
    dist_um(:,find(dist_um==0))=[];
    %put all values in order from low to high um
    dist_um = sort(dist_um, 'ascend');
    
    %in case need to check plot of data distribution
    if false
        %Construct a histogram 
        %and put data into 10 bins
        figure;
        histfit(dist_um,10);
    end
    
    %get the spacing for each percentile from individual stacks here
    %save into relevant experimental folder
    capspacing.label = [50, 60, 70, 75, 80, 90, 95];
    for b = 1:size(capspacing.label,2)
        capspacing.dist_um(b) = prctile(dist_um, capspacing.label(b));
    end
    
    [expDir,~] = fileparts(findmatfile{a});
    
    save([expDir,filesep,'capspacing_prcntiles.mat'],'capspacing','dist_um');
    
end 
 


