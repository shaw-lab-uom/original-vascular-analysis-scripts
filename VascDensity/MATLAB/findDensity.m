

function findDensity(fname)


%Orla Feb 2021
%funtion to find density ONLY. Can input theshold so e.g. can look at the
%density of vessels <5um only.
%must have already ran Devins macro.
prefs.threshFlag= 1 %if want to set prefs.sizeThreshold
prefs.sizeThreshold = 10 % the threshold below which to measure the density e.g 10um
prefs.removeBranches = 0 % the threshold thelow which any diameters are excluded

%find files
find_dirs = findFolders(fname, 'Branch Details.csv');


for n = 1:size(find_dirs,2) %loop all dirs with the data in
    
    clearvars -except fname find_dirs n prefs;
    
    disp(['Looping expDir ', num2str(n), '/', num2str(size(find_dirs,2))]);
    
    %find individual exp dirs for loading and saving data
    [expDir,~] = fileparts(find_dirs{1,n});
    
    %import as a table and as an array
    summaryTable = table2array(readtable([expDir, filesep,'Branch Details.csv']));
    
    %remove any vessels who's radii are below or equal to a specified value (e.e 0um)
    summaryTable(find(summaryTable(:,4)<= prefs.removeBranches),:) = [];
    branchTable_ttt = readtable([expDir, filesep,'Branch Details.csv']);
    
    %volume analysed
    volumeAnalysed = [table2array(branchTable_ttt(1, ...
        find(string(branchTable_ttt.Properties.VariableNames) == 'VolumeAnalysed')))];
    
    
    
    if prefs.threshFlag == 1
        
        summaryTable(find(summaryTable(:,4)> prefs.sizeThreshold/2), :) = [];
        
        
    end % end of setting threshold
    
    
        
    %total length across all branches in m
    branchLength_m = (sum(summaryTable(:,5))/1000)/1000;
    %volume in mm3
    vol_mm3 = volumeAnalysed/1000^3;
    %total density in m per mm3
    density = branchLength_m /vol_mm3;
    
    
    DensityPrefs = prefs;
    matfile = fullfile(expDir, 'VascularDensity');
    save(matfile, 'density', 'branchLength_m', 'vol_mm3', 'summaryTable', 'DensityPrefs')
    
    
    
    
end % end of looping n
end % end of function