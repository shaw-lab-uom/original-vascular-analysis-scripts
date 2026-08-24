function [avePxsz, linePxs]=pixel4ls(expDir, width, pref)
%
%function to find pixel size for line scan
%function written by Dori, March 2018
%
%INPUTS-
%expDir = experimental directory containing ini file and the ROI coord file
%
%OUTPUTS-
%linePxs = the pixel size for the line scan

if nargin < 3
    pref.lsNm = '*_diam.csv';
end

[ttt] = findFolders(expDir, '*.ini');
ini = ini2struct(ttt{1});
file = findFolders(expDir, '*.coord');

if ~isempty(file)
    
    if size(file,2) == 0
        disp(expDir);
        disp('You need to move the ROI coord file into the exp dir');
        disp('this is the picture of the linescan');
        return
    end
    
    coords = csvread(file{1});
    
    d = hypot(diff(coords(:,1)), diff(coords(:,2)));
    d_tot = sum(d);
    
    %code to load the csv file with the start pixel for diameter or RBCV
    %cropped image
    [csvFile] = findFolders(expDir,pref.lsNm);
    Results = readtable(csvFile{1});
    
    %NB/ imageJ starts from 0, matlab starts from 1 - so add 1
    startPixel = Results.BX + 1;
    
    xq = [1:0.1:20];
    xq=round(xq,1);
    zooms = [1,5.3,8];
    pxsz = [2.96875,0.560142,0.371094];
    vq = interp1(zooms,1./pxsz,xq,'linear','extrap');
    pxs = 1./vq(xq==str2num(ini.x_.zoom));
    
    linePxs = d(startPixel:startPixel+width-2).*pxs;
    
    lineLength = pxs*d_tot;
    avePxsz = lineLength/length(coords);
    
else
    avePxsz = NaN;
    linePxs = NaN;
end

end %end of function
