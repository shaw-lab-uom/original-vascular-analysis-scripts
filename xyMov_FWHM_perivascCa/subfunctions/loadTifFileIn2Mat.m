
function [loadedTif] = loadTifFileIn2Mat(tifLoc)

%function created by Kira, March 2019
%Read loop rewritten for speed by Kira Shaw with Claude Code, Aug 2026
%(brought over from the MAPS GUI, which had the same bottleneck)
%For 2P data, loads tif file into workspace

%INPUTS
%tifLoc : a string with the exp dir and the name of the tif file
%OUTPUTS
%loadedTif : the tif frames loaded into the matlab workspace

% SPEED NOTE: the old version called imread(tifLoc, k) once per frame.
% Each call reopens the file and walks the IFD chain from directory 1 up
% to directory k, so for a file with N frames that's an O(N^2) directory
% walk - fine for a handful of frames, but for the thousands of frames in
% a fastScan recording that's what could turn a <1 min ImageJ open into a
% many-minutes MATLAB one. This version opens the file once with the Tiff
% class and steps through directories sequentially with nextDirectory(),
% which is O(N) - the way libtiff is meant to be driven for a multi-page
% file. Output class/values are unchanged (still double) so nothing
% downstream (e.g. thresholdVesselIm's nanmean/nanstd, which can error on
% integer input) is affected.

%get image info (from tif file)
info = imfinfo(tifLoc);
width = info(1).Width;    %width of frame (pixels)
height = info(1).Height;  %height of frame (pixels)
num_frames = numel(info); %number of frames (frames)

%Initialize an empty array to store the tiff file in
loadedTif = zeros(num_frames,width,height); %vessel channel

%inform user loading tif files, slow process...
disp('loading tiff files...');

% open the file once and read directories in sequence (see speed note above)
t = Tiff(tifLoc, 'r');
cleanupObj = onCleanup(@() close(t)); %#ok<NASGU> % closes file even if we error out
for k = 1:num_frames % loop through frames and load each one
    loadedTif(k,:,:) = t.read()';
    if k < num_frames
        t.nextDirectory();
    end
end %end of load frames

end %end of function
