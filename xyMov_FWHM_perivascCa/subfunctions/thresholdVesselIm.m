
function [threshIm, skeleton]=thresholdVesselIm(rawIm, prefs)

%function created by Kira, April 2018
%
%loads in raw image of vessel (2D) and specified threshold (std)
%thresholds image, filters and removes background
%
%INPUTS-
%rawIm       : 2D image of vessel
%thresh      : to remove background (noise) - std
%OUTPUTS-
%threshIm    : 2D image of vessel, now thresholded
%skeleton    : skeleton for thresholded image

if nargin<2
    prefs.imgThresh = 0.5;
    prefs.clean = 1; 
end 

% Threshold image
%find mean intensity value of image
meanim = nanmean(reshape(rawIm, [1, size(rawIm, 1)*size(rawIm,2)]));
%find std intensity value of image
stdim = nanstd(reshape(rawIm, [1, size(rawIm, 1)*size(rawIm,2)]));
%set any background pixels (less than mean+thresh*std) to 0
rawIm(rawIm<(meanim+(prefs.imgThresh*stdim))) = 0;
%set any pixels which are nans to zero 
rawIm(isnan(rawIm)) = 0; 
%set any signal pixels (> than mean+thresh*std) to 255
rawIm(rawIm >=(meanim+(prefs.imgThresh*stdim))) = 255;
%rename variable as thresholded image, as no longer raw
threshIm = rawIm;

if prefs.clean == 1
    %use inbuilt matlab functions to clean thresholded image
    threshIm = bwmorph(threshIm,'majority',100);
    threshIm = bwmorph(threshIm,'close',100);
    threshIm = imfill(threshIm,'holes');
    threshIm = bwareaopen(threshIm, 1000);
end

%create automatic skeleton of raw image
closure = bwmorph(threshIm,'close',Inf); %set at inf to skeletonise
closure = bwmorph(closure,'thin',Inf);
closure = bwmorph(closure,'spur',50); %remove sticky out pts; range50-70
skeleton = bwmorph(closure,'clean',50); %clean up skeleton; range50-100

end