function [dataOut]=avgDataOverSlidingWindow(dataIn,windowsize)
%
%function to average data traces over a sliding window
%function written by Kira, January 2018, updated May 2018
%
%INPUTS-
%dataIn = the continuous data trace, this should be a 2D vector, and the
%time (or frame) dimensions should be 1st!
%windowsize = the size of the window to average data over
%
%OUTPUTS-
%dataOut = data now averaged over the specified sliding window

%the data applies the sliding time window across the time (or frames)
%dimension
%catch to check that this dimension is 1st
if size(dataIn,2)>size(dataIn,1)
    dataIn=dataIn';
end

%define step size for sliding window
%this is because the data will be averaged across the window, then step 1/4
%of the window and average again - this means data is overlapped in
%averaging by multiple windows, and provides a smoothing effect
%this is adapted from Patrick Drew linescan RBCV code
%divide the no of pixels in the window by 4
stepsize = round(.25*windowsize); %pixels
nlines = size(dataIn,1); %time, i.e. no Lines, 1st dimension
%define number of steps to loop the data
%how many averages to take overall considering the the size of the time
%dimension and the sliding window size
nsteps = floor(nlines/stepsize)-3;

for k = 1:nsteps %loop the steps
    
    %get mean intensity of image across each time window and frame
    %this is done the same way as the getVelocityRadon code, so can match
    %the line scan RBCV data
    dataOut(k,:)=squeeze(mean(mean(dataIn(1+(k-1)*stepsize:(k-1) ...
        *stepsize+windowsize,:),1),2));
    
end %end of steps loop

end %end of function
