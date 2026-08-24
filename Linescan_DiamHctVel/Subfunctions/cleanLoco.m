function [locomotion]=cleanLoco(movement,prefs)
%
%function written by Kira, August 2017, updated April 2018
%function to filter and threshold the locomotion data 
%
%INPUTS:
%movement   - this is the raw locomotion trace
%prefs      - prefs.stdThresh is the threshold, if below set the trace
%to zero (i.e. these are 'flickers' not loco) - this is automatically set
%at 2 unless the user specifies otherwise
%OUTPUTS:
%locomotion - this is the locomotion trace, now filtered to remove
%'flickers' in the signal, which are below threshold and thus not deemed
%real loco, and also centred to zero (using diff), and then scaled between
%0 and 1
%requires additional functions: norm_01

if nargin < 2
    %number to multiple std by for the threshold
    %any signal below this value will be set to zero
    prefs.stdThresh = 1.5; 
end

%find differential, i.e. to get rid of any shifts in baseline
%change all the differentials to positive values
locomotion_ttt=abs(diff(movement));
%check the dimensions of the locomotion data are the expected way round
%i.e. with time in the 2nd dimension 
if size(locomotion_ttt,1)>size(locomotion_ttt,2) %check dims
    locomotion_ttt=locomotion_ttt';
end %end of check dims
%finding the diff will lose one value from the vector, add a 0 to the end
%to replace this value
locomotion=[locomotion_ttt,0];
clear locomotion_ttt; % clear temp var

%find threshold for classifying movement as loco or no loco
%threshold value - # * std of the loco trace
%low number for removing little flickers
locoThresh=prefs.stdThresh*nanmean(locomotion);

%changes all values below threshold, and set to 0. 
locomotion(locomotion<=locoThresh)=0; 

%sets the locomotion values between 0 and 1
%requires additional function norm_01
% [locomotion]=norm_01(locomotion);

%check the dimensions of the locomotion data are the expected way round
%i.e. with time in the 2nd dimension 
if size(locomotion,1)>size(locomotion,2) %check dims
   locomotion=locomotion'; 
end %end of check dims

end %end of function
