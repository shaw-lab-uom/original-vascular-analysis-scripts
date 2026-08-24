function [data_out] = norm_01(data_in)
%
%function to normalise data between 0 and 1 
%written by Kira, updated May 2018
%
%INPUTS-
%data_in  = a continuous data trace which you want to scale between 0 and 1
%
%OUTPUTS-
%data_out = the same continuous data trace, now scaled between 0 and 1

%find min value of the data
minVal = min(data_in);
%find max value of the data
maxVal = max(data_in);
%apply formula to scale data between 0 and 1
%data-min/max-min
data_out = (data_in - minVal) / ( maxVal - minVal );

end 

