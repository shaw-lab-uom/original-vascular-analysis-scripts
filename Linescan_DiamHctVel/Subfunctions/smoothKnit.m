
function [smoothData] = smoothKnit(data,prefs)

%function written May 2019 by Kira and Orla
%function will smooth data which has nans in it, by smoothing all the pts
%of the trace with data in separately and stitching them back together -
%i.e.because smoothing does weird things to nans

if nargin < 2
    prefs.smoothFactor = 0.1; %0.1 is 10% and so on
    prefs.plotSmoothKnit = 0;
end

% %for testing function
% if false
%     data = rand(1,1000);
%     data(:,[5:15,80:400,670:750])=NaN;
%     data_ttt = data;
% end

%find all NaNs in data:
nanInd_ttt = find(isnan(data));

if ~isempty(nanInd_ttt)
    
    %find start of block of NaNs:
    nanInd(1,:) = [nanInd_ttt(1) nanInd_ttt(find(diff(nanInd_ttt)>1)+1)];
    %find end of block of NaNs:
    nanInd(2,:) = [nanInd_ttt(find(diff(nanInd_ttt)>1)) nanInd_ttt(end)];
    
    smoothData = data;
    
    %use nanInd to find when not NaN periods are
    for a = 1:size(nanInd,2)
        if a == 1
            if nanInd(1,a) > 1 %there is a block of data before first NaNs
                smoothData(:,1:nanInd(1,a)-1) = smooth(data(:,1:nanInd(1,a)-1), ...
                    prefs.smoothFactor,'loess');
            end
        else %all the in between points
            if a < size(nanInd,2)
                smoothData(:,nanInd(2,a-1)+1:nanInd(1,a)-1) = ...
                    smooth(data(:,nanInd(2,a-1)+1:nanInd(1,a)-1), ...
                    prefs.smoothFactor,'loess');
            else
                smoothData(:,nanInd(2,a-1)+1:nanInd(1,a)-1) = ...
                    smooth(data(:,nanInd(2,a-1)+1:nanInd(1,a)-1), ...
                    prefs.smoothFactor,'loess');
                if nanInd(2,a) < size(data,2) %there is a block of data after final NaNs
                    smoothData(:,nanInd(2,a)+1:end) = smooth(data(:, ...
                        nanInd(2,a)+1:end), prefs.smoothFactor,'loess');
                end
            end
        end
    end
    
else %there aren't any nans, can just do normal smooth
    
    smoothData = smooth(data, prefs.smoothFactor,'loess');
    
end %end of check for nans

if prefs.plotSmoothKnit == 1
    figure;
    plot(data,'b');
    hold on;
    plot(smoothData,'r');
    title(['SmoothFactor = ', num2str(prefs.smoothFactor)]);
    legend('orig','smooth');
end

end