function  dataSpikesRemoved = RemoveSpikes(fname, data, prefs)

% Function to take continuous raw data traces and remove any spike outliers 
% NB can take ls RBCV, ls diam or xyCont data - must be continuous and in
% 2D, but it can have more than 1 in 1st dim (e.g. multiple skel pts for
% xyCont)

%Written by Orla and Kira April 2019

%INPUT-
% data : the 2D continuous data you wish to smooth, e.g. ls RBCV, ls diam
% or contXY trace
% prefs : can specifiy the method you wish matlab to use to detect
% outliers, and whether you want it to plot the graph to show b4 and after
% outlier removal

%OUTPUT-
%dataSpikesRemove : data that has been processed and had spikes removed

%% set prefs:

if nargin < 2
    disp('You havent sent in the data for remove spikes func...');
    return; 
elseif nargin < 3
    %preferences automatically specified:
    %the outlier detection method-
    prefs.method = 'mean'; % 'mean', 'median', 'gsed', 'quartile'
    %for more info follow this link:
    %https://uk.mathworks.com/help/matlab/ref/isoutlier.html
    %plots on or off-
    prefs.plotRemovedSpikes = 0; %1 - to plot all figs in func, 0 to not plot
end

%loop the 1st dim of the data to detect outliers for each cont trace, e.g.
%across skel pts
%% find outliers 
for a = 1:size(data, 1)
    
    %normalise trace - for finding large spikes
    dataNorm(a,:) = data(a,:) - nanmin(data(a,:));
    dataNorm(a,:) = data(a,:) / nanmax(dataNorm(a,:));
    dataNorm(a,:) = dataNorm(a,:) - smooth(dataNorm(a,:),100)';
    
    %finds location of outliers. Puts a 1 if its an outlier
    %uses the method specified in prefs, e.g. mean
    clear findOutliers OutlierLocs;
    findOutliers = isoutlier(dataNorm(a,:), prefs.method);
    
    %find the location of the outliers in the raw trace (i.e. non-normalised)
    %and replace with NaNs
    dataSpikesRemoved(a,:) = data(a,:);
    dataSpikesRemoved(a,findOutliers) = NaN;
    
    
    %% can plot the outlier locations over the trace within a single loop
    %if needed, but dont plot every time as it will crash the computer and
    %create too many plots
    
    %     %finds the location of all the 1s
    %     OutlierLocs = find(findOutliers);
    
    
    %plot the normalised and raw data with outliers plotted
    
    %     %plots the selected outliers on the original trace
    %     figure
    %     subplot (211)
    %     plot(dataNorm(a,:))
    %     hold on
    %     plot(OutlierLocs, dataNorm(a,findOutliers), 'ko');
    %     title('Outliers plotted on Normalised trace')
    %     subplot (212)
    %     plot(data(a,:))
    %     hold on
    %     plot(OutlierLocs, data(a,findOutliers), 'ko');
    %     title('Outliers plotted on Raw trace')

    
end %end of looping 1st dim of data (e.g. skeleton pts in xy data)

%% Plot old trace, new trace and overlap - if plot prefs turned on

if prefs.plotRemovedSpikes == 1 %check if want data plotting
    
    figure;
    
    subplot(3,1,1)
    plot(nanmean(data,1), 'r');
    title('raw data')
    yLimits = get(gca,'YLim');
    
    subplot(312)
    plot(nanmean(dataSpikesRemoved,1), 'b')
    title('Spikes removed from data')
    ylim([yLimits]);
    
    subplot(313)
    plot(nanmean(data,1), 'r');
    hold on
    plot(nanmean(dataSpikesRemoved,1), 'b')
    ylim([yLimits]);
    title('Overlap - red = raw data, blue = spikes removed')
    
    
end %end of check if want it plotting 

end %end of function 


