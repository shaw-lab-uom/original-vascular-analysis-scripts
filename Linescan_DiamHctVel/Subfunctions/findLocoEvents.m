function [locoEvents]=findLocoEvents(locomotion,fps,prefs)
%
%function to find locomotion events (i.e. onsets and offsets)
%created by Kira, August 2017, updated April 2018
%
%INPUTS:
%locomotion - the locomotion trace, this must have been sent through
%cleanLoco function previously, so it is thresholded and scaled
%fps        - the frames per second (for timing info)
%prefs      - requires min distance between separate loco events(in frames)
%and whether user wants a graph plotted to show the locomotion trace with
%onsets and offsets marked
%the min distance between loco events is auto set at 3s (*fps to convert to
%frames), and the plot is switched on
%
%OUTPUTS:
%locoEvents - a 3d vector with the loco onset pts in dim 1, the loco offset
%pts in dim 2, and the difference between onset and offset in dim 3 (i.e.
%so you can work out the length of the loco period)
%this is an index, so these pts can be applied to the raw loco data too if
%want to check the thresholding hasn't removed any meaningful data
%NB also as the indices can only be integers (as cannot ask for a frame
%which isn't an integer) - they may be out by up to 0.5 frames (for a fast
%frame (fps) rate this will be minimal)
%
%PREREQUISITES:
%must have run cleanLoco function on the inputted loco trace before
%for this code to work

%set preferences if user doesn't send them in
if nargin < 2
    %without the fps this function will not work
    %inform user, and exit function
    disp('need to input the frames per second...')
    return
elseif nargin < 3
    %auto set the prefs if not specified
    %3 seconds between distinct loco events
    if fps > 1
        prefs.minDist = round(3 * fps);
    else %linescans
        prefs.minDist = round(3 / fps);
    end
    %plot the loco trace with onsets and offsets marked
    prefs.plotFlag = 0;
    %flag for if want to remove the flickers in the signal (not loco as
    %shorter than 1s) 1 = get rid of flickers, 0 = keep them
    prefs.flickerFlag = 1; 
end

%get index of locomotion periods (i.e. signal above 0)
%this will only give integers, but loco may start at a non integer
locoIndex = find(locomotion>0);
%check the time dimension is 2nd
if size(locoIndex,1)>size(locoIndex,2)
    locoIndex = locoIndex';
end
%find the differential between loco periods
locoDiff = [1,diff(locoIndex)];
%find periods of locomotion, also separated by >min distance threshold
locoInd_onset_md = [1,find(locoDiff>=prefs.minDist)];

%catch for if no loco events
if size(locoInd_onset_md,2) == 1 && size(locoIndex,2) < round(prefs.minDist)
    
    disp('no distinct loco events detected');
    locoEvents = [];
    locoOnset = [];
    locoOffset = [];
 
    %in case only 1 loco event detected
elseif size(locoInd_onset_md,2) == 1 && size(locoIndex,2) >= round(prefs.minDist)

    locoOnset = locoIndex(1);
    locoOffset = locoIndex(end); 
    locoEvents(1,:) = locoOnset;
    locoEvents(2,:) = locoOffset;
    locoEvents(3,:) = locoEvents(2,:)-locoEvents(1,:);
    
else %more than 1 loco event
    
    %find onset of loco
    %compensate for the 1 added in the differentiating
    locoOnset = locoIndex(locoInd_onset_md)-1;
    clear locoInd_onset_md;
    %find the locomotion offset points
    %catch for if the first offset pt is 1
    %recreate the loco index (with points separated by the min distance)
    %for the offsets, i.e. cant start at 1, and add final point to end
    %(rather than beginning) of vector
    %create temporary index
    locoInd_offset_ttt = find(locoDiff>=prefs.minDist);
    for b = 1:size(locoInd_offset_ttt,2) %loop points in temp index
        if locoInd_offset_ttt(b) < 2 %check if the offset pts ask for 1
            %if so replace with nan, as pt 1 cannot be an offset
            locoInd_offset_md(b) = NaN;
        else
            %for all other points subtract 1, as the diff takes the next
            %onset, so want previous value for offset
            locoInd_offset_md(b) = [locoInd_offset_ttt(b)-1];
        end
    end
    %find last point, i.e. final offset
    locoInd_offset_md = [locoInd_offset_md, size(locoDiff, 2)];
    %remove any NaNs from the index
    locoInd_offset_md(isnan(locoInd_offset_md)) = [];
    %only take the indexed points which are separated by the minimum
    %distance
    locoOffset = locoIndex(locoInd_offset_md)+1;
    clear locoInd_offset_md;
    
    %remove any detected loco events from index  that last for less than 1
    %second, as these may just be 'flickers' and not real loco periods
    %which may not be meaningful to cut other data around
    %may not always wish to remove flickers, e.g. in linescan vel code when
    %just looking to detect ALL spikes (as they are noise)
    if prefs.flickerFlag == 1
        if fps > 1 %i.e. xy movies
            %remove anything that lasts less than 1 s
            locoInd_flicker = find(locoOffset-locoOnset<=1*fps);
            locoOnset(locoInd_flicker)  = [];
            locoOffset(locoInd_flicker) = [];
            clear locoInd_flicker;
        else %i.e. line scans
            %remove anything that lasts less than 1 s
            locoInd_flicker = find(locoOffset-locoOnset<=1/fps);
            locoOnset(locoInd_flicker)  = [];
            locoOffset(locoInd_flicker) = [];
            clear locoInd_flicker;
        end
    end
    
    %if user wants to see plot of loco onset/offset detection
    if prefs.plotFlag == 1
        
        %plots the figure to show loco on and off detection
        figure;
        %find size of screen (for figure sizing)
        screenSz=get(0,'Screensize');
        %make fig size of screen
        set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) ...
            screenSz(4)]);
        plot(locomotion);
        hold on;
        %plot the onset and offset pts detected
        %cleaned loco data has been scaled between 0 and 1, so plot the
        %onset and offset around 10% of the max poss height (0.1)
        plot(locoOnset, (ones(size(locoOnset))*0.1),'go', 'LineWidth',2);
        plot(locoOffset, (ones(size(locoOffset))*0.1),'r*', 'LineWidth',2);
        title('locomotion, with onset and offset marked');
        
    end
    
    %put loco onset and offset into one vector for outputting from func
    locoEvents(1,:) = locoOnset;
    locoEvents(2,:) = locoOffset;
    locoEvents(3,:) = locoEvents(2,:)-locoEvents(1,:);
    
end %end of checking if any loco events detected

end %end of function
