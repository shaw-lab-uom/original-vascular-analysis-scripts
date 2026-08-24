function [stim] = findRogueStims(stim,fps)

%function created by Kira, March 2019
%Checks for any increases in stim channel which aren't actually a stim -
%just random spikes in signal 

%INPUTS-
%stim : the binarised stim trace
%fps  : frames per second 
%OUTPUTS- 
%stim : new stim trace with any spikes removed 

%% remove any frames at beginning that are 'ON' b4 1s of time has happened

stimOnInd = find(stim(:,[1:ceil(fps)]));
if ~isempty(stimOnInd)
    disp('stim ch on at very start, removing these'); 
    removeMe = find(stim==0,1)-1; 
    stim(:,[stimOnInd(1):removeMe])=0; 
    disp(['removed frames: ', num2str(stimOnInd(1)), '-', ...
        num2str(removeMe)]);
end
clear stimOnInd removeMe; 


%% remove short stims
stimOnFramesAll = [find(stim==1,1), find(diff(stim))];

%find all odd numbers, except 1st one, and add 1 to find real start of
%stim, as diff takes the point before not after 
if size(stimOnFramesAll,2)>2 %catch to check there are multiple stims in there
    
    %the stim on points are all the odd numbers (NB stim off pts are even)
    %don't need to take first stim on pt, as we have that one correct
    %already, as we used the find func, i.e. find(stim==1,1), so start from
    %2nd odd number (i.e. 3) 
    stimOnInd = [3:2:size(stimOnFramesAll,2)-1];
    
    %loop all the stim on points found that need correcting
    for a = 1:size(stimOnInd,2) 
        %as started from 3, the true index is 3 add which a from loop we're
        %up to
        %add 1 to relevant stim on pts 
        stimOnFramesAll((3-1)+((2*a)-1)) = stimOnFramesAll(stimOnInd(a))+1;
    end
    
    %see how long the stim on and off last - must last >1s to be counted as
    %a proper stim 
    for a = 1:floor(size(stimOnFramesAll,2)/2)
           %find length of stim in frames
           currentStimOnFrames = stimOnFramesAll((a*2)-1):stimOnFramesAll(a*2);
           %use fps to see if stim lasts longer than 1s, if not remove this
           %false stim and replace with zeros 
           if fps > 1 %check if xy movie or linescan
               %xy movie
               if size(currentStimOnFrames,2) < (1*fps)
%                    disp('Stim spike being removed');
                   stim(:,currentStimOnFrames) = 0;
               end
           else
               %linescan
               if size(currentStimOnFrames,2) < (1/fps)
%                    disp('Stim spike being removed');
                   stim(:,currentStimOnFrames) = 0;
               end
           end
    end
    
else 
    
    %inform user the exp dir may be mislabelled as stim 
    disp('This expDir is labelled stim, but there dont seem to be any in trace');
    
end %end of catch to check for multiple stims


end %end of function 
