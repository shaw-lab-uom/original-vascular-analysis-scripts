
% Code written by Kira Shaw, updated October 2024
% Adapted: manual px/fps entry with pixel/frame fallback; smoother perp
%   lines via skeleton coordinate smoothing; fallback angle when perp line
%   cannot be computed; precomputed line locations for speed.
% Summary:
% 1. Get diameter from FWHM of raw vess image (match same skel pts for calc)
% 2. Get out the intensity of ca signal around vessel:
%    Use diameter and FWHM - then expand pixel search either side of
%    vess edges- and get intensity of ca channel

% NOTE: user can remove these if they don't want their MATLAB workspace
% cleared/all open figs closed
clear all; close all;

% CUSTOM SUBFUNCTIONS NEEDED TO RUN THIS SCRIPT:
% findFolders, loadTifFileIn2Mat, thresholdVesselIm, autoDetectAcqParams,
% ini2struct (only the last two are new - see subfunctions/ folder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USER TO EDIT:
%enter individual experimental directory which contains vessel tif file to process:
% NB due to user input requirements can only do 1 vessel tif file at a time
expDir = cd; %exp dir

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preferences
% USER NEEDS TO EDIT THESE PREFS:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PREFS FOR WHICH CHANNELS TO LOAD, AND WHETHER NEED PERIVASC CA

% Set whether there is accompanying perivascular calcium to extract:
prefs.caCh = 0; %1 if want to extract calcium too, 0 if not
% If so, what is the name of the tif file to load into MATLAB
prefs.caChNm = []; %'*methoxy.tif'; %'*calcium.tif'

% What is the name of the vessel diameter tif file to load into MATLAB
% Note that there always needs to be a vessel (in order to expand search
% either side of vessel for accompanying calcium)
prefs.vessChNm = '*vessel.tif';

% PREFS FOR EXTRACTING PIXEL SIZE INFO
% Set exactly ONE of the three flags below to 1.
% If all are 0: values from prefs.pxSz / prefs.fps below are used directly.
%
% NEW (Aug 2026, brought over from the MAPS GUI): whichever of the two
% flags below points at a file that isn't actually there (no .ini found,
% no Experiment.xml found) now falls back automatically to reading the
% vessel tif's own metadata (its ImageJ/Fiji calibration - unit/finterval/
% spacing fields) instead of erroring out. Z-step (for a z-stack tif) is
% always read this way regardless of which flag is set, since neither the
% .ini nor Experiment.xml paths carry a z-step field.

prefs.pxInfoINI    = 0; % Scientifica scope  - reads .ini file
prefs.pxInfoXML    = 0; % ThorLabs scope     - reads Experiment.xml
prefs.pxInfoManual = 1; % any scope          - dialog box prompts for entry
% WARNING: IF LOADING INI FILE WITH PXSZ/FPS INFO - MAY NEED TO EDIT NAME
% OF VARIABLE FROM INI SHEET, LINES 81 & 83 ***

% Used only when all three flags above are 0:
prefs.pxSz = NaN; % pixel size in metres (converted to um internally)
prefs.fps  = NaN; % frames per second (Hz)

% PREFS FOR LENGTH (IN PIXELS) OF SKELETON LINE TO USE FOR GENERATING
% PERPENDICULAR LINE

%set length of skel line to find normal
%the larger the line, the more accurate the perp line, but
%the more info lost at either end (limits) of skel
prefs.skelLineLength = 8; %pixels

% PREFS FOR SIZE (IN PIXELS) OF LINE EITHER SIDE OF VESSEL EDGES - FOR
% SCANNING PERIVASCULAR CALCIUM SIGNAL

% number of pixels to scan from edge of vessel to inside the vessel
prefs.insidePx = 10; %pixels
% number of pixels to scan from edge of vessel to outside the vessel
prefs.outsidePx = 20; %pixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PREFS FOR OPTIONAL EDITING:
%threshold for removing background (std) from raw vessel image (frame1)
prefs.imgThresh = 0.5; %std
%pref for vessel thresh code - to clean the thresholded data
prefs.clean = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load fps / px sz / z-step info
%expDir   : experiment directory with ini file in
%fps      : frames per second
%pxsz     : pixel size
%zstep_um : z-step (microns) - only meaningful if the vessel tif happens to
%           be a z-stack; NaN otherwise. Always read straight from the
%           tif's own metadata (see autoDetectAcqParams), since neither
%           the .ini nor Experiment.xml branches below carry a z-step
%           field.

% locate the vessel tif now (cheap - just a directory lookup) purely so we
% can read its metadata here; the expensive load happens later as before
vessTifPath = findFolders(expDir, prefs.vessChNm);
vessTifPath = vessTifPath{1};
[pxsz_auto, fps_auto, zstep_um, pxSrcAuto, fpsSrcAuto] = ...
    autoDetectAcqParams(vessTifPath, expDir);
if ~isnan(zstep_um)
    disp(['Z-step read from TIF metadata: ' num2str(zstep_um) ' um.']);
end

if prefs.pxInfoINI == 1 % Scientifica - read .ini file

    find_ini_file = findFolders(expDir, '*.ini');
    if isempty(find_ini_file)
        disp(['WARNING: prefs.pxInfoINI is set but no .ini file was found under ' ...
            expDir ' - falling back to the TIF''s own metadata (pxsz: ' ...
            pxSrcAuto ', fps: ' fpsSrcAuto ') instead.']);
        pxsz_um = pxsz_auto;
        fps     = fps_auto;
    else
        ini_file = ini2struct(cell2mat(find_ini_file(1)));
        % NOTE THESE PIXEL SIZE NAMES ARE COMPATIBLE FOR SCISCAN INI FILE NAMES
        % MAY NEED TO EDIT PIXEL SIZE NAME AFTER ini_file. .... TO MATCH
        % INIFILE OUTPUTTED IF USING DIFF HARDWARE
        pxsz     = str2num(ini_file.x_.x0x2epixel0x2esz);   % metres
        pxsz_um  = pxsz * 1000000;                           % convert to um
        fps      = str2num(ini_file.x_.frames0x2ep0x2esec);  % Hz
    end

elseif prefs.pxInfoXML == 1 % ThorLabs - read Experiment.xml

    xmlPath = fullfile(expDir, 'Experiment.xml');
    if ~isfile(xmlPath)
        disp(['WARNING: prefs.pxInfoXML is set but Experiment.xml was not found under ' ...
            expDir ' - falling back to the TIF''s own metadata (pxsz: ' ...
            pxSrcAuto ', fps: ' fpsSrcAuto ') instead.']);
        pxsz_um = pxsz_auto;
        fps     = fps_auto;
    else
        S        = readstruct(xmlPath);
        pxsz_um  = S.LSM.pixelSizeUMAttribute;   % um
        fps      = S.LSM.frameRateAttribute;     % Hz
    end

elseif prefs.pxInfoManual == 1 % any scope - user types values in

    answer = inputdlg( ...
        {'Pixel size (\mum):', 'Frame rate (Hz):'}, ...
        'Enter acquisition parameters', 1, {'', ''});
    if isempty(answer)
        disp('User cancelled parameter entry. Exiting...');
        return;
    end
    pxsz_um = str2double(answer{1}); % um
    fps     = str2double(answer{2}); % Hz

else % hard-coded in prefs

    pxsz    = prefs.pxSz;
    pxsz_um = pxsz * 1000000; % convert metres -> um
    fps     = prefs.fps;

end % end of px/fps loading

% If either value is missing, fall back to pixels / frames so the analysis
% can still run and produce meaningful (unit-labelled) output.
if isnan(pxsz_um) || pxsz_um <= 0
    usePixels = true;
    pxsz_um  = 1;           % diameter will be in pixels
    diamUnit  = 'pixels';
    disp('WARNING: no valid pixel size - diameter will be in pixels.');
else
    usePixels = false;
    diamUnit  = '\mum';
end

if isnan(fps) || fps <= 0
    useFrames = true;
    fps      = 1;            % time vector will be frame numbers
    timeUnit  = 'frames';
    disp('WARNING: no valid frame rate - time axis will be in frames.');
else
    useFrames = false;
    timeUnit  = 's';
end


%% load calcium channel:
if prefs.caCh == 1
    ca_ch = findFolders(expDir, prefs.caChNm);
    ca_ch = ca_ch{1};
    [rawCa] = loadTifFileIn2Mat(ca_ch);
    clear ca_ch;
end

%% load vessel channel:
vess_ch = findFolders(expDir, prefs.vessChNm);
vess_ch = vess_ch{1};
[rawVess] = loadTifFileIn2Mat(vess_ch);
clear vess_ch;

%% use vessel channel to get skeleton out and the perpendicular line:

num_frames = size(rawVess,1);

%% find skeleton from vessel outline (use 50th vessel tif frame)

rawVess_ttt = squeeze(rawVess(50,:,:));
[~, skeleton] = thresholdVesselIm(rawVess_ttt, prefs);
skelandvess = imfuse(skeleton, rawVess_ttt);

d = dialog('Position',[300 300 250 150],'Name','Draw Skeleton');
txt = uicontrol('Parent',d,'Style','text','Position',[20 80 210 40], ...
    'String','Click to draw skeleton');
btn = uicontrol('Parent',d,'Position',[85 20 70 25],'String','OK', ...
    'Callback','delete(gcf)');
clear d txt btn;
clear skeleton;

skeleton = roipoly(skelandvess);
ttt = bwmorph(skeleton,'close',Inf);
ttt = bwmorph(ttt,'thin',Inf);
skeleton = bwmorph(ttt, 'bridge');

skelandvess = imfuse(skeleton, rawVess_ttt);

figure; title('vessel outline with final skeleton');
imagesc(skelandvess);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = dialog('Position',[300 300 250 150],'Name','Define vessel regions');
txt = uicontrol('Parent',d,'Style','text','Position',[20 80 210 40], ...
    'String','Click to draw mask to crop around chosen branch');
btn = uicontrol('Parent',d,'Position',[85 20 70 25],'String','OK', ...
    'Callback','delete(gcf)');
clear d txt btn;

[masks(:,:),~,~] = roipoly(skelandvess);
close all;

disp('Cropping images across all frames...');

rawVess_mask = zeros(num_frames, size(rawVess,2), size(rawVess,3));
for i = 1:num_frames

    rawVess_ttt  = squeeze(rawVess(i,:,:));
    rawVess_mean = mean(mean(rawVess_ttt));
    if rawVess_mean < 2^15
        rawVess_ttt(~masks) = 0;
    else
        rawVess_ttt(~masks) = 2^15;
    end
    rawVess_mask(i,:,:) = rawVess_ttt;

end

skeleton_ttt = skeleton;
skeleton_ttt(~masks) = 0;
skeleton = skeleton_ttt;

clear rawVess_ttt;
rawVess_ttt  = squeeze(rawVess_mask(10,:,:));
[threshIm_ttt, ~] = thresholdVesselIm(rawVess_ttt, prefs);

maxrawVess   = max(size(rawVess_ttt));
vesselArea   = sum(sum(threshIm_ttt));
lengthskel   = sum(sum(skeleton));
roughdiameter = vesselArea / lengthskel;
clear vesselArea threshIm_ttt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xend, yend] = find(bwmorph(skeleton,'endpoints'));
contour_skelpts = bwtraceboundary(skeleton, [xend(1), yend(1)], 'NW');

nSkelPts = floor((size(contour_skelpts,1)/2) - 1);

% Smooth skeleton coordinates before fitting tangent lines.
% A moving average over twice the fitting window removes binary-image
% staircase artefacts and gives a much more stable perpendicular angle
% along the skeleton without losing position accuracy.
smWin    = prefs.skelLineLength * 2;
nWorkPts = min(nSkelPts + prefs.skelLineLength, size(contour_skelpts,1));
contour_skelpts(1:nWorkPts, 1) = movmean(contour_skelpts(1:nWorkPts, 1), smWin);
contour_skelpts(1:nWorkPts, 2) = movmean(contour_skelpts(1:nWorkPts, 2), smWin);

%% find branch normal line

disp('get the perpendicular line from vess image');

vid_name = 'FigFWHM_vesselScan.avi';
writer   = VideoWriter(fullfile(expDir, vid_name));

% Fallback tracking: when a perpendicular angle cannot be computed, shift
% the last valid angle to the current skeleton centre instead of skipping.
prev_perp_angle = NaN;

for k = 1:nSkelPts - prefs.skelLineLength

    skelb4    = k;
    skelafter = k + (prefs.skelLineLength - 1);

    % Centre point of this skeleton window (where the perp line will cross)
    midIdx   = k + round(prefs.skelLineLength / 2);
    x_center = contour_skelpts(midIdx, 2);
    y_center = contour_skelpts(midIdx, 1);

    x = contour_skelpts(skelb4:skelafter, 2);
    y = contour_skelpts(skelb4:skelafter, 1);

    normlength = ceil(roughdiameter) * 2;

    % ---- Compute perpendicular angle via atan2 ------------------------------
    % atan2(dy, dx) is defined for all vessel orientations including vertical
    % (dx=0) and horizontal (dy=0), so no special-casing or flipping needed.
    dx_vec = x(end) - x(1);
    dy_vec = y(end) - y(1);

    if dx_vec == 0 && dy_vec == 0
        % Degenerate window - all skeleton points identical after smoothing
        lineValid = false;
    else
        vessel_angle    = atan2(dy_vec, dx_vec);
        perp_angle      = vessel_angle + pi/2;
        prev_perp_angle = perp_angle; % update only on genuine success
        lineValid       = true;
    end

    % ---- Fallback: reuse previous angle at current skeleton centre ----------
    if ~lineValid && ~isnan(prev_perp_angle)
        perp_angle = prev_perp_angle; % prev_perp_angle intentionally not updated
        lineValid  = true;
        disp(['Skelpt ' num2str(k) ': degenerate window - using fallback angle.']);
    end

    % ---- Build perp line from angle and centre point -----------------------
    if lineValid
        normx = x_center + linspace(-normlength, normlength, 10000) .* cos(perp_angle);
        normy = y_center + linspace(-normlength, normlength, 10000) .* sin(perp_angle);

        % Clip to image bounds
        xInd = [find(normx < 1), find(normx > maxrawVess)];
        normx(xInd) = [];  normy(xInd) = [];
        yInd = [find(normy < 1), find(normy > maxrawVess)];
        normx(yInd) = [];  normy(yInd) = [];

        if isempty(normx) || isempty(normy)
            lineValid = false;
        end
    end

    % ---- Build binary line image and save -----------------------------------
    if lineValid
        x_ttt = round(normx);
        y_ttt = round(normy);

        lineIm_ttt = zeros(size(squeeze(rawVess(1,:,:))));
        for m = 1:size(x_ttt, 2)
            lineIm_ttt(y_ttt(m), x_ttt(m)) = 1;
        end
        clear m;

        if size(lineIm_ttt,1) > size(squeeze(rawVess(10,:,:)),1)
            lineIm_ttt = lineIm_ttt(1:size(squeeze(rawVess(10,:,:)),1), :);
        elseif size(lineIm_ttt,2) > size(squeeze(rawVess(10,:,:)),2)
            lineIm_ttt = lineIm_ttt(:, 1:size(squeeze(rawVess(10,:,:)),2));
        end

        lineIm(k,:,:) = lineIm_ttt;

    else
        lineIm(k,:,:) = NaN;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create video to show scan along vessel
    if k == 1
        open(writer);
        figure;
    end

    imshow(imfuse(skeleton, squeeze(rawVess_mask(1,:,:))));
    hold on;
    plot(x, y, 'r', 'LineWidth', 2);  % smoothed skeleton window (tangent)
    if lineValid
        plot(normx, normy, 'b', 'LineWidth', 2); % perpendicular line
    end
    axis equal;
    title(['Skelpt:' num2str(k)]);
    frame_ttt = getframe(gcf);
    writeVideo(writer, frame_ttt.cdata);
    hold off;

end % end of skel loop

close(writer);
close;

clearvars -except lineIm rawVess rawCa skeleton masks expDir prefs fps ...
    pxsz_um usePixels useFrames diamUnit timeUnit;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get out diam and perivascular calc intensity:

if prefs.caCh == 1
    maskCa = zeros(size(rawCa));
end
maskVess = zeros(size(rawVess));
for i = 1:size(rawVess,1)
    if prefs.caCh == 1
        rawIm_ttt = squeeze(rawCa(i,:,:));
        rawIm_ttt(~masks) = NaN;
        maskCa(i,:,:) = rawIm_ttt; clear rawIm_ttt;
    end
    rawIm_ttt = squeeze(rawVess(i,:,:));
    rawIm_ttt(~masks) = NaN;
    maskVess(i,:,:) = rawIm_ttt; clear rawIm_ttt;
end
skeleton(~masks) = 0;

if prefs.caCh == 1
    skelandca  = imfuse(skeleton, squeeze(maskCa(10,:,:)));
end
skelandvess = imfuse(skeleton, squeeze(maskVess(10,:,:)));

if prefs.caCh == 1
    figure;
    subplot(211);
    imagesc(imfuse(skelandca, squeeze(lineIm(10,:,:))));
    axis equal;
    title('Calcium channel w/ skel, and perp line for skel pt 10');
    subplot(212);
    imagesc(imfuse(skelandvess, squeeze(lineIm(10,:,:))));
    axis equal;
    title('Vess channel w/ skel, and perp line for skel pt 10');
    pause(0.5);
    saveas(gcf, fullfile(expDir, filesep, 'FigFWHM_EGperpLineScan.fig'));
    close;
end

% Precompute pixel locations along each perp line.
% These are identical across every frame so computing them once here
% (rather than inside the double loop) gives a significant speed-up.
disp('Precomputing perpendicular line pixel locations...');
nSkelLines = size(lineIm, 1);
lineLocations_pre = cell(nSkelLines, 1);
for j = 1:nSkelLines
    lineLocations_pre{j} = find(squeeze(lineIm(j,:,:)));
end

disp('calculation: diameter and (optional) intensity profile of calcium ch per skel pt/frame...');
if prefs.caCh == 1
    calcium = zeros(size(lineIm,1), size(maskCa,1));
end
cont_diam = zeros(size(lineIm,1), size(maskVess,1));

for i = 1:size(maskVess,1) % loop frames

    if mod(i,100) == 0
        disp(['frame: ', num2str(i), '/', num2str(size(maskVess,1))]);
    end

    ImVess_ttt = squeeze(maskVess(i,:,:));
    if prefs.caCh == 1
        ImCa_ttt = squeeze(maskCa(i,:,:));
    end

    for j = 1:size(lineIm,1) % loop skeleton points

        % Use precomputed locations (same for every frame)
        lineLocations = lineLocations_pre{j};

        % FWHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        profile_ttt = ImVess_ttt(lineLocations);
        halfMax = (min(profile_ttt) + max(profile_ttt)) / 2;
        index1  = find(profile_ttt >= halfMax, 1, 'first');
        index2  = find(profile_ttt >= halfMax, 1, 'last');

        if false
            figure; plot(profile_ttt); hold on;
            plot(index1, profile_ttt(index1), 'go');
            plot(index2, profile_ttt(index2), 'ro');
            ind_ttt = find(~isnan(profile_ttt));
            plot([ind_ttt(1), ind_ttt(end)], [halfMax, halfMax], 'k');
            title('intensity profile with FWHM');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ~isempty(index1) && ~isempty(index2)
            cont_diam(j,i) = (index2 - index1) * pxsz_um;
        else
            cont_diam(j,i) = NaN;
        end

        if prefs.caCh == 1
            lineInd  = [[index1-prefs.outsidePx : index1+prefs.insidePx], ...
                        [index2-prefs.insidePx  : index2+prefs.outsidePx]];
            removeMe = find(lineInd > size(lineLocations,1));
            lineInd(removeMe) = [];
            lineInd(lineInd <= 0) = [];
            calcLocations    = lineLocations(lineInd);
            calcprofile_ttt  = ImCa_ttt(calcLocations);
            calcium(j,i)     = nanmean(calcprofile_ttt);

            if false
                newperpLine = zeros(size(ImCa_ttt));
                newperpLine(calcLocations) = 1;
                figure; imagesc(imfuse(newperpLine, ImCa_ttt)); axis equal;
            end

            clear calcprofile_ttt lineInd calcLocations;
        end

    end % end loop skel pts
end % end loop frames

if prefs.caCh == 1
    for b = 1:size(calcium,1)
        calcium_norm(b,:) = calcium(b,:);
        calcium_norm(b,:) = calcium_norm(b,:) - nanmin(calcium_norm(b,:));
        calcium_norm(b,:) = calcium_norm(b,:) / nanmax(calcium_norm(b,:));
        calcium_norm(b,:) = calcium_norm(b,:) - smooth(calcium_norm(b,:),100)';
    end
end

%% remove nans from cont trace
nanInd = find(isnan(cont_diam(:,1)));
cont_diam(nanInd,:) = [];

% time vector - in seconds if fps known, in frame numbers if not
time = [0 : size(cont_diam,2)-1] / fps;

%% plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

screenSz = get(0,'Screensize');

if prefs.caCh == 1

    figure;
    set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) screenSz(4)]);

    subplot(3,4,[1,2,5,6]);
    imagesc(calcium);
    title('Calc Raw: all data');
    xlabel('Frames'); ylabel('skeleton pt');
    subplot(3,4,[9,10]);
    plot(nanmean(calcium,1));
    xlabel('Frames'); ylabel('Intensity, AU');
    xlim([1 size(calcium,2)]);
    title('Calc Raw: average across skel pts');

    subplot(3,4,[3,4,7,8]);
    imagesc(calcium_norm);
    title('Calc Norm: all data');
    xlabel('Frames'); ylabel('skeleton pt');
    subplot(3,4,[11,12]);
    plot(nanmean(calcium_norm,1));
    xlabel('Frames'); ylabel('Intensity, AU');
    xlim([1 size(calcium_norm,2)]);
    title('Calc Norm: average across skel pts');

    pause(0.5);
    saveas(gcf, fullfile(expDir, filesep, 'FigFWHM_calcFluor.fig'));
    close;

end

figure;
set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) screenSz(4)]);
subplot(3,2,[1:4]);
imagesc(cont_diam);
title('Diameter: all data');
xlabel('Frames'); ylabel('skeleton pt');
subplot(3,2,[5:6]);
plot(time, nanmean(cont_diam,1), 'r');
xlabel(['Time (' timeUnit ')']);
ylabel(['Avg diam (' diamUnit ')']);
xlim([time(1) time(end)]);
title('Diam: average across skel pts');
pause(0.5);
saveas(gcf, fullfile(expDir, filesep, 'FigFWHM_vessDiam.fig'));
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save data into exp dir %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Saving as matfile...');
if prefs.caCh == 1
    matfile = fullfile(expDir, 'contData_xyFWHM.mat');
    save(matfile, 'lineIm', 'rawVess', 'rawCa', 'skeleton', 'masks', 'time', ...
        'calcium', 'calcium_norm', 'cont_diam', 'pxsz_um', 'fps', 'zstep_um', 'nanInd', ...
        'prefs', 'usePixels', 'useFrames', 'diamUnit', 'timeUnit', '-v7.3');
else
    matfile = fullfile(expDir, 'contData_xyFWHM.mat');
    save(matfile, 'lineIm', 'rawVess', 'skeleton', 'masks', 'time', 'nanInd', ...
        'cont_diam', 'pxsz_um', 'fps', 'zstep_um', 'prefs', ...
        'usePixels', 'useFrames', 'diamUnit', 'timeUnit', '-v7.3');
end
