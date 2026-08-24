function extractLinescanVelocity_adapted(fname, prefs)

% Line scan analysis - use Radon transform on linescan stripes (from RBC
% shadow) to detect RBCV.
% Patrick Drew code from github, adapted by Dori, Kira & Orla, July 2017.
% Updated May 2018 by Kira.
% _adapted version: bug fixes and improvements, August 2026.
%
% Key changes from original:
%   - BUG FIX: deltax now correctly applies pxsz when pxFlag=0, so
%     velocities have correct physical units (mm/s) in all code paths.
%   - BUG FIX: smoothKnit prefs field name corrected (see smoothKnit.m).
%   - Robust percentile normalisation replaces global min-max, so a single
%     bright artefact no longer compresses the whole dynamic range.
%   - Noise angle threshold and minimum velocity cutoff are now prefs
%     (prefs.noiseAngleThresh, prefs.minVel) rather than hardcoded.
%   - wrongWay detection can now be set via prefs.wrongWay instead of
%     requiring a txt file on disk (txt file check is still a fallback).
%   - Frame loading shows percentage progress, not just every 500 frames.
%   - Calls GetVelocityRadon_adapted (improved background subtraction,
%     finer angle resolution, degenerate window check).
%
% Executes the Radon code from:
% https://sites.esm.psu.edu/~pjd17/Drew_Lab/Resources.html
% Reference:
% Drew PJ, Blinder P, Cauwenberghs G, Shih AY, Kleinfeld D, Rapid
% determination of particle velocity from space-time line-scan data
% using the Radon transform, Journal of Computational Neuroscience,
% 29(1-2):5-11
%
% INPUTS
%   fname  - experimental directory(s) containing RBCV.tif file
%   prefs  - optional struct of preferences (defaults set below)
%
% OUTPUTS
%   No variables returned. Results saved to expDir as .mat and figures.
%
% Required functions on path:
%   GetVelocityRadon_adapted, avgDataOverSlidingWindow, findFolders,
%   pixel4ls, findLocoEvents, RemoveSpikes, smoothKnit, cleanLoco,
%   norm_01, findRogueStims

if nargin < 2

    %%%% TIMING PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prefs.loadINI   = 0;   % 1 = read mspline/fps from ini file, 0 = manual
    prefs.pxInfoXML = 1;   % 1 = pixel info from ThorLabs Experiment.xml
    prefs.mspline   = [];  % ms per line (leave [] to read from file)
    prefs.fps       = [];  % frames per second (leave [] to read from file)

    % Window size must be divisible by 4 (Drew paper requirement).
    % 40 ms recommended by Drew; 200 ms used by Rungta for SNR.
    prefs.windowSz  = 200; % ms

    % Border crop applied to each side of each frame (pixels)
    prefs.borderSz  = 0;

    %%%% PIXEL SIZE PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prefs.pxSz      = [];  % um/pixel, set manually if known
    prefs.pxFlag    = 0;   % 1 = call pixel4ls to recompute px size
    prefs.csvFileNm = '*_RBCV.csv';

    %%%% CHANNEL PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prefs.tifNm     = '*RBCV.tif';
    prefs.locoChNm  = [];  % e.g. '*ch_3.tif'; [] = not loaded
    prefs.stimChNm  = [];  % e.g. '*ch_4.tif'; [] = not loaded

    %%%% IMAGE CONDITIONING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Percentile clipping for robust normalisation.
    % Values outside [prefs.normLo, prefs.normHi] percentile are clipped
    % before scaling to [0,1].  Prevents a single bright artefact from
    % compressing the dynamic range of the whole recording.
    prefs.normLo    = 1;   % lower percentile for normalisation
    prefs.normHi    = 99;  % upper percentile for normalisation

    % Image thresholding (background suppression - off by default)
    prefs.Thresh    = 0;
    prefs.imgThresh = 0.5;

    %%%% WRONG-WAY SCAN CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 'auto' - detect from streak slant (default, no user input needed).
    %          A quick Radon sample is taken from the middle of the
    %          recording; if streaks are right-tilting the spatial axis is
    %          flipped (flipud) and wrongWayDetected is set true.
    % true   - always flip (manual override).
    % false  - never flip (manual override).
    % Legacy: a wrongWay.txt file in expDir is still honoured as true.
    prefs.wrongWay  = 'auto';

    %%%% VELOCITY CALCULATION PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Angles within +/- noiseAngleThresh degrees of 90 degrees are treated
    % as noise (near-horizontal lines = unreliably high apparent velocity).
    prefs.noiseAngleThresh = 25; % degrees (original hardcoded value)

    % Velocity values below minVel are set to NaN.
    % Adjust for capillary recordings where true flow < 0.1 mm/s.
    prefs.minVel    = 0.1; % mm/s

    % Charpak scan-direction correction (Rungta et al. 2019, Front. Neurosci.)
    % The Drew formula computes APPARENT velocity (V_app = deltax/deltat * cot(theta)).
    % Because each pixel is acquired sequentially, V_app differs from V_real
    % depending on whether the scanner moves with or against blood flow.
    %
    %   'retrograde'  - scanner moves OPPOSITE to flow direction (recommended).
    %                   Underestimates: V_real = V_app * Vscan / (Vscan - V_app)
    %   'anterograde' - scanner moves in SAME direction as flow.
    %                   Overestimates: V_real = V_app * Vscan / (Vscan + V_app)
    %   'uncorrected' - apply Drew formula without correction (original behaviour).
    %
    % The correction magnitude is V_app/Vscan.  For fast scanners
    % (Vscan >> V_app) the error is small; for slow scanners or fast vessels
    % it can reach tens of percent.
    %
    % NOTE: this is about scan direction RELATIVE TO FLOW, not the spatial
    % flip corrected by prefs.wrongWay.  If wrongWay is applied the effective
    % direction may be swapped — check your setup.
    % 'auto'        - detect from streak slant direction (recommended)
    % 'retrograde'  - override: scanner against flow
    % 'anterograde' - override: scanner with flow
    % 'uncorrected' - no correction (original Drew behaviour)
    prefs.scanDirection = 'auto';

    %%%% TRACE CLEANING PREFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prefs.removeSpikes      = 1;
    prefs.method            = 'mean'; % 'mean','median','gsed','quartile'
    prefs.plotRemovedSpikes = 1;

    prefs.doSmooth          = 1;
    prefs.smoothFactor      = 0.05;
    prefs.plotSmoothKnit    = 0;

    prefs.plotFlag   = 0;
    prefs.minDist    = 100;
    prefs.flickerFlag = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find the RBCV tif file

find_tif_file = findFolders(fname, prefs.tifNm);

for a = 1:size(find_tif_file, 2)

    clearvars -except a find_tif_file prefs fname
    disp([num2str(a), '\', num2str(size(find_tif_file, 2))]);

    [expDir, ~] = fileparts(find_tif_file{1, a});

    info       = imfinfo(find_tif_file{1, a});
    width      = info(1).Width;
    height     = info(1).Height;
    num_images = numel(info);

    %% TIMING INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if prefs.loadINI == 1
        find_ini_file = findFolders(expDir, '*.ini');
        x = cell2mat(find_ini_file);
        ini_file = ini2struct(x);
    end

    if prefs.loadINI == 1 && isempty(prefs.fps)
        fps = str2num(ini_file.x_.frames0x2ep0x2esec);
    elseif prefs.pxInfoXML == 1
        S   = readstruct([expDir, filesep, 'Experiment.xml']);
        fps = 1 / S.LSM.frameRateAttribute;
    else
        fps = prefs.fps;
    end

    if prefs.loadINI == 1 && isempty(prefs.mspline)
        mspline = str2num(ini_file.x_.ms0x2ep0x2eline);
    elseif prefs.pxInfoXML == 1
        mspline = ((1 / S.LSM.frameRateAttribute) / S.LSM.pixelXAttribute) * 1000;
    else
        mspline = prefs.mspline;
    end
    lps = 1000 / mspline;

    % TIMING DIAGNOSTIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sanity check: height * mspline should equal 1000/fps (one frame).
    % If they differ by more than 1%, mspline or fps is wrong and velocity
    % and time vectors will be off.  Common cause: pixelXAttribute in the
    % ThorLabs XML refers to scan-line length (width) rather than the
    % number of scan repetitions per frame (height), or fps is rounded.
    expected_frame_ms = 1000 / fps;
    actual_frame_ms   = height * mspline;
    timing_error_pct  = abs(actual_frame_ms - expected_frame_ms) / expected_frame_ms * 100;
    if timing_error_pct > 1
        warning(['TIMING MISMATCH: height*mspline = %.3f ms, but 1000/fps = %.3f ms ' ...
            '(%.1f%% error).  Check that pixelXAttribute in Experiment.xml ' ...
            'represents scan repetitions per frame (= height = %d), not scan-line ' ...
            'pixel count (= width = %d).  If wrong, set mspline manually in prefs.'], ...
            actual_frame_ms, expected_frame_ms, timing_error_pct, height, width);
    else
        fprintf('  Timing check OK: height*mspline = %.3f ms, 1/fps = %.3f ms (%.2f%% error)\n', ...
            actual_frame_ms, expected_frame_ms, timing_error_pct);
    end

    %% LOADING PIXEL SIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if prefs.loadINI == 1 && prefs.pxFlag == 0
        pxsz = str2num(ini_file.x_.x0x2epixel0x2esz) * 1000000;
    elseif prefs.loadINI == 1 && prefs.pxFlag == 1
        disp('Calling pixel4ls func to calc px sz');
        [pxsz, linePxs] = pixel4ls(expDir, width, prefs.csvFileNm);
    elseif prefs.pxInfoXML == 1
        pxsz = S.LSM.pixelSizeUMAttribute;
    else
        pxsz = prefs.pxSz;
    end

    % BUG FIX: ensure linePxs is always in micrometres so that the
    % downstream deltax calculation is unit-consistent.
    % pixel4ls already returns linePxs in um.  When not called, set
    % linePxs to the physical length of the line in um.
    if ~exist('linePxs', 'var')
        linePxs = (width - 1) * pxsz; % um (scalar)
    end

    %% LOAD TIF FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    otherChFlag = 0;

    if ~isempty(prefs.locoChNm)
        loco_ch       = cell2mat(findFolders(expDir, prefs.locoChNm));
        info_loco     = imfinfo(loco_ch);
        height_ttt    = info_loco(1).Height;
        num_images_ttt = numel(info_loco);
        width_raw     = info_loco(1).Width;
        clear info_loco;
        otherChFlag   = 1;
    else
        disp('No loco channel loaded');
        loco_ch = NaN;
    end

    if ~isempty(prefs.stimChNm)
        stim_ch = cell2mat(findFolders(expDir, prefs.stimChNm));
        if isnan(loco_ch)
            info_stim      = imfinfo(stim_ch);
            height_ttt     = info_stim(1).Height;
            num_images_ttt = numel(info_stim);
            width_raw      = info_stim(1).Width;
            clear info_stim;
            otherChFlag    = 1;
        end
    else
        disp('No stim channel loaded');
        stim_ch = NaN;
    end

    if otherChFlag == 1
        if height ~= height_ttt || num_images ~= num_images_ttt
            disp(['Exp dir: ', expDir]);
            disp('Check num frames in diam vs loco ch... skipping...');
            continue;
        end
    else
        width_raw = NaN;
    end
    clear num_images_ttt height_ttt;

    disp('Loading images...');

    % Pre-allocate arrays (avoids repeated reallocation inside the loop)
    A = zeros(width, height, num_images, 'single');
    if otherChFlag == 1
        B = zeros(width_raw, height, num_images, 'single');
        C = zeros(width_raw, height, num_images, 'single');
    end
    clear otherChFlag;

    reportInterval = max(1, round(num_images / 10)); % report every ~10%
    for k = 1:num_images
        if mod(k, reportInterval) == 0 || k == 1
            disp(sprintf('  Frame %d / %d  (%.0f%%)', k, num_images, ...
                100 * k / num_images));
        end
        A(:, :, k) = single(imread(find_tif_file{1, a}, k)');
        if ~isnan(loco_ch)
            B(:, :, k) = single(imread(loco_ch, k)');
        end
        if ~isnan(stim_ch)
            C(:, :, k) = single(imread(stim_ch, k)');
        end
    end

    rawLine = reshape(A, [width, height * num_images]);
    clear A;
    if ~isnan(loco_ch)
        locoLine = reshape(B, [width_raw, height * num_images]);
        clear B;
    end
    if ~isnan(stim_ch)
        stimLine = reshape(C, [width_raw, height * num_images]);
        clear C;
    end

    %% Extract velocity trace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    windowsize = round((prefs.windowSz / mspline) / 4) * 4; % pixels

    % NORMALISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Robust percentile-based normalisation: clip at prefs.normLo /
    % prefs.normHi before scaling to [0,1].  A single saturated pixel or
    % bleached frame no longer distorts the whole dynamic range.
    p_lo   = double(prctile(rawLine(:), prefs.normLo));
    p_hi   = double(prctile(rawLine(:), prefs.normHi));
    rawLine = double(rawLine);
    rawLine = (rawLine - p_lo) / (p_hi - p_lo);
    rawLine = max(0, min(1, rawLine)); % clip extremes to [0, 1]

    % Diagnostic figure: first 100 time-points after normalisation
    figure;
    imagesc(rawLine(:, 1:min(100, size(rawLine, 2))));
    title('First 100 points of linescan after percentile normalisation');
    colormap gray;
    saveas(gcf, fullfile(expDir, 'vel_EGnormalisation.fig'));
    close;

    % WRONG-WAY AUTO-DETECTION AND CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convention: after any correction, streaks should be LEFT-tilting in
    % the space-time image (retrograde orientation, recommended by Charpak
    % 2019).  LEFT-tilting -> thetas > 0 -> cot(thetas) > 0 -> positive
    % velocity.  RIGHT-tilting means the spatial axis is flipped ("wrong
    % way") and a flipud is needed to restore the expected orientation.
    %
    % Detection uses a QUICK single-window Radon on normalised rawLine so
    % no deltax/deltat are needed — only the SIGN of the angle matters.
    % The sample is taken from the middle of the recording to avoid edge
    % artefacts, with a minimum of one full windowsize of data.

    % --- build sample window ---
    nCols      = size(rawLine, 2);
    samp_len   = min(nCols, windowsize * 5);          % up to 5 windows
    samp_start = max(1, round((nCols - samp_len) / 2) + 1);
    samp_end   = samp_start + samp_len - 1;
    samp_data  = rawLine(:, samp_start:samp_end)';    % transpose: time x space
    samp_data  = samp_data - mean(samp_data(:));      % mean-subtract

    % --- quick Radon angle estimate (coarse + fine, same convention) ---
    angles_qc      = 0:179;
    angles_qf      = -3:0.1:3;
    radon_qc       = radon(samp_data, angles_qc);
    [~, idx_qc]    = max(var(radon_qc));
    coarse_q       = angles_qc(idx_qc);
    radon_qf       = radon(samp_data, coarse_q + angles_qf);
    [~, idx_qf]    = max(var(radon_qf));
    quick_theta    = -1 * ((coarse_q + angles_qf(idx_qf)) - 90);  % same as GetVelocityRadon_adapted

    % --- resolve legacy txt file ---
    find_ww = findFolders(expDir, 'wrongWay.txt');
    legacy_ww = size(find_ww, 2) >= 1;
    if legacy_ww
        disp('  wrongWay.txt found — treating as wrongWay = true (legacy).');
    end

    % --- decide whether to flip ---
    if isequal(prefs.wrongWay, true) || isequal(prefs.wrongWay, 1) || legacy_ww
        % manual override: always flip
        wrongWayDetected = true;
        wrongWayMethod   = 'manual_override';
    elseif isequal(prefs.wrongWay, false) || isequal(prefs.wrongWay, 0)
        % manual override: never flip
        wrongWayDetected = false;
        wrongWayMethod   = 'manual_override';
    else
        % auto: flip if quick_theta is negative (right-tilting = wrong way)
        if abs(quick_theta) < 5
            % angle too close to zero (near-vertical streaks = very slow flow
            % or noisy sample); cannot determine direction reliably
            wrongWayDetected = false;
            wrongWayMethod   = 'auto_uncertain';
            warning(['  Wrong-way auto-detection uncertain (quick_theta = %.1f deg). ' ...
                'Set prefs.wrongWay manually if needed.'], quick_theta);
        elseif quick_theta < 0
            wrongWayDetected = true;
            wrongWayMethod   = 'auto';
        else
            wrongWayDetected = false;
            wrongWayMethod   = 'auto';
        end
    end

    % --- apply flip and report ---
    if wrongWayDetected
        rawLine = flipud(rawLine);
        fprintf('  WRONG WAY DETECTED (method: %s, quick_theta = %.1f deg).\n', ...
            wrongWayMethod, quick_theta);
        disp('  Spatial axis flipped (flipud). Temporal order preserved.');
    else
        fprintf('  Correct orientation confirmed (method: %s, quick_theta = %.1f deg).\n', ...
            wrongWayMethod, quick_theta);
    end

    % GET RBC ANGLE (Radon transform) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Running adapted Radon transform function...');
    [thetasz32, the_tz32, ~] = GetVelocityRadon_adapted(rawLine', windowsize);

    % TIMING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spw    = windowsize * mspline / 1000; % seconds per window
    deltat = spw;

    % REMOVE NOISE SPIKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Angles within +/- noiseAngleThresh of 90 degrees correspond to
    % near-horizontal lines (ambiguous / unreliably high velocity).
    [thetaEvents] = findLocoEvents( ...
        abs(abs(thetasz32) - 90) < prefs.noiseAngleThresh, spw, prefs);

    if ~isempty(thetaEvents)
        thetaOnset  = thetaEvents(1, :);
        thetaOffset = thetaEvents(2, :);
        thetaOnset(thetaOnset == 0) = 1;
        for ev = 1:size(thetaEvents, 2)
            thetasz32(thetaOnset(ev):thetaOffset(ev) + 1) = NaN;
        end
    end

    % VELOCITY CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 1: apparent velocity from Drew et al. (2010) Eq. 4.
    %   V_app = (deltax / deltat) * cot(theta)
    % deltax: total physical length of the linescan line in mm.
    %   sum(linePxs) is always in um (bug fix above ensures this).
    deltax = sum(linePxs) / 1000; % mm

    vel_signed = (deltax / deltat) * cot(deg2rad(thetasz32));

    % Step 2: AUTO-DETECT SCAN DIRECTION FROM STREAK SLANT %%%%%%%%%%%%%%%%
    % The sign of the raw (signed) velocity encodes the relationship between
    % scanner direction and flow direction:
    %
    %   Geometry of rawLine [space x time]:
    %     Right-tilting streaks (bottom-left -> top-right):
    %       RBCs move in SAME direction as scanner -> anterograde -> vel_signed < 0
    %     Left-tilting streaks  (bottom-right -> top-left):
    %       RBCs move AGAINST scanner direction    -> retrograde -> vel_signed > 0
    %
    % We use the median of non-noise, finite values to determine the dominant
    % direction.  prefs.scanDirection = 'auto' enables this; setting it to
    % 'anterograde', 'retrograde', or 'uncorrected' overrides detection.
    %
    % After detection the velocity is expressed as a positive magnitude
    % (physiological convention) so that RemoveSpikes and smoothKnit work
    % correctly.  The detected direction is saved in the output .mat file.

    noise_mask = isnan(thetasz32); % windows already flagged as noise
    finite_vel = vel_signed(~noise_mask & isfinite(vel_signed));

    if strcmpi(prefs.scanDirection, 'auto')
        med_sign = sign(nanmedian(finite_vel));
        if med_sign < 0
            detected_scanDir = 'anterograde';
        elseif med_sign > 0
            detected_scanDir = 'retrograde';
        else
            detected_scanDir = 'uncorrected';
            warning('Could not determine scan direction from streak slant; no correction applied.');
        end
        fprintf('  Auto-detected scan direction: %s (median raw vel sign = %+d)\n', ...
            detected_scanDir, med_sign);
    else
        detected_scanDir = lower(prefs.scanDirection);
    end

    % Convert signed velocity to positive magnitude using detected direction.
    % Anterograde streaks produce negative vel_signed — flip to positive.
    % Retrograde streaks are already positive.
    if strcmpi(detected_scanDir, 'anterograde')
        vel_app = -vel_signed;  % negate so magnitude is positive
    else
        vel_app = vel_signed;   % retrograde or uncorrected: already positive
    end

    % Step 3: scanner speed — needed for Charpak correction.
    % Vscan = scan line length (mm) / time per line (s)
    Vscan = deltax / (mspline / 1000); % mm/s
    fprintf('  Scanner speed Vscan = %.1f mm/s\n', Vscan);

    % Step 4: Charpak scan-direction correction (Rungta et al. 2019).
    % Each pixel is acquired sequentially, so scanner motion biases V_app.
    % Correction magnitude scales with V_app / Vscan — small when Vscan >> V_app.
    switch detected_scanDir
        case 'retrograde'
            % Scanner moves against flow: V_app underestimates V_real.
            % Rungta 2019 Eq. 6: V_real = V_app * Vscan / (Vscan - V_app)
            % Formula diverges when V_app >= Vscan (physically impossible for
            % retrograde); set those windows to NaN.
            vel = vel_app .* Vscan ./ (Vscan - vel_app);
            vel(vel_app >= Vscan) = NaN;
            disp('  Charpak retrograde correction applied.');
        case 'anterograde'
            % Scanner moves with flow: V_app overestimates V_real.
            % V_real = V_app * Vscan / (Vscan + V_app)
            vel = vel_app .* Vscan ./ (Vscan + vel_app);
            disp('  Charpak anterograde correction applied.');
        otherwise % 'uncorrected'
            vel = vel_app;
    end

    % Remove sub-threshold velocities — use MAGNITUDE threshold so negative
    % raw velocities are not silently discarded before direction detection.
    vel(abs(vel) < prefs.minVel) = NaN;

    % Sigma-clipping: remove values > 3 SD from the mean
    m   = nanmean(vel);
    s   = nanstd(vel);
    vel(abs(vel - m) > 3 * s) = NaN;

    velocity = vel';
    vel_app  = vel_app';
    clear vel vel_signed;
    if size(velocity, 1) > size(velocity, 2)
        velocity = velocity';
        vel_app  = vel_app';
    end

    % CLEAN VELOCITY TRACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if prefs.removeSpikes == 1
        disp('Removing spikes from velocity trace...');
        velocity = RemoveSpikes(expDir, velocity, prefs);
        if prefs.plotRemovedSpikes
            saveas(gcf, [expDir, filesep, 'vel_SpikesRemovedTrace.fig']);
            close;
        end
    end

    if prefs.doSmooth == 1
        [velocity_smooth] = smoothKnit(velocity, prefs);
        if prefs.plotSmoothKnit
            saveas(gcf, fullfile(expDir, 'vel_smoothTrace.fig'));
            close;
        end
    else
        velocity_smooth = NaN;
    end

    lineangle = thetasz32'; clear thetasz32;
    timepts   = the_tz32';  clear the_tz32;

    % CREATE TIME VECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    time = (timepts * mspline) / 1000; % seconds

    % Guard against off-by-one length mismatches
    if size(lineangle, 2) > size(time, 2)
        n = size(time, 2);
        lineangle       = lineangle(:, 1:n);
        velocity        = velocity(:, 1:n);
        velocity_smooth = velocity_smooth(:, 1:n);
    end

    %% MATCH LOCO / STIM CHANNELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isnan(loco_ch)
        [rawLoco_tw] = avgDataOverSlidingWindow(locoLine', windowsize);
        locomotion   = rawLoco_tw'; clear rawLoco_tw;
        [locomotion] = cleanLoco(locomotion);
        [locomotion] = norm_01(locomotion);
        locomotion(:, locomotion <= 0.01) = 0;
        if size(locomotion, 2) ~= size(time, 2)
            disp('WARNING: loco trace different size to vel trace');
        end
    else
        locomotion = NaN;
    end

    if ~isnan(stim_ch)
        [rawStim_tw] = avgDataOverSlidingWindow(stimLine', windowsize);
        stim         = rawStim_tw'; clear rawStim_tw;
        stim         = stim > (mean(stim) + std(stim));
        [stim]       = findRogueStims(stim, spw);
        if size(stim, 2) ~= size(time, 2)
            disp('WARNING: stim trace different size to vel trace');
        end
    else
        stim = NaN;
    end

    %% Plot continuous traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure;
    screenSz = get(0, 'Screensize');
    set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) screenSz(4)]);

    p1 = subplot(4, 1, 1);
    imagesc(time, [], rawLine);
    xlim([time(1), time(end)]);
    xlabel('time (s)');
    ylabel('distance along vessel');
    colormap gray;
    title('Line Scan Image');

    p2 = subplot(4, 1, 2);
    plot(time, lineangle, 'b');
    xlim([time(1), time(end)]);
    xlabel('time (s)');
    ylabel('angle (degrees)');

    p3 = subplot(4, 1, 3);
    plot(time, velocity, 'r');
    hold on;
    plot(time, velocity_smooth, 'k');
    legend('velocity', 'velocity smooth', 'AutoUpdate', 'off');
    xlabel('time (s)');
    ylabel('velocity (mm/s)');
    xlim([time(1), time(end)]);

    p4 = subplot(4, 1, 4);
    plot(time, locomotion);
    xlim([time(1), time(end)]);
    xlabel('time (s)');
    ylabel('AU');
    title('Locomotion');

    linkaxes([p1, p2, p3, p4], 'x');
    saveas(gcf, fullfile(expDir, 'vel_contTrace.fig'));
    close(gcf);

    %% Save variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prefs2output = prefs;
    matfile = fullfile(expDir, 'contData_ls_RBCV_200msWindow');
    save(matfile, 'velocity', 'velocity_smooth', 'vel_app', 'Vscan', ...
        'detected_scanDir', 'wrongWayDetected', 'wrongWayMethod', ...
        'lineangle', 'locomotion', 'stim', ...
        'fps', 'mspline', 'time', 'timepts', 'rawLine', 'spw', ...
        'pxsz', 'lps', 'prefs2output', '-v7.3');

end % end of looping tif files

end % end of function
