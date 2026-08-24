function [pxsz_um, fps, zstep_um, pxSource, fpsSource] = autoDetectAcqParams(tifPath, expDir)
% autoDetectAcqParams  Try to recover pixel size (microns), frame rate (Hz)
% and z-step (microns) for a vessel TIF without asking the user to type
% them in.
%
% Written by Kira Shaw with Claude Code, Aug 2026 - adapted from the MAPS
% GUI's version of this function (same TIF/ini/xml cascade), extended
% here to also read z-step, since neither the .ini nor Experiment.xml
% paths below carry one.
%
% Checks, in order, stopping as soon as both pixel size and fps are found:
%   1. The TIF's own metadata - XResolution + the ImageJ-style
%      ImageDescription block (unit=..., finterval=..., spacing=...) that
%      Fiji writes when it saves a calibrated (hyper)stack.
%   2. A .ini file anywhere under expDir (Scientifica / SciScan rig).
%   3. An Experiment.xml file in expDir (ThorLabs rig).
% Whatever isn't found is left as NaN so the calling script can fall back
% to its own manual-entry / hard-coded-prefs path.
%
%INPUTS
% tifPath : full path to the vessel tif file
% expDir  : experiment directory (folder containing the tif)
%OUTPUTS
% pxsz_um   : pixel size in microns, or NaN if not found anywhere
% fps       : frame rate in Hz, or NaN if not found anywhere
% zstep_um  : z-step in microns - only present for a z-stack tif (i.e.
%             ImageJ's "spacing" field); NaN for a plain xy time-lapse,
%             and NaN if the .ini/xml fallbacks below are used instead
%             (neither carries a z-step field in this script)
% pxSource  : where pxsz_um came from - 'TIF metadata' / 'ini file' /
%             'Experiment.xml' / 'not found'
% fpsSource : same, for fps

pxsz_um   = NaN;
fps       = NaN;
zstep_um  = NaN;
pxSource  = '';
fpsSource = '';

%% 1) the tif's own metadata -------------------------------------------------
try
    info1 = imfinfo(tifPath);
    info1 = info1(1);

    descr = '';
    if isfield(info1, 'ImageDescription')
        descr = info1.ImageDescription;
    end

    % --- pixel size: XResolution is px/unit, "unit" is stated in descr -----
    if isfield(info1, 'XResolution') && info1.XResolution > 0
        if ~isempty(regexpi(descr, 'unit\s*=\s*micron', 'once')) || ...
           ~isempty(regexpi(descr, 'unit\s*=\s*um\b',    'once'))
            pxsz_um  = 1 / info1.XResolution;      % px/um -> um/px
            pxSource = 'TIF metadata';
        elseif isfield(info1, 'ResolutionUnit')
            % no explicit "unit=micron" - fall back to the TIFF standard
            % ResolutionUnit tag (2 = inch, 3 = centimeter)
            switch info1.ResolutionUnit
                case 2
                    pxsz_um  = 25400 / info1.XResolution;
                    pxSource = 'TIF metadata (inch)';
                case 3
                    pxsz_um  = 10000 / info1.XResolution;
                    pxSource = 'TIF metadata (cm)';
            end
        end
    end

    % --- frame rate: ImageJ's finterval is the time between frames (s) -----
    tok = regexpi(descr, 'finterval\s*=\s*([\d.eE+-]+)', 'tokens', 'once');
    if ~isempty(tok)
        fi = str2double(tok{1});
        if fi > 0
            fps       = 1 / fi;
            fpsSource = 'TIF metadata';
        end
    end

    % --- z-step: ImageJ's "spacing" field, microns (z-stacks only) ---------
    tok = regexpi(descr, 'spacing\s*=\s*([\d.eE+-]+)', 'tokens', 'once');
    if ~isempty(tok)
        zs = str2double(tok{1});
        if zs > 0
            zstep_um = zs;
        end
    end
catch ME
    disp(['autoDetectAcqParams: could not read TIF metadata (' ME.message ')']);
end

if ~isnan(pxsz_um) && ~isnan(fps)
    return;
end

%% 2) a Scientifica / SciScan .ini file somewhere under expDir ---------------
if isnan(pxsz_um) || isnan(fps)
    try
        iniFiles = findFolders(expDir, '*.ini');
        if ~isempty(iniFiles)
            ini_file = ini2struct(iniFiles{1});
            % NOTE: field names below match SciScan .ini output. If your
            % rig writes a .ini file with different variable names, edit
            % the two field names below to match (see ini_file in the
            % workspace to check what they're actually called).
            if isnan(pxsz_um) && isfield(ini_file,'x_') && ...
                    isfield(ini_file.x_,'x0x2epixel0x2esz')
                v = str2double(ini_file.x_.x0x2epixel0x2esz) * 1e6; % m -> um
                if v > 0
                    pxsz_um  = v;
                    pxSource = 'ini file';
                end
            end
            if isnan(fps) && isfield(ini_file,'x_') && ...
                    isfield(ini_file.x_,'frames0x2ep0x2esec')
                v = str2double(ini_file.x_.frames0x2ep0x2esec);
                if v > 0
                    fps       = v;
                    fpsSource = 'ini file';
                end
            end
        end
    catch ME
        disp(['autoDetectAcqParams: could not read ini file (' ME.message ')']);
    end
end

%% 3) a ThorLabs Experiment.xml in expDir ------------------------------------
if isnan(pxsz_um) || isnan(fps)
    xmlPath = fullfile(expDir, 'Experiment.xml');
    if isfile(xmlPath)
        try
            S = readstruct(xmlPath);
            % NOTE: as above, edit the field paths below if your rig's
            % Experiment.xml is laid out differently.
            if isnan(pxsz_um) && isfield(S,'LSM') && isfield(S.LSM,'pixelSizeUMAttribute')
                v = S.LSM.pixelSizeUMAttribute;
                if v > 0
                    pxsz_um  = v;
                    pxSource = 'Experiment.xml';
                end
            end
            if isnan(fps) && isfield(S,'LSM') && isfield(S.LSM,'frameRateAttribute')
                v = S.LSM.frameRateAttribute;
                if v > 0
                    fps       = v;
                    fpsSource = 'Experiment.xml';
                end
            end
        catch ME
            disp(['autoDetectAcqParams: could not read Experiment.xml (' ME.message ')']);
        end
    end
end

if isempty(pxSource),  pxSource  = 'not found'; end
if isempty(fpsSource), fpsSource = 'not found'; end

end
