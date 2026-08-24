function [thetas, the_t, spread_matrix] = GetVelocityRadon_adapted(data, windowsize)
%
% Applies a sliding time window over RBCV linescan data.  Within each
% window a Radon transform finds the angle of the RBC streak lines.
% Windows overlap by 3/4 to provide a smoothed angle estimate.
%
% Adapted from Drew et al. (2010), J Comput Neurosci 29:5-11.
% Original adaptation by Kira, Dori and Orla (2017), updated May 2018.
% _adapted version August 2026.
%
% Changes from GetVelocityRadon:
%   - Two-way mean subtraction replaces scalar global-mean subtraction.
%     Removing per-column (spatial) means corrects for inhomogeneous
%     illumination across the vessel width, which the global mean cannot.
%   - Fine angle search resolution increased from 0.25 deg to 0.1 deg
%     and search range widened from +/-2 deg to +/-3 deg.  This reduces
%     velocity error at high flow rates where cot(theta) is very sensitive
%     to small angle changes.
%   - Degenerate window check: if a window has near-zero variance (e.g.
%     saturated frame, no signal), thetas(k) is set to NaN rather than
%     returning a spurious angle.
%   - spread_matrix is now returned as a double array (was implicitly so).
%
% INPUTS
%   data       - matrix of size (nlines x npoints), time in 1st dimension
%   windowsize - number of lines per window, must be a multiple of 4
%
% OUTPUTS
%   thetas        - time-varying angle of RBC streaks (degrees), nsteps x 1
%   the_t         - time-point of each estimate (in lines), nsteps x 1
%   spread_matrix - variance of Radon projections vs angle, nsteps x 180

stepsize = round(0.25 * windowsize);
nlines   = size(data, 1);
npoints  = size(data, 2);
nsteps   = floor(nlines / stepsize) - 3;

angles      = 0:179;
angles_fine = -3 : 0.1 : 3; % finer resolution, wider range than original

spread_matrix      = zeros(nsteps, length(angles));
spread_matrix_fine = zeros(nsteps, length(angles_fine));
thetas             = NaN(nsteps, 1);
the_t              = NaN(nsteps, 1);

for k = 1:nsteps

    the_t(k) = 1 + (k - 1) * stepsize + windowsize / 2;

    row_start = 1 + (k - 1) * stepsize;
    row_end   = (k - 1) * stepsize + windowsize;
    data_hold = data(row_start:row_end, :);

    % Degenerate window check: skip windows with no usable signal
    if var(data_hold(:)) < eps
        continue; % thetas(k) stays NaN
    end

    % Two-way mean subtraction for illumination correction.
    % col_means (1 x npoints) captures the steady spatial intensity
    % gradient across the vessel width.
    % row_means (windowsize x 1) captures the temporal baseline.
    % Adding the global mean back preserves the two-way ANOVA identity:
    %   corrected = data - col_means - row_means + global_mean
    col_means   = mean(data_hold, 1);               % 1 x npoints
    row_means   = mean(data_hold, 2);               % windowsize x 1
    global_mean = mean(data_hold(:));
    data_hold   = data_hold ...
                  - repmat(col_means, windowsize, 1) ...
                  - repmat(row_means, 1, npoints) ...
                  + global_mean;

    % Coarse Radon transform over all 180 integer degrees
    radon_hold         = radon(data_hold, angles);
    spread_matrix(k,:) = var(radon_hold);
    [~, the_theta]     = max(spread_matrix(k, :));
    coarse_angle       = angles(the_theta);

    % Fine Radon transform around the coarse estimate
    radon_hold_fine         = radon(data_hold, coarse_angle + angles_fine);
    spread_matrix_fine(k,:) = var(radon_hold_fine);
    [~, the_theta_fine]     = max(spread_matrix_fine(k, :));
    thetas(k) = coarse_angle + angles_fine(the_theta_fine);

end

% Rotate from Radon convention (0-179, max variance perpendicular to
% streaks) to physical angle measured from vertical (same convention as
% original GetVelocityRadon).
thetas = -1 * (thetas - 90);

end
