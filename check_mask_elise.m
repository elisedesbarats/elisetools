function check_mask_elise(data_file, mask_file, varargin)
% First version 05/05/2026 by Elise Desbarats at McGill university using Claude. 
% General use function for overlaying and comparing a mask with your data by performing a % overlap check and plotting the difference.
% Operates purely in the data's native voxel space via spm_vol/spm_read_vol, so will also check if the data is bigger than the mask. 
% Probably most useful for checking wholebrain masks, but can be used to overlay a grey matter mask to see where it is on the brain 
% and check the dimensions.
% 
% Mask is resliced to data space using voxel-to-world interpolation. 
%
% By default, will also peform a check of the headers of the two images to test voxel dimensions. This could also be used
% to compare the headers of two sets of data by passing in a second data file as the mask (but the plot will look odd and 
% the coverage will not be correct so you should pass 'headeronly').
%
%
% See notes for more information
%
% Usage:
%   check_mask_coverage(data_file, mask_file) 
%   check_mask_coverage(data_file, mask_file, 'noheadercheck') % plots only
%   check_mask_coverage(data_file, mask_file, 'noorthviews') % only fmri_display figure
%   check_mask_coverage(data_file1, data_file2, 'headeronly') % compare two data files, dont plot
%
%
% Inputs:
%   data_file      - path to a 3D NIfTI beta/contrast image
%   mask_file      - path to a NIfTI mask image (will be binarized)
%
% Optional:
%   'noheadercheck' - skip the header diagnostic (dimensions, voxel size,
%                     origin, orientation, extra-brain signal, and mask
%                     extent checks). Useful if you have already run the
%                     diagnostic or want a faster visual check only.
%
%
%   'noorthviews'   - skip SPM orthviews (canlab flat display still shown)
%
%   'headeronly'    - only perform the check for dimensions, size etc. 
%                     useful for using this function to compare two data 
%                     files because in that case the other diagnostics will be off
%
%
% Outputs:
%   - printed coverage report (data voxels, mask voxels, covered, missed)
%   - printed header diagnostic (unless 'noheadercheck' is passed)
%   - SPM orthviews: beta image underlay with green (covered) and
%                    red (missed) overlays (unless 'noorthviews' is passed)
%   - canlab_results_fmridisplay flat slice panel that displays the coverage % 
%      result and the paths to the files you inputted.
%     (grey), covered voxels (green), and missed voxels (red) (unless 'headeronly' is passed)
%
% Requirements: SPM on path, CanlabCore tools 

% Example:
%   check_mask_coverage('beta_001.nii', 'gray_matter_mask.nii')
%   check_mask_coverage('beta_001.nii', 'gray_matter_mask.nii', 'noheadercheck')
%
% Notes:
%   - Handles both NaN-coded and zero-coded out-of-brain voxels
%
%   - Coverage is defined as the fraction of data voxels covered by the
%     mask. A data voxel is any finite, nonzero voxel in the beta image
%     regardless of whether it falls within any brain boundary.
%
%   - Header diagnostic (unless 'noheadercheck') runs the following checks
%       ONLY on the file passed as data:
%       1. Dimensions      : voxel grid size must match after reslicing
%       2. Voxel size      : flags mismatches > 0.01 mm
%       3. Origin          : flags shifts > 0.5 mm in any direction
%       4. Orientation     : flags rotation/flip differences > 0.001
%       5. Extra-brain signal (Check 2): reports what fraction of data
%          voxels fall outside the CANlab default brainmask. A high fraction
%          suggests your beta image has substantial non-brain signal which
%          will inflate n_data and make coverage appear worse than it is.
%       6. Mask extent (Check 4): compares the spatial extent of your
%          candidate mask against the CANlab brainmask in mm along each
%          axis (X: L-R, Y: A-P, Z: I-S). Flags axes where the candidate
%          mask does not reach the boundary of the CANlab brain by more
%          than 4mm (~1-2 voxels). This can reveal if your mask was
%          defined on a different MNI152 template and is systematically
%          undersized in one direction (e.g. missing inferior cerebellum).
%
%   - Both checks 2 and 4 use brainmask_canlab.nii as the MNI152
%     reference. If not found via which('brainmask_canlab.nii'), the
%     function attempts to locate it relative to the CanlabCore root.
%     Ensure CanlabCore appears before SPM toolboxes on the MATLAB path
%     to avoid picking up SPM's own brainmask.nii instead.
%
%   - Reslicing uses nearest-neighbour interpolation (appropriate for
%     binary masks). If your mask is probabilistic, it will be binarized
%     at any value > 0 after reslicing.
%
%
% I created this using Claude as a way to sanity check wholebrain masks, and haven't extensively tested other use cases. 
% Feel free to make further edits to fit with your analysis needs. 


%% -----------------------------------------------------------------------
%  1. READ DATA AND MASK WITH SPM (no CANlab implicit masking) + check varargin
%
% ------------------------------------------------------------------------
% parse optional arguments
run_headercheck = true;  % default on
run_orthviews   = true;

if any(strcmpi(varargin, 'noheadercheck')),  run_headercheck = false; end
if any(strcmpi(varargin, 'noorthviews')),    run_orthviews   = false; end

fprintf('\nReading data:  %s\n', data_file);
V_dat  = spm_vol(data_file);
Y_dat  = spm_read_vols(V_dat);          % native voxel space, no masking

fprintf('Reading mask:  %s\n', mask_file);
V_mask = spm_vol(mask_file);




%% -----------------------------------------------------------------------
%  2. RESLICE MASK INTO DATA SPACE IF NEEDED
% ------------------------------------------------------------------------
same_space = isequal(V_dat.dim(1:3), V_mask.dim(1:3)) && ...
             isequal(V_dat.mat,      V_mask.mat);

if same_space
    fprintf('Mask already in data space — no reslicing needed.\n');
    Y_mask = spm_read_vols(V_mask);
else
    fprintf('Reslicing mask to data space via voxel-to-world interpolation...\n');

    % Build a grid of voxel coordinates in data space
    [X, Y, Z]  = ndgrid(1:V_dat.dim(1), 1:V_dat.dim(2), 1:V_dat.dim(3));
    xyz_vox    = [X(:), Y(:), Z(:), ones(numel(X), 1)]';   % 4 x N

    % Transform: data voxel coords -> world -> mask voxel coords
    xyz_world      = V_dat.mat  * xyz_vox;          % world coords
    xyz_mask_vox   = V_mask.mat \ xyz_world;        % mask voxel coords

    % Sample mask at those coordinates (nearest neighbour — binary mask)
    Y_mask_flat = spm_sample_vol(V_mask, ...
        xyz_mask_vox(1,:), ...
        xyz_mask_vox(2,:), ...
        xyz_mask_vox(3,:), ...
        0);                     % 0 = nearest neighbour

    Y_mask = reshape(Y_mask_flat, V_dat.dim(1:3));
end

%% -----------------------------------------------------------------------
%  3. BUILD BINARY MASKS
%     Handle both NaN-coded and 0-coded out-of-brain voxels
%       EXTRA: gives dignostics about the header 
% ------------------------------------------------------------------------
% Data: voxel is "present" if it is finite and nonzero
d_binary = isfinite(Y_dat) & (Y_dat ~= 0);   % logical 3D volume

% Mask: binarize at any positive value
m_binary = isfinite(Y_mask) & (Y_mask > 0);  % logical 3D volume

if run_headercheck
    fprintf('Running header diagnostics');
    diagnose_mask_header(V_dat, V_mask,  d_binary); % see function at the bottom of the file.
end 

%% -----------------------------------------------------------------------
%  4. COVERAGE STATISTICS
%     % covered = fraction of DATA voxels that the mask covers
% ------------------------------------------------------------------------
n_data    = sum(d_binary(:));
n_mask    = sum(m_binary(:));
n_covered = sum(d_binary(:) & m_binary(:));
n_missed  = sum(d_binary(:) & ~m_binary(:));   % data voxels outside mask
n_extra   = sum(m_binary(:) & ~d_binary(:));   % mask voxels with no data

pct_covered = 100 * n_covered / n_data;
pct_missed  = 100 * n_missed  / n_data;

fprintf('\n========================================\n');
fprintf('  MASK COVERAGE REPORT\n');
fprintf('========================================\n');
fprintf('  Data voxels              : %6d\n', n_data);
fprintf('  Mask voxels              : %6d\n', n_mask);
fprintf('  Covered (data & mask)    : %6d  (%5.1f%% of data)\n', n_covered, pct_covered);
fprintf('  Missed  (data, no mask)  : %6d  (%5.1f%% of data)\n', n_missed,  pct_missed);
fprintf('  Extra   (mask, no data)  : %6d  (%5.1f%% of mask)\n', n_extra,   100*n_extra/n_mask);
fprintf('========================================\n\n');

if n_missed == 0
    fprintf('[OK]  Mask fully covers all data voxels.\n\n');
else
    fprintf('[!!]  Mask is MISSING %d data voxels (%.1f%%).\n\n', n_missed, pct_missed);
end

%% -----------------------------------------------------------------------
%  5. WRITE TEMPORARY NIFTIS
%     data outline, covered, missed — all in data voxel space
% ------------------------------------------------------------------------
    function write_tmp(fname, vol_binary, V_ref)
        V_out       = V_ref;
        V_out.fname = fname;
        V_out.dt    = [spm_type('uint8') 0];
        V_out.pinfo = [1 0 0]';
        spm_write_vol(V_out, double(vol_binary));
    end

tmp_data    = fullfile(pwd, 'tmp_data_outline.nii');
tmp_covered = fullfile(pwd, 'tmp_covered.nii');
tmp_missed  = fullfile(pwd, 'tmp_missed.nii');

write_tmp(tmp_data,    d_binary,             V_dat);
write_tmp(tmp_covered, d_binary & m_binary,  V_dat);
write_tmp(tmp_missed,  d_binary & ~m_binary, V_dat);

%% -----------------------------------------------------------------------
%  6. SPM ORTHVIEWS - skip this by passing 'noorthviews' as varargin
%     Underlay : raw beta image
%     Green    : covered voxels
%     Red      : missed voxels
% ------------------------------------------------------------------------
if run_orthviews
    spm_figure('GetWin', 'Graphics');
    spm_orthviews('Reset');
    h = spm_orthviews('Image', data_file, [0.0 0.5 1.0 0.5]);
    
    spm_orthviews('AddColouredImage', h, tmp_covered, [0 0.8 0]);
    
    if n_missed > 0
        spm_orthviews('AddColouredImage', h, tmp_missed, [1 0 0]);
    end
    
    spm_orthviews('Redraw');
end 

%% -----------------------------------------------------------------------
%  7. CANLAB RESULTS DISPLAY
%     Pass d_binary as explicit mask so CANlab doesn't impose its brainmask.
%     Load covered/missed with that same mask so voxel spaces are consistent.
% ------------------------------------------------------------------------

% Load all three images using the data outline as the explicit mask
% This forces CANlab to operate in your data's voxel space, not its own
data_obj    = fmri_data(tmp_data,    tmp_data);   % mask = itself = full data space
covered_obj = fmri_data(tmp_covered, tmp_data);   % masked to data space
missed_obj  = fmri_data(tmp_missed,  tmp_data);

% Display data outline as grey underlay in the flat slice panel
o2 = canlab_results_fmridisplay(data_obj, 'compact2', 'noverbose', ...
    'pos', [0.0 0.0 1.0 0.5]);

% Overlay covered in green
o2 = addblobs(o2, region(covered_obj), 'color', [0 0.8 0], 'alpha', 0.6);

% Overlay missed in red
if n_missed > 0
    o2 = addblobs(o2, region(missed_obj), 'color', [1 0 0], 'alpha', 0.9);
end

%% -----------------------------------------------------------------------
%  8. TITLES AND ANNOTATION
% ------------------------------------------------------------------------
sgtitle(sprintf('Mask coverage | Covered: %.1f%% | Missed: %d voxels', ...
    pct_covered, n_missed));

annotation('textbox', [0.65 0.75 0.3 0.2], ...
    'String', sprintf('Data:\n%s\n\nMask:\n%s', data_file, mask_file), ...
    'FitBoxToText', 'on', ...
    'Interpreter', 'none', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 12);

%% -----------------------------------------------------------------------
%  9. CLEAN UP
% ------------------------------------------------------------------------
delete(tmp_data);
delete(tmp_covered);
delete(tmp_missed);

%% -----------------------------------------------------------------------
%  EXTRA. HEADER DIAGNOSTIC SUBFUNCTION
% ------------------------------------------------------------------------

function diagnose_mask_header(V_dat, V_mask, d_binary)
issues = {}; 

vox_dat  = sqrt(sum(V_dat.mat(1:3,1:3).^2));
vox_mask = sqrt(sum(V_mask.mat(1:3,1:3).^2));

origin_dat  = -V_dat.mat(1:3,1:3)  \ V_dat.mat(1:3,4);
origin_mask = -V_mask.mat(1:3,1:3) \ V_mask.mat(1:3,4);

R_dat  = V_dat.mat(1:3,1:3)  ./ vox_dat;
R_mask = V_mask.mat(1:3,1:3) ./ vox_mask;

tol_vox    = 0.01;
tol_origin = 0.5;
tol_rot    = 0.001;

fprintf('\n========================================\n');
fprintf('  MASK HEADER DIAGNOSTIC\n');
fprintf('========================================\n');

% Dimensions
if ~isequal(V_dat.dim(1:3), V_mask.dim(1:3))
    fprintf('  [!!] DIMENSION MISMATCH: data [%d %d %d] vs mask [%d %d %d]\n', ...
        V_dat.dim(1:3), V_mask.dim(1:3));
else
    fprintf('  [OK] Dimensions match: [%d %d %d]\n', V_dat.dim(1:3));
end

% Voxel size
if any(abs(vox_dat - vox_mask) > tol_vox)
    fprintf('  [!!] VOXEL SIZE MISMATCH: data [%.3f %.3f %.3f] vs mask [%.3f %.3f %.3f] mm\n', ...
        vox_dat, vox_mask);
else
    fprintf('  [OK] Voxel size: [%.3f %.3f %.3f] mm\n', vox_dat);
end

% Origin
origin_diff_mm = abs((origin_dat - origin_mask) .* vox_dat);
if any(origin_diff_mm > tol_origin)
    fprintf('  [!!] ORIGIN MISMATCH: difference [%.3f %.3f %.3f] mm\n', origin_diff_mm);
else
    fprintf('  [OK] Origin match within tolerance\n');
end

% Orientation
rot_diff = abs(R_dat - R_mask);
if any(rot_diff(:) > tol_rot)
    fprintf('  [!!] ORIENTATION MISMATCH: max element difference %.6f\n', max(rot_diff(:)));
else
    fprintf('  [OK] Orientation matches\n');
end


% Load CANlab brainmask as reference
canlab_brain = which('brainmask_canlab.nii');
if isempty(canlab_brain)
    fprintf('  [!!] CANlab brainmask not found on path — check CanlabCore on path?');
else
    V_brain   = spm_vol(canlab_brain);
    Y_brain   = spm_read_vols(V_brain);
    
    % Reslice CANlab mask to data space using same interpolation as main function
    [X, Y, Z]      = ndgrid(1:V_dat.dim(1), 1:V_dat.dim(2), 1:V_dat.dim(3));
    xyz_vox        = [X(:), Y(:), Z(:), ones(numel(X), 1)]';
    xyz_world      = V_dat.mat * xyz_vox;
    xyz_brain_vox  = V_brain.mat \ xyz_world;
    brain_flat     = spm_sample_vol(V_brain, ...
        xyz_brain_vox(1,:), xyz_brain_vox(2,:), xyz_brain_vox(3,:), 0);
    brain_binary   = reshape(brain_flat, V_dat.dim(1:3)) > 0;

    % Check 2: data voxels outside CANlab brain
    n_data_outside_brain = sum(d_binary(:) & ~brain_binary(:));
    pct_outside          = 100 * n_data_outside_brain / sum(d_binary(:));
    fprintf('\n[ Check 2: Data voxels outside CANlab brainmask ]\n');
    fprintf('  Note: CANlab brainmask used as brain boundary reference (%s)\n', canlab_brain);
    fprintf('  Data voxels outside brain : %d (%.1f%% of data)\n', n_data_outside_brain, pct_outside);
    if pct_outside > 5
        fprintf('  [!!] SUBSTANTIAL EXTRA-BRAIN SIGNAL — this may inflate n_data and\n');
        fprintf('       make coverage look worse than it is. Consider masking your\n');
        fprintf('       data before running check_mask_coverage.\n');
    elseif pct_outside > 0
        fprintf('  [~]  Minor extra-brain signal (%.1f%%) — unlikely to affect results.\n', pct_outside);
    else
        fprintf('  [OK] All data voxels fall within CANlab brainmask.\n');
    end

    % Check 4: mask extent vs CANlab brain extent
    % Reslice candidate mask to CANlab brain space for a fair extent comparison
    [Xb, Yb, Zb]    = ndgrid(1:V_brain.dim(1), 1:V_brain.dim(2), 1:V_brain.dim(3));
    xyz_vox_b        = [Xb(:), Yb(:), Zb(:), ones(numel(Xb), 1)]';
    xyz_world_b      = V_brain.mat * xyz_vox_b;
    xyz_mask_vox_b   = V_mask.mat \ xyz_world_b;
    mask_in_brain_sp = spm_sample_vol(V_mask, ...
        xyz_mask_vox_b(1,:), xyz_mask_vox_b(2,:), xyz_mask_vox_b(3,:), 0);
    mask_brain_sp    = reshape(mask_in_brain_sp, V_brain.dim(1:3)) > 0;

    % Compare extents in mm along each axis
    brain_world = V_brain.mat * xyz_vox_b;
    mask_world  = V_brain.mat * xyz_vox_b;   % same grid, different occupancy

    labels = {'X (L-R)', 'Y (A-P)', 'Z (I-S)'};
    fprintf('\n[ Check 4: Mask extent vs CANlab brain extent (mm) ]\n');
    fprintf('  Note: using CANlab brainmask as MNI152 reference\n');

    extent_issues = {};
    for ax = 1:3
        brain_coords = brain_world(ax, Y_brain(:) > 0);
        mask_coords  = brain_world(ax, mask_brain_sp(:));

        brain_range  = [min(brain_coords) max(brain_coords)];
        mask_range   = [min(mask_coords)  max(mask_coords)];

        undershoot_lo = brain_range(1) - mask_range(1);  % +ve = mask doesn't extend as far
        undershoot_hi = mask_range(2)  - brain_range(2); % +ve = mask extends further

        fprintf('  %s:\n', labels{ax});
        fprintf('    Brain : [%.1f  %.1f] mm\n', brain_range(1),  brain_range(2));
        fprintf('    Mask  : [%.1f  %.1f] mm\n', mask_range(1),   mask_range(2));

        tol_extent = 4;  % mm — one voxel at 4mm, two at 2mm
        if undershoot_lo > tol_extent
            fprintf('    [!!] Mask does not extend to inferior/posterior/left boundary by %.1f mm\n', undershoot_lo);
            extent_issues{end+1} = sprintf('%s low end', labels{ax});
        end
        if -undershoot_hi > tol_extent
            fprintf('    [!!] Mask does not extend to superior/anterior/right boundary by %.1f mm\n', -undershoot_hi);
            extent_issues{end+1} = sprintf('%s high end', labels{ax});
        end
        if isempty(extent_issues)
            fprintf('    [OK]\n');
        end
    end

    if ~isempty(extent_issues)
        issues{end+1} = 'mask extent';
    end
end


% Verdict
issues = {};
if ~isequal(V_dat.dim(1:3), V_mask.dim(1:3)),  issues{end+1} = 'dimensions';   end
if any(abs(vox_dat - vox_mask) > tol_vox),      issues{end+1} = 'voxel size';   end
if any(origin_diff_mm > tol_origin),            issues{end+1} = 'origin';       end
if any(rot_diff(:) > tol_rot),                  issues{end+1} = 'orientation';  end

fprintf('----------------------------------------\n');
if isempty(issues)
    fprintf('  [OK] No header issues — mask and data are in the same space.\n');
else
    fprintf('  [!!] Issues: %s\n', strjoin(issues, ', '));
    fprintf('       Reslicing will correct for this but consider resampling offline if severe.\n');
end
fprintf('========================================\n\n');

end




end
