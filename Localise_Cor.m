clearvars; clc;

%% Paths and pixel sizes
avg_psf_path     = 'C:\Users\zhaor\Desktop\PSF_avg2.tif';
sample_psf_paths = { ...
    'C:\Users\zhaor\Desktop\PSF_3.tif'
};

px_x = 0.0451379;          % μm/pixel (lateral)
px_z = 0.0799050;          % μm/slice (axial)

%% Load data
avg_all = double(tiffreadVolume(avg_psf_path));
[H, W, Zall] = size(avg_all);

N = numel(sample_psf_paths);
samp_all = cell(1, N);
for n = 1:N
    S = double(tiffreadVolume(sample_psf_paths{n}));
    [Hs, Ws, Zs] = size(S);
    assert(Hs == H && Ws == W, 'XY dimensions mismatch (sample %d)', n);

    if Zs == Zall
        samp_all{n} = S;

    elseif Zs > Zall
        % — Sample has more z-slices than avg: center-crop around peak to match Zall
        zprof = squeeze(S(round((H+1)/2), round((W+1)/2), :));
        [~, id_s] = max(zprof);
        half = floor(Zall/2);
        z1 = max(1, id_s - half);
        z2 = min(Zs, z1 + Zall - 1);
        z1 = max(1, z2 - Zall + 1);
        samp_all{n} = S(:,:, z1:z2);

    else
        % — Sample has fewer z-slices: symmetrically pad along z without interpolation
        pad_total = Zall - Zs;
        pad_pre  = floor(pad_total/2);
        pad_post = pad_total - pad_pre;
        samp_all{n} = cat(3, repmat(S(:,:,1),1,1,pad_pre), S, repmat(S(:,:,end),1,1,pad_post));
    end
end

%% — Crop to Zkeep slices centered at each peak (+ optional manual offset) —
Zkeep = 30;                 % Increase for a wider window (e.g., 60)
offset_avg  = 0;            % <— Shift avg PSF window by this many slices (negative: toward smaller z)
offset_samp = 0;            % <— Shift sample window (usually 0)

% 1) AVG
[~, id_avg]  = max(squeeze(avg_all(round((H+1)/2), round((W+1)/2), :)));
z_idx_avg    = center_to_idx(id_avg,  Zall, Zkeep, offset_avg);
avg_psf      = avg_all(:,:,z_idx_avg);

% 2) SAMPLES
sample_stacks = cell(1, N);
for n = 1:N
    zprof = squeeze(samp_all{n}(round((H+1)/2), round((W+1)/2), :));
    [~, id_s]   = max(zprof);
    z_idx_samp  = center_to_idx(id_s,  Zall, Zkeep, offset_samp);
    sample_stacks{n} = samp_all{n}(:,:,z_idx_samp);

    fprintf('Sample %d: avg [%d..%d], samp [%d..%d]\n', ...
        n, z_idx_avg(1), z_idx_avg(end), z_idx_samp(1), z_idx_samp(end));
end

%% Parameter settings (stable + efficient)
opts.xy_method = 'gauss';

% =========================== ROI selection ===========================
% Two modes: fixed center vs centroid-based (recenter)
%
% - 'fixed'    : ROI fixed at geometric image center; no per-layer XY alignment
%                 → preserves real lateral (XY) motion.
%                 e.g., center = [(W+1)/2, (H+1)/2]
%
% - 'recenter' : ROI follows brightest point or centroid; each layer is shifted
%                 to align center → produces cleaner axial matching but removes true XY drift.
%
% Guidelines:
%   1) To preserve true XY motion → use 'fixed' and set do_internal_xy_align = false.
%   2) For cleaner axial matching → use 'recenter' and set do_internal_xy_align = true.
%   3) For reporting trajectories → estimate XY with 'fixed', combine Z from 'recenter'.
%
opts.roi_center = [];              % Automatically use 3D intensity peak
opts.do_internal_xy_align = true;

opts.roi_halfsize = 16;
opts.min_corr = 0.5;
opts.sigma = 8;
opts.z_upsample = 1;               % No slice interpolation
opts.z_window   = 14;              % Half-width of sliding window
opts.mask_radius_ratio = 0.35;
opts.use_zscore = false;

opts.enforce_monotonic = true;     % Enforce monotonic z-mapping
opts.compute_full_corr = true;     % Compute full Z×Z correlation matrix for diagnostics

assert( 2*opts.z_window + 1 <= size(avg_psf,3), ...
    'z_window too large: must satisfy 2*z_window+1 <= Zkeep (current Zkeep=%d, z_window=%d)', ...
    size(avg_psf,3), opts.z_window);

%% Processing
spList_all = cell(1,N); map_all = cell(1,N);
corr_all   = cell(1,N); info_all= cell(1,N);

for n=1:N
    fprintf('Processing sample %d/%d (centered %d slices)...\n', n, N, Zkeep);
    [spList, map_z, corr_vec, info] = psf_corr_map_to_avg(avg_psf, sample_stacks{n}, px_x, px_z, opts);
    spList_all{n} = spList; map_all{n} = map_z;
    corr_all{n} = corr_vec; info_all{n} = info;
end

%% Plot 3D trajectory (enforced 3D view)
sample_id = 1;
spList = spList_all{sample_id};
xs = []; ys = []; zs = []; cs = [];
for r = 1:numel(spList)
    v = spList{r};
    if ~isempty(v) && numel(v) >= 3 && all(isfinite(v(1:3)))
        xs(end+1) = v(1); ys(end+1) = v(2); zs(end+1) = v(3);
        cs(end+1) = (numel(v) >= 4)*v(4) + (numel(v) < 4)*NaN;
    end
end

% Convert x, y from pixels to μm and re-center (for intuitive visualization)
cx = (W+1)/2; cy = (H+1)/2;
x_um = (xs - cx) * px_x;
y_um = (ys - cy) * px_x;  % If anisotropic pixels, change to px_y

figure('Name','PSF fitted positions (single sample)'); hold on;
if ~isempty(xs)
    plot3(x_um, y_um, zs, 'o-');
    scatter3(x_um, y_um, zs, 22, cs, 'filled');
end
grid on; axis vis3d; view(45,20);
xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('z (\mum)');
title(sprintf('Sample %d fitted locations (centered %d slices)', sample_id, Zkeep));
hold off;

%% Correlation heatmap (full Z×Z)
C = info_all{sample_id}.corr_full;   % Z x Z full alignment map
figure('Name','Correlation (full ZxZ)');
imagesc(C); axis image; colorbar;
xlabel('sample z index'); ylabel('avg z index');
title('Pearson correlation (full search map)');

%% Z-layer mapping curve
map_z = map_all{sample_id};
figure('Name','Layer mapping');
plot(1:numel(map_z), map_z, '-', 'LineWidth', 1.5);
xlabel('sample z index'); ylabel('mapped avg z index');
title('z-layer mapping (sample → avg)'); grid on;

%% Helper function
function idx = center_to_idx(center, Zall, Zkeep, offset)
center = round(center) + offset;
half   = floor(Zkeep/2);
z1     = max(1, center - half);
z2     = min(Zall, z1 + Zkeep - 1);
z1     = max(1, z2 - Zkeep + 1);
idx    = z1:z2;
end
