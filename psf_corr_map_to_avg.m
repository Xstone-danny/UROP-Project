function [sortedProbeList, map_z, corr_vec, info] = psf_corr_map_to_avg(avg_psf, sample_psf, px_x, px_z, opts)

% ===== Default options =====
if nargin < 5, opts = struct; end
def = struct('sigma',6,'z_upsample',1,'use_zscore',false,'mask_radius_ratio',0.35, ...
             'z_window',[],'do_internal_xy_align',true, ...
             'roi_center',[],'roi_halfsize',24,'xy_method','gauss','min_corr',0.45, ...
             'enforce_monotonic',true, 'compute_full_corr', true);
fn = fieldnames(def);
for k = 1:numel(fn)
    f = fn{k}; if ~isfield(opts, f), opts.(f) = def.(f); end
end

% ===== Basic checks / conversions =====
avg_psf    = double(avg_psf);
sample_psf = double(sample_psf);
[H, W, Z] = size(avg_psf);
assert(isequal(size(sample_psf), [H, W, Z]), 'avg_psf and sample_psf must have identical sizes');

% ===== Smoothing =====
if opts.sigma > 0
    for k = 1:Z
        avg_psf(:,:,k)    = imgaussfilt(avg_psf(:,:,k),    opts.sigma, 'padding', 'replicate');
        sample_psf(:,:,k) = imgaussfilt(sample_psf(:,:,k), opts.sigma, 'padding', 'replicate');
    end
end

% ===== z upsampling (off by default) =====
if opts.z_upsample > 1
    [avg_psf, Z2]   = z_interp_linear(avg_psf, opts.z_upsample);
    [sample_psf, ~] = z_interp_linear(sample_psf, opts.z_upsample);
    Z = Z2;  px_z = px_z / opts.z_upsample;
end

% ===== Mask and ROI =====
cx = (W+1)/2; cy = (H+1)/2;
R  = opts.mask_radius_ratio * min(H, W);
mask = make_circular_mask(H, W, cx, cy, R);
mask_lin = find(mask);

if isempty(opts.roi_center)
    % Use the global 3D maximum of the sample as the ROI center
    [~, id3] = max(sample_psf(:));
    [ry0, cx0, ~] = ind2sub([H, W, Z], id3);
    x0 = cx0; y0 = ry0;
else
    x0 = opts.roi_center(1); y0 = opts.roi_center(2);
end
hs = opts.roi_halfsize;
x1 = max(1, round(x0 - hs)); x2 = min(W, round(x0 + hs));
y1 = max(1, round(y0 - hs)); y2 = min(H, round(y0 + hs));
roi_sz = [y2 - y1 + 1, x2 - x1 + 1];

% ===== (x,y) localisation within a fixed ROI (gauss/centroid/peak) =====
xy_rc = nan(Z, 2);
for s = 1:Z
    patch = sample_psf(y1:y2, x1:x2, s);
    switch lower(opts.xy_method)
        case 'peak'
            [~, idx] = max(patch(:));
            [rr, cc] = ind2sub(size(patch), idx);
            xs = cc + x1 - 1; ys = rr + y1 - 1;

        case 'centroid'
            [xx, yy] = meshgrid(x1:x2, y1:y2);
            w = patch; w = w - min(w(:)); w = max(w, 0);
            ssum = sum(w(:)); if ssum == 0, ssum = eps; end
            xs = sum(xx(:).*w(:)) / ssum; ys = sum(yy(:).*w(:)) / ssum;

        otherwise % 'gauss'
            if exist('lateral_localisation','file') == 2
                probe_image = patch(:); probe_list = (1:numel(probe_image)).';
                s_estimate = 0; cutLength = max(roi_sz); mode = 2;
                [xg, yg] = lateral_localisation(probe_list, probe_image, roi_sz, s_estimate, cutLength, mode);
                xs = xg + x1 - 1; ys = yg + y1 - 1;
            else
                % Fallback to centroid if the Gaussian localiser is unavailable
                [xx, yy] = meshgrid(x1:x2, y1:y2);
                w = patch; w = w - min(w(:)); w = max(w, 0);
                ssum = sum(w(:)); if ssum == 0, ssum = eps; end
                xs = sum(xx(:).*w(:)) / ssum; ys = sum(yy(:).*w(:)) / ssum;
            end
    end
    xy_rc(s, :) = [ys, xs];
end

% ===== Internal copy alignment (for correlation only) =====
avg_work    = avg_psf;
sample_work = sample_psf;
if opts.do_internal_xy_align
    for s = 1:Z
        dr = cy - xy_rc(s,1);
        dc = cx - xy_rc(s,2);
        sample_work(:,:,s) = imtranslate(sample_work(:,:,s), [dc, dr], 'linear', 'FillValues', 0);
    end
end

% ===== Axial peak indices + window size =====
zprof_avg  = squeeze(avg_work(round(cy), round(cx), :));
zprof_samp = squeeze(sample_work(round(cy), round(cx), :));
[~, id_avg]  = max(zprof_avg);
[~, id_samp] = max(zprof_samp);

if isempty(opts.z_window)
    zwin = max(15, round(0.2 * Z));
else
    zwin = opts.z_window;
end

% ===== Pack masked pixels into column vectors =====
A = reshape(avg_work,    [], Z); A = A(mask_lin, :);   % [M, Z]
S = reshape(sample_work, [], Z); S = S(mask_lin, :);   % [M, Z]

% ===== Mapping: sliding-window Pearson correlation (fast and stable) =====
map_z    = nan(Z, 1);
corr_vec = nan(Z, 1);
for s = 1:Z
    a0 = id_avg + (s - id_samp);
    a1 = max(1, round(a0 - zwin)); a2 = min(Z, round(a0 + zwin));
    win = a1:a2;

    Aw  = A(:, win); muA = mean(Aw, 1); Awc = Aw - muA; denA = sqrt(sum(Awc.^2, 1));
    Sv  = S(:, s);   muS = mean(Sv);    Svc = Sv - muS;  denS = sqrt(sum(Svc.^2));
    if denS < eps
        rvec = nan(numel(win), 1);
    else
        rvec = (Awc' * Svc) ./ (denA.' * denS);
    end

    % Integer peak (max correlation within the window)
    [corr_peak, imax] = max(rvec);
    map_int = a1 + imax - 1;   % Integer-layer index mapped onto avg z
    % If you only need the best matched z layer, delete everything below down to the monotonicity section.

    % ---- sub-z: 3-point parabolic refinement around imax (-1, 0, +1) ----
    map_sub = map_int;         % Default to the integer solution
    if imax > 1 && imax < numel(rvec)
        y1 = rvec(imax-1);     % left
        y2 = rvec(imax);       % center (peak)
        y3 = rvec(imax+1);     % right
        denom = (y1 - 2*y2 + y3);
        if abs(denom) > eps
            d = 0.5 * (y1 - y3) / denom;   % fractional offset (unit: slice)
            % Optional clamp to [-0.5, 0.5] to avoid outliers
            d = max(min(d, 0.5), -0.5);
            map_sub = map_int + d;         % continuous index (sub-z)
        end
    end

    % Write back: correlation peak and sub-z index
    corr_vec(s) = corr_peak;
    map_z(s)    = map_sub;

    % (Optional) If you also want the integer layer for comparison, preallocate at top:
    %   map_z_int = nan(Z,1);
    % and here add:
    %   map_z_int(s) = map_int;
end

% ===== Monotonicity (optional) =====
if opts.enforce_monotonic
    for s = 2:Z
        if map_z(s) < map_z(s-1), map_z(s) = map_z(s-1); end
    end
end

% ===== Full Z×Z correlation map (for diagnostic visualization only) =====
if opts.compute_full_corr
    % Zero-mean / unit-norm columns; dot products approximate Pearson coefficients
    A0 = A - mean(A, 1); SA = sqrt(sum(A0.^2, 1)); SA(SA < eps) = 1; A0 = A0 ./ SA;
    S0 = S - mean(S, 1); SS = sqrt(sum(S0.^2, 1)); SS(SS < eps) = 1; S0 = S0 ./ SS;
    corr_full = A0' * S0;   % Z x Z
else
    corr_full = [];
end

% ===== Output list (filter by minimum correlation) =====
sortedProbeList = cell(Z, 1);
for s = 1:Z
    a = map_z(s);
    if ~isnan(a) && corr_vec(s) >= opts.min_corr
        z_um = (a - 1) * px_z;
        xs = xy_rc(s,2); ys = xy_rc(s,1);
        sortedProbeList{s,1} = [xs, ys, z_um, corr_vec(s)];
    else
        sortedProbeList{s,1} = [];
    end
end

% ===== Diagnostic info =====
info = struct();
info.mask = mask;
info.corr_full = corr_full;   % Full Z×Z correlation map
info.xy_rc = xy_rc;
info.id_avg = id_avg; info.id_samp = id_samp;
info.roi = struct('x1',x1,'x2',x2,'y1',y1,'y2',y2,'center',[x0,y0],'halfsize',hs);
info.px_x = px_x; info.px_z = px_z;

end

% ==== Helper functions ====
function [V2, Z2] = z_interp_linear(V, up)
[H, W, Z] = size(V);
if up <= 1, V2 = V; Z2 = Z; return; end
Z2 = (Z - 1) * up + 1;
V2 = zeros(H, W, Z2, 'like', V);
V2(:,:,1:up:end) = V;
for k = 1:Z-1
    for t = 1:up-1
        alpha = t / up;
        V2(:,:, (k-1)*up + 1 + t ) = (1 - alpha) * V(:,:,k) + alpha * V(:,:,k+1);
    end
end
end

function M = make_circular_mask(H, W, cx, cy, R)
[xg, yg] = meshgrid(1:W, 1:H);
M = (xg - cx).^2 + (yg - cy).^2 <= R.^2;
end

