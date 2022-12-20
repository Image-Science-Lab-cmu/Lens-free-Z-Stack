clear all; close all;
addpath('funcs');
addpath('util');
addpath('util/vol3d');
warning('off', 'all');

%% Setup
data_dir = 'data';
exp_name = 'points_on_a_line';

% masks to use
pinhole_mask = struct('use_mask', 1, 'name', 'pinhole', ...
                    'mcode', [0 0 0; 0 1 0; 0 0 0], ...
                    'mask_pitch', 1e-1);
stereo_arr = zeros(3, round(1.89 / 1e-1)+2);
stereo_arr(2, 2) = 1; stereo_arr(2, end-1) = 1;
stereo_mask = struct('use_mask', 2, 'name', 'stereo', ...
                    'mcode', stereo_arr, ...
                    'mask_pitch', 1e-1);
% quad_arr = zeros(round(1.89 / 1e-1)+2); 
% quad_arr(2, 2) = 1; quad_arr(end-1, 2) = 1; quad_arr(2, end-1) = 1; quad_arr(end-1, end-1) = 1;
% quad_mask = struct('use_mask', 3, 'name', 'quad', ...
%                    'mcode', quad_arr, ...
%                    'mask_pitch', 1e-1);
mseq_24_100 = struct('use_mask', 4, 'name', 'mseq 24 100', ...
                    'mcode', mseq(2, 4) * mseq(2, 4)', ...
                    'mask_pitch', 100 * 1e-3);
mseq_24_100p = struct('use_mask', 5, 'name', 'mseq 24 100p', ...
                    'mcode', max(0, mseq(2, 4) * mseq(2, 4)'), ...
                    'mask_pitch', 100 * 1e-3);
mseq_28_30 = struct('use_mask', 6, 'name', 'mseq 28 30', ...
                    'mcode', mseq(2, 8) * mseq(2, 8)', ...
                    'mask_pitch', 30 * 1e-3);
mseq_28_30p = struct('use_mask', 7, 'name', 'mseq 28 30p', ...
                    'mcode', max(0, mseq(2, 8) * mseq(2, 8)'), ...
                    'mask_pitch', 30 * 1e-3);
% mseq_27_100 = struct('use_mask', 8, 'name', 'mseq 27 100', ...
%                      'mcode', mseq(2, 7) * mseq(2, 7)', ...
%                      'mask_pitch', 100 * 1e-3);
% mseq_27_100p = struct('use_mask', 9, 'name', 'mseq 27 100p', ...
%                      'mcode', max(0, mseq(2, 7) * mseq(2, 7)'), ...
%                      'mask_pitch', 100 * 1e-3);

% reconstruction algorithms & parameters
cgs_high_reg = struct('name', 'cgs 1e-1', 'recons_type', 'cgs-fft', 'lambda', 1e-1, 'tol', 1e-6, 'maxit', 50);
cgs_med_reg = struct('name', 'cgs 1e-2', 'recons_type', 'cgs-fft', 'lambda', 1e-2, 'tol', 1e-6, 'maxit', 50);
fista_med_reg = struct('name', 'fista 1e-2', 'recons_type', 'fista', 'lambda', 1e-2, 'maxit', 50);
fista_low_reg = struct('name', 'fista 1e-3', 'recons_type', 'fista', 'lambda', 1e-3, 'maxit', 50);

% filter on measurements
no_hpf = struct('name', '', 'type', 'no_filter');

% experiments to run
exps = {struct('mask', pinhole_mask, 'recons_param', cgs_high_reg, 'filter', no_hpf),... 
        struct('mask', pinhole_mask, 'recons_param', fista_med_reg, 'filter', no_hpf), ...
        struct('mask', stereo_mask,  'recons_param', cgs_high_reg,  'filter', no_hpf), ...
        struct('mask', stereo_mask,  'recons_param', fista_low_reg, 'filter', no_hpf), ...
        struct('mask', mseq_24_100p, 'recons_param', cgs_high_reg,  'filter', no_hpf), ...
        struct('mask', mseq_24_100p, 'recons_param', fista_low_reg, 'filter', no_hpf), ...
        struct('mask', mseq_28_30p,  'recons_param', cgs_med_reg,   'filter', no_hpf), ...
        struct('mask', mseq_28_30p, 'recons_param', fista_low_reg, 'filter', no_hpf)};

% use gpu
use_gpu = 1;
for i = 1:length(exps), exps{i}.use_gpu = use_gpu; end

% misc setup
use_id = 64; %use last frame to be consistant in angle coordinates
nsamples = 20;
bin_size = 1;
sensor_pitch = 2.4 * 1e-3 * 2 * 16 * bin_size;
xyz_range = [-5 5 -5 5 5 10];
show_recons_rep = @(vol)plot_volume_lct_style(permute(mat2gray(vol, [0, max(vol(:))]), [1 3 2]), 'angle_diopter');  
should_filter = @(mask_name)strcmp(mask_name(1:4), 'mseq') && strcmp(mask_name(end), 'p');

%% Reconstruct volume from single measurement
folder = sprintf('%s/%s/', data_dir, exp_name);
filename = sprintf('%s/zcam_intensity_0.mat', folder);
load(filename, 'gt_intensity', 'ds', 'gt_v', 'gt_f', 'gt_a');
M = (size(gt_intensity, 1) / length(ds)); % number of masks rendered
C = 1;
N = (size(gt_intensity, 1)/M); % number of frames
Hs = round(size(gt_intensity, 2)/bin_size);
Ws = round(size(gt_intensity, 3)/bin_size);
siz_i = [N Hs Ws C];

zs = -ds(end:-1:1); % scene z values

figure('units', 'normalized', 'outerposition', [0 0 0.5 1])
for jj = 1:length(exps),
    mask_name = exps{jj}.mask.name;
    use_mask = exps{jj}.mask.use_mask;
    mcode = exps{jj}.mask.mcode;
    mask_pitch = exps{jj}.mask.mask_pitch;

    meas_z = zeros(1, Hs, Ws, C);
    meas_z(1, :, :, :) = imresize(squeeze(gt_intensity((use_id-1)*M+use_mask, :, :, :)), [Hs Ws]);

    % reparameterize measurement
    [i_t, p_t, q_t] = i2i_t(meas_z, [ds(use_id)], sensor_pitch);

    % get 3D kernel (write or reads from cache)
    [kernel, k_range, tt_y, tt_x, tt_z] = zcam_kernel_extended(...
        mcode, mask_pitch, p_t, q_t, [ds(use_id)], zs, nsamples, ...
        sprintf('cache/kernel_%s_id%d.mat', mask_name, use_id));

    % plot kernel
    if mod(jj, 2) == 1,
        subplot(4, 3, (ceil(jj/2)-1)*3+1), 
        plot_volume_lct_style(permute(kernel, [1 3 2]), ... %(z, x, y)
                'angle_diopter', ...
                cellfun(@(x)sprintf('%.2f', x), {k_range(1) k_range(2)}, 'UniformOutput', false), ... % xlabels
                cellfun(@(x)sprintf('%.2f', x), {k_range(3) k_range(4)}, 'UniformOutput', false), ... % ylabels
                cellfun(@(x)sprintf('%.2f', x), {k_range(5) k_range(6)}, 'UniformOutput', false) ... % zlabels
                );
         title(sprintf('%s kernel', mask_name));
    end
     
    % reconstruct volume
    rp = exps{jj}.recons_param;
    hp = exps{jj}.filter;
    [i_tf, kernel_f] = hpf_meas_and_kernel(i_t, kernel, should_filter(mask_name), hp);
    rp.lambda = rp.lambda * norm(i_tf(:), 'fro');
    [tt_hat, recon_time] = deconv_3d_extended(i_tf, kernel_f, rp);
    
    % plot reconstruction
    subplot(4, 3, (ceil(jj/2)-1)*3+2+mod(jj+1, 2)),
    show_recons_rep(tt_hat);
    title(rp.name);
    
    fprintf('Done reconstruction for %s mask %s recons parameters!', ...
        mask_name, rp.name);
end

fprintf('Done reconstruction for %d experiments!\n',...
    length(exps));