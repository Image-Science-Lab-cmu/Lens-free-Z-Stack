function [i_t, p_t, q_t] = i2i_t(meas, ds, pitch, p_t, q_t)
% i_t(p_t, q_t, d_t) = 1/d_t^2 i(p_t/d_t, q_t/d_t)
    % meas is (N, H, W, C)
    N = size(meas, 1); assert(length(ds) == N);
    H = size(meas, 2); 
    W = size(meas, 3); 
    C = size(meas, 4);
    
    % choose cone given by furthest d
    if nargin < 4,
        [val, id_max] = max(abs(ds));
        d_max = ds(id_max);
        p0s = linspace(-0.5 * W * pitch, 0.5 * W * pitch, W);
        q0s = linspace(-0.5 * H * pitch, 0.5 * H * pitch, H);
        p_t = sort(p0s/d_max, 'ascend');
        q_t = sort(q0s/d_max, 'ascend'); % make p_t and q_t ascending
    end

    [tmp_x, tmp_y] = meshgrid(p_t, q_t);
    i_t = zeros(N, H, W, C);
    for i = 1:N, for c = 1:C,
        dt = 1/ds(i);
        i_t(i, :, :, c) = dt.^(-2) * reshape(trianglefilt_and_interp2(...
            p0s, q0s, squeeze(meas(i, :, :, c)), tmp_x/dt, tmp_y/dt, ...
            pitch_of(p0s), pitch_of(q0s), pitch_of(p_t/dt), pitch_of(q_t/dt)),...
            1, length(q_t), length(p_t), 1);
    end; end
end