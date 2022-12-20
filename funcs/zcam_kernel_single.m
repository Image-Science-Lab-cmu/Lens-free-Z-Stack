function [kernel, kernel_range, tt_y, tt_x, tt_s] = zcam_kernel_single(mcode, mask_pitch, p_t, q_t, ds, use_id, nsamples, matname)
    % load or build zcam kernel
    if nargin < 6,
        nsamples = 5;
        matname = '';
    end
    
    if exist(matname, 'file') % attempt to load kernel
        load(matname, 'kernel', 'kernel_range', 'tt_y', 'tt_x', 'tt_s');
        % verify dimensions match
        if all(size(tt_y) == size(q_t)) && all(size(tt_x) == size(p_t)) && all(size(tt_s) == size(ds)),
            return;
        end
    end
        
    dt = 1 ./ ds; %ds is assume to be linear in 1./ds
    xx_pitch = pitch_of(p_t);
    yy_pitch = pitch_of(q_t);
    ss_pitch = pitch_of(dt);
    tt_y = q_t;
    tt_x = p_t;
    tt_s = -dt(end:-1:1);
    
    mt_y = min(q_t) - max(tt_y) : yy_pitch : max(q_t) - min(tt_y) + 1e-8;
    mt_x = min(p_t) - max(tt_x) : xx_pitch : max(p_t) - min(tt_x) + 1e-8;
    mt_s = dt(use_id) - max(tt_s)  : ss_pitch : dt(use_id) - min(tt_s) + 1e-8;
    
    mask_height = size(mcode, 1) * mask_pitch; 
    mask_width  = size(mcode, 2) * mask_pitch;
    mask_yy = linspace(-(mask_height/2 - mask_pitch/2), ...
                        mask_height/2 - mask_pitch/2, size(mcode, 1));
    mask_xx = linspace(-(mask_width/2 - mask_pitch/2), ...
                        mask_width/2 - mask_pitch/2,  size(mcode, 2));

    % Monte carlo render to avoid aliasing
    % 1. compute low pass filter widths
    xpd = makedist('Triangular', 'a', -0.5 * pitch_of(mt_x), 'b', 0, 'c', 0.5 * pitch_of(mt_x));
    ypd = makedist('Triangular', 'a', -0.5 * pitch_of(mt_y), 'b', 0, 'c', 0.5 * pitch_of(mt_y));
    zpd = makedist('Triangular', 'a', -0.5 * pitch_of(mt_s), 'b', 0, 'c', 0.5 * pitch_of(mt_s));
    % 2. generate samples based on low pass filters
    [x, y, z] = meshgrid(mt_x, mt_y, mt_s);
    x = repmat(reshape(x, 1, []), nsamples, 1);
    y = repmat(reshape(y, 1, []), nsamples, 1);
    z = repmat(reshape(z, 1, []), nsamples, 1);
    x = x + random(xpd, size(x));
    y = y + random(ypd, size(y));
    z = z + random(zpd, size(z));
    % 3. sample 
    vals = z.^(-2) .* interp2(mask_xx, mask_yy, mcode, x ./ z, y ./ z, 'nearest', 0);
    % 3. average samples
    vals = reshape(mean(vals, 1), [length(mt_y) length(mt_x) length(mt_s)]);
    kernel = permute(vals, [3 1 2]); %(z, y, x)        

    kernel = kernel * xx_pitch * yy_pitch * ss_pitch;
    
    kernel_range = [min(mt_x) max(mt_x) min(mt_y) max(mt_y) min(mt_s) max(mt_s)];
    
    if ~isempty(matname)
        save(matname, 'kernel', 'kernel_range', 'tt_y', 'tt_x', 'tt_s')
    end
end