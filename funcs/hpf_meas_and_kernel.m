function [i_tf, kernel_f] = hpf_meas_and_kernel(i_t, kernel, use_filter, hpf_params)
    if use_filter == 1
        switch hpf_params.type
            case '2d'
                fprintf('2d filter meas & kernel with %s\n', hpf_params.name);
                hpf = hpf_params.hpf;
                % filter i_t (z, y, x)
                i_tf = zeros(size(i_t));
                for i = 1:size(i_t, 1)
                    i_tf(i, :, :, :) = imfilter(squeeze(i_t(i,:,:,:)), hpf, 'symmetric', 'same');
                    % zero padding creates high values on border/corners
                end
                % filter kernel (z, y, x)
                kernel_f = zeros(size(kernel));
                for i = 1:size(kernel, 1)
                    kernel_f(i, :, :) = imfilter(squeeze(kernel(i,:,:)), hpf, 0, 'same');
                end
                return
            case 'mean_removal'
                fprintf('remove col & row mean\n');
                % remove col & row mean in i_t
                i_tf = i_t - mean(i_t, 2) - mean(i_t, 3);
                kernel_f = kernel - mean(kernel, 2) - mean(kernel, 3);
                return
            case 'no_filter'
                % valid type, do nothing
            otherwise
                fprintf('Invalid hpf type %s\n', hpf_params.type);
        end
    end
    i_tf = i_t;
    kernel_f = kernel;
end