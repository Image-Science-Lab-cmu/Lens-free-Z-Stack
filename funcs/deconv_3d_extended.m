function [tt_hat, recons_time] = deconv_3d_extended(i_t, kernel, params)
recons_type = params.recons_type;
% i_t (1/z, y/z, x/z)
% kernel (1/z, y/z, x/z)

% normalize kernel first
k_norm = sqrt(sum(kernel.^2, 'all'));
k = kernel / k_norm;

if isfield(params, 'use_gpu'), use_gpu = params.use_gpu; else, use_gpu = 0; end

switch recons_type        
    case {'cgs-fft', 'fista'}
        % setup operators
        siz_i = size(i_t);
        siz_t = size(kernel) - siz_i + 1;
        if use_gpu,
            gpu = gpuDevice();
            reset(gpu);
            tmp = gpuArray(zeros(siz_t));
            k = gpuArray(k);
            [y0, Fk0, siz] = cgs_operator_fftconvn(tmp, 'A', [], k, []);
        else
            [y0, Fk0, siz] = cgs_operator_fftconvn(zeros(siz_t), 'A', [], k, []);
        end
        fA = @(tt) cgs_operator_fftconvn(reshape(tt, siz_t), ...
                        'A', Fk0, [], siz);
        fAadj = @(i_t) cgs_operator_fftconvn(reshape(i_t, siz_i), ...
                        'Adj', Fk0, [], siz);
        
        lambda = params.lambda;
        maxit = params.maxit; 
        if isfield(params, 't_init'), 
            t_init = params.t_init; % start from this solution;
        else
%             t_init = zeros(siz_t);
            t_init = fAadj(i_t);
        end

        tstart = tic;
        switch recons_type
            case 'cgs-fft'
                tol = params.tol;
                fCGS = @(x) vec(fAadj(fA(x))) + lambda .* vec(x);
                if use_gpu,
                    t_init = gpuArray(vec(t_init));
                    rhs = vec(fAadj(gpuArray(i_t)));
                else
                    t_init = vec(t_init);
                    rhs = vec(fAadj(i_t));
                end
                [tt_hat, flag, relres, iter] = cgs(...
                    fCGS, rhs, tol, maxit, [], [], t_init);
                fprintf('cgs exited with flag %d relres %.2e after %d iter\n', ...
                        flag, relres, iter);
            case 'fista'
                fA = @(tt) vec(cgs_operator_fftconvn(reshape(tt, siz_t), ...
                        'A', Fk0, [], siz));
                fAadj = @(i_t) vec(cgs_operator_fftconvn(reshape(i_t, siz_i), ...
                        'Adj', Fk0, [], siz)); 
                psi = @(x) vec(x);
                psi_adj = @(x) vec(x);
                i_t = vec(i_t);
                t_init = vec(t_init);
                if use_gpu,
                    i_t = gpuArray(i_t);
                    t_init = gpuArray(t_init);
                end
                tt_hat = fista_backtracking(i_t, lambda, psi, psi_adj, ...
                    fA, fAadj, maxit, t_init);

        end
        tt_hat = reshape(tt_hat, siz_t) / k_norm;
        if use_gpu
            tt_hat = gather(tt_hat);
        end

        recons_time = toc(tstart);
        total_time = seconds(recons_time);
        total_time.Format = 'hh:mm:ss';
        fprintf('%s took %s\n', ...
            recons_type, total_time);
        
    otherwise
        fprintf('Invalid deconv method %s\n', recons_type);
        recons = [];
        recons_time = 0;
        return;
end
end