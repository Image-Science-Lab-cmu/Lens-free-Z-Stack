function [x_hat, obj] = fista_backtracking(y, beta, psi, psi_adj, A, A_adj, max_iter, x_hat0)
%Solves \min \| y - A(x) \|^2 +beta \| Psi(x) \|_1 
%
% Input
%   y -- measurements
%   A -- function handle of A
%   A_adj -- function handle of A adjoint
%   psi -- function handle of psi
%   psi_adj -- function handle op psi adjoint
%   beta --- scaling parameter
%   max_itere --- maximum iterations
%
% Output
%   x -- solution to optimization problem


if nargin < 8,
%     x_hat0 = zeros(size(A_adj(y)));
    x_hat0 = A_adj(y);
end
y_hat = x_hat0;

gradient = @(x) -A_adj(y-A(x));
least_square = @(x) 0.5*sum(vec(y-A(x)).^2);

st0 = 1;

t0 = 2;
t_beta = 0.5;
obj_old = 1e100;

for iter = 1:max_iter
    t = t0;
    g = gradient(y_hat);
    while true
        Gt = (y_hat - psi_adj(soft_thresholding(psi(y_hat-t*g), t*beta)))/t;
        test_lhs = least_square(y_hat - t*Gt);
        test_rhs = least_square(y_hat) - t*vec(g)'*vec(Gt) + 0.5*t*sum(vec(Gt).^2);
        if test_lhs > test_rhs
            t = t_beta * t;
        else
            break;
        end
    end
    % update x_hat 
    x_hat = psi_adj(soft_thresholding(psi(y_hat-t*g), t*beta));

    st1 = (1+sqrt(1+4*st0^2))/2;
    
    y_hat = x_hat + (x_hat - x_hat0)*(st0-1)/st1;
    
    st0 = st1;
    x_hat0 = x_hat;
    
%    imshow(x_hat); drawnow
    
    % print status
%     obj = 0.5 * norm(vec(y-A(x_hat)))^2 + beta * norm(vec(psi(x_hat)),1);
    data_term = 0.5 * norm(vec(y-A(x_hat)))^2;
    prior_term = beta * norm(vec(psi(x_hat)),1);
    obj = data_term + prior_term;
    fprintf('iter = %d data = %f + prior = %f = obj = %f\n', iter, data_term, prior_term, obj);
    if (2*abs(obj_old - obj)/(1e-8+obj+obj_old)) < 1e-6
        break;
    end
    obj_old = obj;
end
end

function y = vec(x)
y = x(:);
end

function y = soft_thresholding(x, beta)
y = max(0, x - beta) - max(0, -x-beta);
end
