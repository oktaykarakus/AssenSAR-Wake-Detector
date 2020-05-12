function [x, v] = GMC_regularisation(y, A, AH, rho, lam, gamma)

% [x, v] = GMC_regularisation(y, A, AH, rho, lam, gamma)
%
% Original name:
% srls_GMC: Sparse-Regularized Least Squares with generalized MC (GMC) penalty
%
% Saddle point problem:
%
% argmin_x  argmax_v { F(x,v) =
%  1/2 ||y - A x||^2 + lam ||x||_1 - gamma/2 ||A(x-v)||_2^2 - lam ||v||_1 }
%
% INPUT
%   y 	    data
%   A, AH   operators for A and A^H
%   rho     rho >= maximum eigenvalue of A^H A
%   lam     regularization parameter, lam > 0
%   gamma   0 <= gamma < 1
%
% OUTPUT
%   x, v
% Ivan Selesnick
% May 2016
% Revised: July 2016
% Reference:
% I. Selesnick, 'Sparse Regularization via Convex Analysis'
% IEEE Transactions on Signal Processing, 2017.
%
% This version is slightly modified by O.Karakus in October 2018. For
% original function please refer to Ivan Selesnick's web site.
%
% Algorithm: Forward-backward, Theorem 25.8 in Bauschke and Combettes (2004)

MAX_ITER = 500;
TOL_STOP = 5e-3;

% soft thresholding for complex data
% soft = @(x, T) max(1 - T./abs(x), 0) .* x;

% soft thresholding for real data
soft = @(t, T) max(t - T, 0) + min(t + T, 0);

% rho = max(eig(A'*A));
mu = 1.9 / ( rho * max( 1,  gamma / (1-gamma) ) );

AHy = AH(y);

% initialization
x = zeros(size(AHy));
v = zeros(size(AHy));

iter = 1;
old_x = x;
delta_x = inf;
while (delta_x(iter) > TOL_STOP) && (iter < MAX_ITER)    
    iter = iter + 1;
    
    % update x
    zx = x - mu * ( AH(A(x + gamma*(v-x))) - AHy );
    zv = v - mu * ( gamma * AH(A(v-x)) );

    % update v
    x = soft(zx, mu * lam);
    v = soft(zv, mu * lam);
    delta_x(iter) = max(abs( x(:) - old_x(:) )) / max(abs(old_x(:)));
    old_x = x;
    if or((iter > 49 && delta_x(iter) > 1e-1), (iter > 119 && delta_x(iter) > 1e-2))
        break
    end
end
