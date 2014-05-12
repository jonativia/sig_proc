function c_m_n = get_c_m_n_exp(alpha_m, n, phi, t_phi, T, t_0)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2011
%
% Date        : 21/11/2011
% Supervisor  : Dr Pier Luigi Dragotti
% Author      : Jon Onativia
%
% File        : get_c_m_n_exp.m
% -------------------------------------------------------------------------
% Compute the c_m_n coefficients to reproduce exponentials with parameters
% given by the alpha_m vector using the exponential reproducing kernel phi:
%   exp(alpha_m*t) = sum_n ( c_m_n * phi(t-n) )
%
% USAGE:
%  c_m_n = get_c_m_n_exp(alpha_m, n, phi, t[, T, t_0])
%
% INPUT:
%  - alpha_m : Vector of size M with the parameters of the exponentials to 
%              be reproduced.
%  - n       : Vector of size N with the values where the summation will be
%              evaluated.
%  - phi     : Exponential reproducing kernel.
%  - t_phi   : Time stamps of the kernel.
%  - T       : Optional argument. Scale factor. Default T = 1.
%  - t_0     : Optional argument. t value where c_m_0 will be evaluated. Default t_0 = 0.
%
% OUTPUT:
%  - c_m_n   : Coefficients to reproduce the exponentials.
%

if nargin < 4 || nargin > 6
    error('get_c_m_n_exp:err_arg', 'The number of input arguments is incorrect.')
elseif nargin < 6
    t_0 = 0;
    if nargin < 5
        T = 1;
    end
end

% Rearrange the arguments (n row vector, alpha_m column vector)
n       = n(:).';
n_len   = length(n);
alpha_m = alpha_m(:);
T_s     = t_phi(2) - t_phi(1);

% Kernel's boundaries
t_1 = t_phi(1) / T;
t_2 = t_phi(end) / T;

% Compute c_m_0 vector
l     = ceil(t_0/T - t_2) : floor(t_0/T - t_1);
idx   = round( (t_0 - T * (t_1 + l)) / T_s) + 1;
phi_l = phi(idx);
num   = exp(alpha_m * t_0 / T);
den   = exp(alpha_m * l) * phi_l;
c_m_0 = num ./ den;

% Compute the remaining c_m_n from c_m_0
exp_mat = exp(alpha_m * n);
c_m_n   = exp_mat .* repmat(c_m_0, 1, n_len);
