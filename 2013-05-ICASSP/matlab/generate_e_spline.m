function [phi, t] = generate_e_spline(alpha_vec, T_s, T, mode)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2011
%
% Date        : 21/11/2011
% Supervisor  : Dr Pier Luigi Dragotti
% Author      : Jon Onativia
%
% File        : generate_e_spline.m
% -------------------------------------------------------------------------
% Generate the exponential spline of order P+1 corresponding to a vector of
% alpha values and with a given temporal resolution. The resulting spline 
% is obtained in time domain computing the P convolutions of the P+1 zero 
% order E-splines:
%   phi_a_vec(t) = phi_a_0(t) * phi_a_1(t) * ... * phi_a_N(t)
%
% USAGE:
%  [phi, t] = generate_e_spline(alpha_vec, T_s[, T, mode])
%
% INPUT:
%  - alpha_vec : Vector of P+1 alpha values of the E=spline.
%  - T_s       : Time resolution of the spline.
%  - T         : Optional argument. Scale factor. Default T = 1.
%  - mode      : Optional argument. 'causal', 'symmetric' or 'anticausal'. 
%                Default mode = 'causal'.
%
% OUTPUT:
%  - phi       : Vector of size (P+1)/T + 1 with the values of the
%                E-spline.
%  - t         : Time stamps of the corresponding values of the phi vector.
%

if nargin < 2 || nargin > 4
    error('generate_e_spline:err_arg', 'The number of input arguments is incorrect.')
elseif nargin < 4
    mode = 'causal';
    if nargin < 3
        T = 1;
    end
end

N = length(alpha_vec) - 1;

% Convert alpha_vec into a row vector
alpha_vec = alpha_vec(:).';

% Apply scaling factor
T_s = T_s / T;

t_phi_1        = 0;
t_phi_2        = 1;
t_phi          = (t_phi_1:T_s:t_phi_2)';
sub_phi        = exp(t_phi * alpha_vec);
sub_phi(end,:) = 0;

phi   = [0; sub_phi(1:end-1,1)];
t_0   = t_phi(1);
t_end = t_phi(end);
for i = 1:N
    t_0   = t_0 + t_phi(1);
    t_end = t_end + t_phi(end);
    phi   = T_s * conv(phi, sub_phi(:,i+1));
end

t = (t_0:T_s:t_end)';
t = t * T;

if strcmp(mode, 'symmetric')
    t_mid      = (t(end) - t(1)) / 2;
    t          = t - t_mid;
    [~, i_max] = max(phi);
    if phi(i_max) ~= phi(t == 0)
        t = t - t(i_max);
    end
elseif strcmp(mode, 'anticausal')
    phi = phi(end:-1:1);
    t   = -t(end:-1:1);
end
