function [phi, t, d_l] = generate_e_spline_freq(alpha_vec, T_s, T, moms, c0, mode)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2011
%
% Date        : 19/09/2012
% Supervisor  : Dr Pier Luigi Dragotti
% Author      : Jose Antonio Uriguen
%
% File        : generate_e_spline.m
% -------------------------------------------------------------------------
% Generate the exponential spline of order P+1 corresponding to a vector of
% alpha values and with a given temporal resolution. The resulting spline 
% is obtained in frequency domain by implementing the Fourier transform  
% and using the IFFT:
%   phi_a_vec(t) = ifft((1-exp(alpha_vec-1j*w))./(1j*w-alpha_vec))
%
% USAGE:
%  [phi, t] = generate_e_spline_freq(alpha_vec, T_s[, T, mode])
%
% INPUT:
%  - alpha_vec : Vector of P+1 alpha values of the E-spline.
%  - T_s       : Time resolution of the spline.
%  - T         : Optional argument. Scale factor. Default T = 1.
%  - moms      : Optional argument. Default ESpline (0) or exponential MOMS (1).
%  - c0        : Optional argument. Value for the eMOMS varphi(wm) = c0.
%  - mode      : Optional argument. 'causal', 'symmetric' or 'anticausal'. 
%                Default mode = 'causal'.
%
% OUTPUT:
%  - phi       : Vector of size (P+1)/T + 1 with the values of the
%                E-spline.
%  - t         : Time stamps of the corresponding values of the phi vector.
%  - d_l       : Coefficients d_l of the eMOMS phi(w) = beta(w) gamma(w), 
%                for gamma(w)=sum(d_l w^l) to interpolate the set of points
%                (wm,beta^-1(wm)). Here the parameters satisfy alpha_vec = jwm.
%

if nargin < 2 || nargin > 6
    error('generate_e_spline:err_arg', 'The number of input arguments is incorrect.')
elseif nargin < 6
    mode = 'causal';
    if nargin < 5
        c0 = exp(alpha_vec * floor(length(alpha_vec)/2));
        if nargin < 4
            moms = 0;
            if nargin< 3
                T = 1;
            end
        end
    end
end

P = length(alpha_vec) - 1;

% Convert alpha_vec into a row vector
alpha_vec = alpha_vec(:).';

% Apply scaling factor
T_s = T_s / T;

t_0   = 0;
t_end = P+1;
t     = (t_0:T_s:t_end)' * T;


% Default coefficients
d_l = 1;

length_Phi = (P+1)/T_s+1;
M_fft      = 2^nextpow2(length_Phi);
w          = 2*pi*[0:M_fft/2-1, -M_fft/2:1:-1]/(M_fft*T_s);

[X, Y]            = meshgrid(alpha_vec, 1j*w);
temp              = (1-exp(X-Y))./(Y-X);
temp(isnan(temp)) = 1;

if numel(alpha_vec) > 1
    beta_w = prod(temp.');
else
    beta_w = temp.';
end


if moms
    d_l = zeros(1,P+1);
    for m = 1:P+1
        m0     = setdiff(1:P+1, m);
        lambda = prod(1-exp(-alpha_vec(m)+alpha_vec(m0)));
        d_l    = d_l + 1/c0(m)/lambda * poly(alpha_vec(m0));
    end
    d_l = fliplr(d_l);
    
    gamma_w = 0;
    for i = 0:P
        gamma_w = gamma_w + d_l(i+1)*(1j*w).^i;
    end
    beta_w  = beta_w .* gamma_w;
end

phi = M_fft*ifft(beta_w,M_fft);
phi = phi(1:length_Phi).';


area   = sum(real(phi)*T_s);
beta_0 = d_l(1)*calculate_laplace(alpha_vec, 0);
phi    = phi/area*beta_0;


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
