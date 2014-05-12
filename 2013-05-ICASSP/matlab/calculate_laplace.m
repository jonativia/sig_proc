function [beta_alpha_vec] = calculate_laplace(alpha_vec, gamma_vec)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2011
%
% Date        : 19/09/2012
% Supervisor  : Dr Pier Luigi Dragotti
% Author      : Jose Antonio Uriguen
%
% File        : calculate_laplace.m
% -------------------------------------------------------------------------
% Calculate the Laplace transform of the default ESpline at values s =
% gamma_vec:
%   prod((1-exp(alpha_vec-gamma_vec))./(gamma_vec-alpha_vec))
%
% USAGE:
%  [beta_alpha_vec] = calculate_laplace(alpha_vec, gamma_vec)
%
% INPUT:
%  - alpha_vec      : Vector of P+1 alpha values of the E-spline.
%  - gamma_vec      : Vector of s values of the E-spline.
%
% OUTPUT:
%  - beta_alpha_vec : Laplace transform of the ESpline at s = gamma_vec.
%

[X, Y] = meshgrid(alpha_vec, gamma_vec);
temp = (1-exp(X-Y))./(Y-X);
temp(isnan(temp)) = 1;
if numel(alpha_vec)>1
    beta_alpha_vec = prod(temp.');
else
    beta_alpha_vec = temp.';
end


end