function h = annihilating_filter(tau, K)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2011
%
% Date        : 21/11/2011
% Supervisor  : Dr Pier Luigi Dragotti
% Authors     : Jon Onativia
%
% File        : annihilating_filter.m
% -------------------------------------------------------------------------
%
% h = annihilating_filter(tau, K)
%
% INPUT:
%  - tau : vector of length N + 1 with the first moments of signal x(t)
%  - K   : number of Diracs
%
% OUTPUT:
%  - h   : coefficients of the annihilating filter
%
% It is assumed that the moments come from a signal x(t) of the form:
%         K-1
%  x(t) = sum ( a_k * delta(t - t_k) )
%         k=0
%
% The length of the input tau has to satisfy N+1 >= 2K.
%

% Check input consistency
N = length(tau) - 1;
if N < 2*K-1
    error('Invalid length moments vector, at least 2K moments are needed.')
end

A = toeplitz(tau(K:N), tau(K:-1:1));
B = -tau(K+1:N+1);

h = linsolve(A, B);

h = [1; h];

