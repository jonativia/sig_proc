function [t_k a_k] = locate_diracs(tau, K)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2011
%
% Date        : 21/11/2011
% Supervisor  : Dr Pier Luigi Dragotti
% Authors     : Jon Onativia
%
% File        : locate_diracs.m
% -------------------------------------------------------------------------
%
% [t_k a_k] = locate_diracs(tau, K)
%
% INPUT:
%  - tau : vector of length N + 1 with the first moments of signal x(t)
%  - K   : number of Diracs
%
% OUTPUT:
%  - t_k : location of the diracs
%  - a_k : amplitude of the diracs
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

% Obtain the annihilating filer, zeros of the filter correspond to the t_k
h = annihilating_filter(tau, K);

% Find the locations t_k
t_k = roots(h);

% Find the amplitudes a_k
A = zeros(K, K);
for i = 0:K-1
    A(i+1,:) = t_k(1:K).^i;
end
B = tau(1:K);
B = B(:);
a_k = linsolve(A, B);

% Return results as row vectors
t_k = t_k.';
a_k = a_k.';
