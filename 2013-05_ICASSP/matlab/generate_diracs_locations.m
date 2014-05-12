function [locs] = generate_diracs_locations(K, tau, T_s, duration, min_distance)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2012
%
% Date        : 26/10/2012
% Supervisor  : Dr Pier Luigi Dragotti
% Author      : Jon Onativia
%
% File        : generate_diracs_locations.m
% -------------------------------------------------------------------------
% 
%
% USAGE:
%  [locs] = generate_diracs_locations(K, tau, T, length)
%
% INPUT:
%  - K        : Dirac rate per tau time interval.
%  - tau      : Time interval with K Diracs.
%  - T_s      : Temporal resolution.
%  - duration : Total length of the stream.
%
% OUTPUT:
%  - locs     : Locations of the Diracs.
%

if nargin == 4
    min_distance = 0;
end

% Locations of first time interval
locs = sort(rand(K, 1) * tau);
while min(abs(diff(locs))) < min_distance*T_s
    locs = sort(rand(K, 1) * tau);
end

% Generate Diracs for the rest of the interval
for t_1 = tau : T_s : duration
    t_0 = t_1 + T_s - tau;
    
    % Count Diracs in previous tau-T_s seconds and randomly generate new
    % Diracs
    k        = round(rand(1) * (K - sum(locs > t_0)));
    new_locs = sort(t_1 + rand(k, 1)*T_s);
    while min(abs(diff([locs(end); new_locs]))) < min_distance*T_s
        k        = round(rand(1) * (K - sum(locs > t_0)));
        new_locs = sort(t_1 + rand(k, 1)*T_s);
    end
    locs     = [locs; new_locs];    
end

locs(locs>duration) = duration;
