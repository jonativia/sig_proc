function [locations, amplitudes, histogram, t_c] = get_locations_from_histogram(t_k, a_k, t0, tF, N, threshold, overlap)

if (nargin < 6 || nargin > 7)
    error('get_locations_from_histogram:err_arg', 'The number of input arguments is incorrect.')
elseif (nargin == 7 && (overlap < 0 || overlap > 1))
    error('get_locations_from_histogram:err_param', 'The overlap value must be a float in [0,1].')
elseif nargin == 6
    overlap = 0;
end
    
locations  = [];
amplitudes = [];

% Compute bin centers
t_c = linspace(t0, tF, N);

% Compute bin width
delta = (tF - t0) / (N-1);
delta = (1+overlap) * delta;

% Compute the histogram
histogram  = zeros(1, length(t_c));
for ith = 1 : N
    histogram(ith) = sum( (t_k > t_c(ith)-delta/2) & (t_k <= t_c(ith)+delta/2));
end

% Detect locations from peaks of the histogram above the threshold
for ith = 1 : N
    if histogram(ith) > threshold
        if ( ith < N && (histogram(ith) >= histogram(ith+1)) ) ...
        && ( ith > 1 && (histogram(ith) >  histogram(ith-1)) )
            t_i = t_k( (t_k > t_c(ith)-delta/2) & (t_k <= t_c(ith)+delta/2));
            a_i = a_k( (t_k > t_c(ith)-delta/2) & (t_k <= t_c(ith)+delta/2));
            locations  = [locations; mean(t_i)];
            amplitudes = [amplitudes; mean(a_i)];
        end
    end
end
