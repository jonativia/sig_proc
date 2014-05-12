
clc, clear all, close all
format longg

% Signal's characteristics
N_tot = 1000;
K_tot = 100;
K     = N_tot / (2*K_tot);
N     = 2 * K^2;
T     = 1 / 16;
tau   = N * T;
P     = 2*K-1;
L     = P + 1;

% Temporal resolution
TTs   = 64;
T_s   = T / TTs;
len   = N_tot * T;
t_vec = (0 : T_s : len)';
L_t   = length(t_vec);

err_threshold = 10^-6;

% % Generate stream of Diracs
% locs = generate_diracs_locations(K, tau, 16*T_s, len, 1);
% 
% % Remove Diracs located near the borders
% locs(locs<tau)            = [];
% locs(locs>t_vec(end)-tau) = [];
% 
% % Generate amplitude of the Diracs
% amps = .6 + .6 * rand(size(locs));
% 
% % Bind locations to temporal grid
% i_locs    = round(locs/T_s) + 1;
% locs      = t_vec(i_locs);

load('thousand_diracs.mat', 'locs', 'amps')

n0 = floor((locs(1)-tau)/T)-1;
nF = ceil((locs(end)+tau)/T+L);

t_vec = ( n0*T : T_s : nF*T )';
n_vec = ( n0 : nF )';
t_n   = n_vec*T;

%% Count the rate of Diracs per tau interval
K_tau = zeros(length(t_n(1:end-N)), 1);
it = 1;
for t0 = t_n(1:end-N)'
    tF = t0 + tau;
    
    idx = find(locs>t0 & locs<=tF);
    while length(idx) > K
        locs(idx(end)) = [];
        amps(idx(end)) = [];
        idx = find(locs>t0 & locs<=tF);
    end
    
    K_tau(it) = sum(locs>t0 & locs<=tF);
    it        = it + 1;
end

plot_figs = true;
if plot_figs
    font_size = 16;
    
    i0 = 1;
    iF = 50000;
    % Plot Diracs and number of Diracs within 'tau' intervals
    figure
    set(gcf, 'Position', [50 50 560 420])
    subplot(211)
    stem(locs, amps, '^k', 'fill')
    hold on
    plot([4 4+tau], [1.4 1.4], '+-k', 'LineWidth', 2)
    text(5, 1.6, '\tau', 'FontSize', 20);
    hdl = xlabel('time (s)');
    set(hdl, 'FontSize', font_size)
%     axis([t_vec(i0) t_vec(iF) -.2 1.8])
    hdl = legend('Diracs', 'Location', 'NorthEast');
    set(hdl, 'FontSize', font_size)
    subplot(212)
    plot(t_n(1:end-N)', K_tau, 'k')
    hdl = xlabel('t_i');
    set(hdl, 'FontSize', font_size)
    hdl = legend('# of Diracs in [t_i, t_i + \tau]', 'Location', 'SouthEast');
    set(hdl, 'FontSize', font_size)
%     axis([t_vec(i0) t_vec(iF) -1 1.2*K])
end

%% Generate temporal signal
i_locs    = round((locs-t_vec(1))/T_s) + 1;
x         = zeros(size(t_vec));
x(i_locs) = amps;

% Construct the sampling kernel
m            = 0:P;
gamma        = (2*pi) / (N-P);
alpha_0      = -1j * gamma * P / 2;
alpha_m      = alpha_0 + 1j * m * gamma;
[phi, t_phi] = generate_e_spline(alpha_m, T_s, T, 'anticausal');
L_phi        = length(phi);
phi          = real(phi);

% Exponential reproducing coefficients
c_m_n = get_c_m_n_exp(alpha_m, 1:N, phi, t_phi, T);

% Sampling kernel corresponds to the time reversed version of phi
h   = phi(end:-1:1);
t_h = -t_phi(end:-1:1);

% Compute y_n convolving x(t) and phi(-t/T) and sampling at t=nT
y_t   = conv(x, h);
y_t   = y_t(1:end-L_phi+1);
y_n   = y_t(1:T/T_s:end);

% Sequential processing
tk_vec     = [];
ak_vec     = [];
win_idx    = [];
locations  = [];
histogram  = [];
t_c        = [];
for ith = 1:length(t_n)-N
    
    t0 = t_n(ith);
    tF = t_n(ith) + tau;
    tT = tF - P*T; % => perfect recovery interval goes from t0 to tT
    
    true_locs = locs(locs>t0 & locs<tF);
    true_amps = amps(locs>t0 & locs<tF);
    sub_amps  = amps(locs>t0 & locs<=tT);
    
%     disp(['perf reconstr. (' num2str(t0) ', ' num2str(tT) '], tau = (' num2str(t0) ', ' num2str(tF) ']'])
%     [true_locs true_amps]
    
    % Obtain samples that correspond to this time interval
    yi_n = y_n(ith+1:ith+N);
    
    % Remove previously detected Diracs' contribution
    nA   = n_vec(ith+1);
    nB   = n_vec(ith) + N;
    tA   = nA*T + t_phi(1);
    tt   = round(tk_vec(tk_vec>tA & tk_vec<tF)/T_s)*T_s;
    aa   = ak_vec(tk_vec>tA & tk_vec<tF);
    yc_n = yi_n;
    if ~isempty(tt)
        yy              = get_phi_tk_n_mat(phi,t_phi,tt,nA:nB,T,T_s) * aa;
        yc_n(1:nB-nA+1) = yi_n(1:nB-nA+1) - yy;
    end
    Pyi = yc_n' * yc_n / length(yc_n);
    
    % Compute the moments and construct matrix
    sc_m  = c_m_n * yc_n;
%     S     = toeplitz(sc_m(round(P/2)+1:end), sc_m(round(P/2)+1:-1:1));
%     k_est = rank(S);
%     si_m  = c_m_n * yi_n;
%     S2    = toeplitz(si_m(round(P/2)+1:end), si_m(round(P/2)+1:-1:1));
%     k2    = rank(S2);
    
    % Oracle to get number of Diracs within the time interval
    k_tot = sum(locs>t0 & locs<tF);
    k_int = sum(locs>t0 & locs<=tT);
    k_det = sum(tk_vec>t0 & tk_vec<=tT);
    k     = k_int - k_det;

    if Pyi > err_threshold
        for k = 1 : K
            % Retrieve locations with annihilating filter
            [u_k b_k]      = locate_diracs(sc_m, k);
            p_u_k          = angle(u_k);
            p_u_k(p_u_k<0) = p_u_k(p_u_k<0) + 2*pi;
            t_k            = t_n(ith) + T * p_u_k(:) / gamma;

            % Retrieve amplitudes with retrieved locations and samples y_n
            t_k     = round(t_k/T_s) * T_s;
            phi_mat = get_phi_tk_n_mat(phi, t_phi, t_k, ...
                                       nA:nB, T, T_s);
            a_k     = phi_mat \ yc_n;
            
            % Remove Diracs with close to zero amplitude
            t_k(a_k<0.1) = [];
            a_k(a_k<0.1) = [];
            if isempty(t_k)
                continue;
            end
            phi_mat = get_phi_tk_n_mat(phi, t_phi, t_k, ...
                                       nA:nB, T, T_s);
            
            % Reconstruct samples from retrieved locations to check if it is
            % correct
            yyii_n = phi_mat * a_k;
            err    = sqrt((yc_n-yyii_n)' * (yc_n-yyii_n));
            
            % Sort locations
            [~, idx] = sort(t_k);
            t_k      = t_k(idx);
            a_k      = a_k(idx);
        
%            if tF > 243.7
%                disp(['perf reconstr. (' num2str(t0) ', ' num2str(tT) '], tau = (' num2str(t0) ', ' num2str(tF) ']'])
%                [true_locs [zeros(k_det,1);t_k;zeros(k_tot-k_int,1)] true_amps [zeros(k_det,1);a_k;zeros(k_tot-k_int,1)]]
%             end

            if err < err_threshold
                tk_vec = [tk_vec; t_k];
                ak_vec = [ak_vec; a_k];
            
%                if tF > 243.7
%                     disp(['Dirac recovered. t_k = [' num2str(t_k') '], a_k = [' num2str(true_amps(k_det+1:end)') ']'])
%                 end
                break;
            end
        end

    end
    
end

%% Compare the detected spikes with the real spikes
num_sp   = length(locs);
hit_sp   = false(num_sp, 1);
sspp_ids = [];
sspp_cpy = tk_vec;
delta_t  = 50*T_s;
for ith_sp = 1 : num_sp
    t_i  = locs(ith_sp);
    inds = find(sspp_cpy > (t_i - delta_t/2) & sspp_cpy < (t_i + delta_t/2));
    
    if ~isempty(inds)
        hit_sp(ith_sp) = true;
        sspp_ids       = [sspp_ids; find(tk_vec == sspp_cpy(inds(1)))];
        
        % Remove this spike from detected spikes
        sspp_cpy(inds(1)) = [];
        
        if length(inds) > 1
            warning(['More than one spike detected in the neighbourhood of Dirac at ' num2str(t_i)]);
        end
    end
end

% Accuracy of detected spikes
hit_rate  = sum(hit_sp) / length(hit_sp) * 100;
false_pos = length(sspp_cpy);
if ~isempty(tk_vec)
    mse       = mean( (locs(hit_sp) - tk_vec(sspp_ids)).^2 );
    std_dev   = sqrt(mse);
end
disp('')
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp(['Total number of real Diracs     : ' num2str(num_sp)])
disp(['Total number of detected Diracs : ' num2str(length(tk_vec))])
disp(['Real Diracs detected            : ' num2str(sum(hit_sp))])
if ~isempty(tk_vec)
    disp(['MSE of locations                : ' num2str(mse)])
    disp(['Std deviation of locations      : ' num2str(std_dev)])
end
disp(['Dirac detection rate            : ' num2str(hit_rate) '%'])
disp(['False positives                 : ' num2str(false_pos)])
disp(' ')


%% Plots
if plot_figs
    font_size = 16;
    marker_size = 10;
    
    t1 = tau;
    t2 = 12*tau;
    dt = t_vec(end) - t_vec(1);
    
    % Plot signal and samples
    figure
    set(gcf, 'Position', [100 100 560*4/5 420*4/5])
    stem(locs, amps, '^r', 'fill', 'MarkerSize', marker_size)
    hold on
    stem(tk_vec, ak_vec, 'pb', 'fill', 'MarkerSize', marker_size)
    plot([6 6+tau], [1.4 1.4], '+-k', 'LineWidth', 2)
    text(7.2, 1.5, '\tau', 'FontSize', 20);
    hdl = legend('Original Diracs', 'Reconstructed Diracs');
    set(hdl, 'FontSize', font_size)
    hdl = xlabel('Time (s)');
    set(hdl, 'FontSize', 16)
    axis([t1 t2 -.2 1.6])    
end
