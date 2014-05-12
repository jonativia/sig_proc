
clc, clear all, close all

% Signal's characteristics
K   = 5;
tau = K^2 / 8;
N   = 2 * K^2;
T   = 1 / 16;
SNR = 10;

% Parameters
% P     = 2*K-1;
P = 22;
TTs   = 64;
T_s   = T / TTs;
len   = 640;
t_vec = (0 : T_s : len)';
L_t   = length(t_vec);
thres = 0.25;

% Generate stream of Diracs
locs = generate_diracs_locations(K, tau, 16*T_s, len, 6);
disp('Diracs locations generated')
locs(locs<tau)            = [];
locs(locs>t_vec(end)-tau) = [];
amps   = .6 + .6 * rand(size(locs));
i_locs = round(locs/T_s) + 1;
locs   = t_vec(i_locs);

% load('noisy_diracs.mat', 'locs', 'amps', 'y_n', 'e_n')

compute_post_histogram = false;

% Count the rate of Diracs per tau interval
K_tau = zeros((len-tau)/T, 1);
it = 1;
for t0 = 0:T:len-tau
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

% Generate the continuous-time signal
i_locs    = round(locs/T_s) + 1;
x         = zeros(size(t_vec));
% x(i_locs) = 1;
x(i_locs) = amps;


plot_figs = false;
if plot_figs
    font_size = 16;
    
    % Plot Diracs and number of Diracs within 'tau' intervals
    figure
    set(gcf, 'Position', [50 50 560 420])
    subplot(211)
    stem(locs, amps, '.')
    hold on
    plot([4 4+tau], [1.2 1.2], '+-k', 'LineWidth', 2)
    text(5, 1.35, '\tau', 'FontSize', 20);
    hdl = xlabel('time (s)');
    set(hdl, 'FontSize', font_size)
    axis([t_vec(1) t_vec(end) -.2 1.5])
    hdl = legend('Diracs', 'Location', 'NorthEast');
    set(hdl, 'FontSize', font_size)
    subplot(212)
    plot(0:T:len-tau, K_tau)
    hdl = xlabel('t_i (s)');
    set(hdl, 'FontSize', font_size)
    hdl = legend('# of Diracs in [t_i, t_i + \tau]', 'Location', 'SouthEast');
    set(hdl, 'FontSize', font_size)
    axis([t_vec(1) t_vec(end) -1 1.2*K])
end

% Construct the sampling kernel
m            = 0:P;
% gamma        = (2*pi) / (N-P);
% alpha_0      = -1j * gamma * P / 2;
% alpha_m      = alpha_0 + 1j * m * gamma;
% [phi, t_phi] = generate_e_spline(alpha_m, T_s, T, 'anticausal');
% phi          = real(phi);
gamma        = 2 * pi / (P+1);
alpha_0      = -1j * P * pi / (P + 1);
alpha_m      = alpha_0 + 1j * gamma * m;
c_0          = exp(alpha_m * floor((P+1)/2));
[phi, t_phi] = generate_e_spline_freq(alpha_m, T_s, T, 1, c_0, 'anticausal');
Lden = P+1;
mid  = floor(Lden/2);
t_diric = ((0:T_s/T:P+1)-mid)*2*pi/(P+1);
b   = diric(t_diric, (P+1));
phi = real(b.');

% Exponential reproducing coefficients
c_m_n = get_c_m_n_exp(alpha_m, 1:N, phi, t_phi, T);

% Sampling kernel corresponds to the time reversed version of phi
h     = phi(end:-1:1);
t_h   = -t_phi(end:-1:1);
L_phi = length(phi);

% figure
% set(gcf, 'Position', [100 100 280 210])
% plot(t_h, h, 'k')

% Compute the weigthing mask
x1   = [0; ones(size((T_s:T_s:(N-P)*T-T_s)'))];
y1   = conv(x1, h);
y1   = y1 / max(y1);
mask = y1(1+TTs:TTs:end);
half_mask        = mask;
L_2              = floor(length(half_mask) / 2);
half_mask(1:L_2) = 1;

% Compute y_n convolving x(t) and phi(-t/T) and sampling at t=nT
y_t   = conv(x, h);
y_t   = y_t(1:end-L_phi+1);
y_n   = y_t(1:T/T_s:end);
t_n   = (0 : T : len)'; % time stamps of samples y_n
n_vec = (0 : length(t_n))';

% Add noise
P_y   = y_n' * y_n / length(y_n);
e_n   = randn(size(y_n));
P_e   = e_n' * e_n / length(e_n);
A_e   = sqrt(10^(-SNR/10) * P_y / P_e);
e_n   = A_e * e_n;
sigma = sqrt((e_n' * e_n) / length(e_n));

tic
% Sequential processing
hist_res   = 16; % in terms of T_s
bins_per_T = TTs / hist_res;
M          = round((P+1)/2);
L          = P - M;
tk_vec     = [];
ak_vec     = [];
win_idx    = [];
seq_locs   = [];
seq_amps   = [];
seq_hist   = [];
seq_hist_t = [];
max_detect = N-P;
threshold  = thres * max_detect;
for ith = 1:length(t_n)-N
    
    t0 = t_n(ith);
    tF = t_n(ith) + tau;
    tT = tF - P*T; % => perfect recovery interval goes from t0 to tT
    
    true_locs = locs(locs>t0 & locs<tF);
    true_amps = amps(locs>t0 & locs<tF);
    
    % Obtain samples that correspond to this time interval
    yi_n = y_n(ith+1:ith+N) + e_n(ith+1:ith+N);
    
    % Remove previous Diracs' contribution
    nA = n_vec(ith+1);
    nB = n_vec(ith) + N;
%     tA = nA*T + t_phi(1);
%     tt = round(locs(locs>tA & locs<=t0)/T_s)*T_s;
%     yy = zeros(nB-nA+1, 1);
%     for n = nA : nB
%         [~,idx,~] = intersect(t_phi+n*T, tt);
%         if ~isempty(idx)
%             yy(n-nA+1) = sum(phi(idx));
%         end
%     end
%     yi_n(1:nB-nA+1) = yi_n(1:nB-nA+1) - yy;
    
    % Apply mask and compute moments
%     opt = 'no_mask';
    opt = 'mask';
%     opt = 'half_mask';
    switch opt
        case 'mask'
            yi_n = yi_n .* mask;
        case 'half_mask'
            yi_n = yi_n .* half_mask;
    end
    si_m = c_m_n * yi_n;
    
    % Oracle to get number of Diracs within the time interval
%     k = sum(locs>t0 & locs<=tT);
    k = K;

    % Retrieve locations with matrix pencil denoising
    u_k            = acmp_p(si_m, k, round(P/2), P, 1);
%     R   = calculate_r(M, P, c_m_n, sigma);
%     u_k = acmp_weightQ(si_m, R, K, L, M, P, 1);
    p_u_k          = angle(u_k);
    p_u_k(p_u_k<0) = p_u_k(p_u_k<0) + 2*pi;
    t_k            = t_n(ith) + T * p_u_k(:) / gamma;
    
    % Retrieve locations with Cadzow denoising
%     cad_it  = 10;
%     A_m     = toeplitz(si_m(floor(length(si_m)/2)+1:end), ...
%                        si_m(floor(length(si_m)/2)+1:-1:1));
%     A_m     = cadzow(A_m, K, cad_it);
%     sc_m    = [fliplr(A_m(1,:)) A_m(2:end,1).'];
%     sc_m    = sc_m.';
%     [u_k ~] = locate_diracs(sc_m, K);
%     p_u_k          = angle(u_k);
%     p_u_k(p_u_k<0) = p_u_k(p_u_k<0) + 2*pi;
%     t_k            = t_n(ith) + T * p_u_k(:) / gamma;

    % Retrieve amplitudes with retrieved locations and samples y_n
    t_k     = round(t_k/T_s) * T_s;
    phi_mat = get_phi_tk_n_mat(phi, t_phi, t_k, ...
                               nA:nB, T, T_s);
    a_k     = phi_mat \ yi_n;

    % Sort locations
    [~, idx] = sort(t_k);
    t_k      = t_k(idx);
    a_k      = a_k(idx);
    
    % Reject Diracs with wrong amplitudes
%     t_k(a_k<0.3 | a_k>1.7) = [];
%     a_k(a_k<0.3 | a_k>1.7) = [];            
    t_k(abs(a_k)<0.25) = [];
    a_k(abs(a_k)<0.25) = [];
    
%     % Remove retrieved locations in second half of the time interval (noisy
%     % locations)
%     t_k(t_k > (t0+tF)/2) = [];

    if ~isempty(t_k)
        tk_vec  = [tk_vec; t_k];
        ak_vec  = [ak_vec; a_k];
        win_idx = [win_idx; ones(length(t_k),1)*ith];
    end

    % Retrieve Dirac in ]t0,t0+T] from histogram
    [locats, amplis, histo, t_bins] = ...
        get_locations_from_histogram(tk_vec, ak_vec, ...
        t0-T/bins_per_T, t0+T, bins_per_T+2, threshold);
    seq_locs   = [seq_locs; locats];
    seq_amps   = [seq_amps; amplis];
    seq_hist   = [seq_hist histo(2:end-1)];
    seq_hist_t = [seq_hist_t t_bins(2:end-1)];
end

% Retrieve Diracs of last interval
t_last = t0;
for t0 = t0 : T : len-T
    [locats, amplis, histo, t_bins] = get_locations_from_histogram(tk_vec, ak_vec, ...
                                             t0-T/bins_per_T, t0+T, bins_per_T+2, threshold);
    seq_locs   = [seq_locs; locats];
    seq_amps   = [seq_amps; amplis];
    seq_hist   = [seq_hist histo(2:end-1)];
    seq_hist_t = [seq_hist_t t_bins(2:end-1)];
end
toc
%% Compare the detected spikes with the real spikes
num_sp   = length(locs);
hit_sp   = false(num_sp, 1);
sspp_ids = [];
sspp_cpy = seq_locs;
delta_t  = 1*T;
for ith_sp = 1 : num_sp
    t_i  = locs(ith_sp);
    inds = find(sspp_cpy > (t_i - delta_t/2) & sspp_cpy < (t_i + delta_t/2));
    
    if ~isempty(inds)
        hit_sp(ith_sp) = true;
        
        sspp_ids       = [sspp_ids; find(seq_locs == sspp_cpy(inds(1)))];
        
        % Remove this spike from detected spikes
        sspp_cpy(inds(1)) = [];
        
        if length(inds) > 1
            warning(['More than one spike detected in the neighbourhood of Dirac at ' num2str(t_i)]);
        end
    end
end
seq_hit_locs = seq_locs(sspp_ids);
seq_fp_locs  = sspp_cpy;

% Accuracy of detected spikes
hit_rate  = sum(hit_sp) / length(hit_sp) * 100;
false_pos = length(seq_fp_locs);
if ~isempty(seq_locs)
    mse       = mean( (locs(hit_sp) - seq_hit_locs).^2 );
    std_dev   = sqrt(mse);
end
disp('')
disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
disp(['Total number of real Diracs     : ' num2str(num_sp)])
disp(['Total number of detected Diracs : ' num2str(length(seq_locs))])
disp(['Real Diracs detected            : ' num2str(sum(hit_sp))])
if ~isempty(seq_locs)
    disp(['MSE of locations                : ' num2str(mse)])
    disp(['Std deviation of locations      : ' num2str(std_dev)])
end
disp(['Dirac detection rate            : ' num2str(hit_rate) '%'])
disp(['False positives                 : ' num2str(false_pos)])
disp(' ')

%% A posteriori histogram (should give same resutl as sequential one)
if compute_post_histogram
    % Retrieve locations from peaks of the locations histogram
    t0         = t_vec(1);
    tF         = t_vec(end);
    bins       = (tF-t0)/T * bins_per_T + 1;
    [post_locs, post_amps, post_hist, post_hist_t] = get_locations_from_histogram(tk_vec, ak_vec, ...
                                                 t0, tF, bins, threshold);

    % Compare the detected spikes with the real spikes
    num_sp   = length(locs);
    hit_sp   = false(num_sp, 1);
    sspp_ids = [];
    sspp_cpy = post_locs;
    delta_t  = 64*T_s;
    for ith_sp = 1 : num_sp
        t_i  = locs(ith_sp);
        inds = find(sspp_cpy > (t_i - delta_t/2) & sspp_cpy < (t_i + delta_t/2));

        if ~isempty(inds)
            hit_sp(ith_sp) = true;
            sspp_ids       = [sspp_ids; find(post_locs == sspp_cpy(inds(1)))];

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
    if ~isempty(post_locs)
        mse       = mean( (locs(hit_sp) - post_locs(sspp_ids)).^2 );
        std_dev   = sqrt(mse);
    end
    disp('')
    disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    disp(['Total number of real Diracs     : ' num2str(num_sp)])
    disp(['Total number of detected Diracs : ' num2str(length(post_locs))])
    disp(['Real Diracs detected            : ' num2str(sum(hit_sp))])
    if ~isempty(post_locs)
        disp(['MSE of locations                : ' num2str(mse)])
        disp(['Std deviation of locations      : ' num2str(std_dev)])
    end
    disp(['Dirac detection rate            : ' num2str(hit_rate) '%'])
    disp(['False positives                 : ' num2str(false_pos)])
    disp(' ')
end

%% Plots
plot_figs = true;
if plot_figs
    font_size   = 16;
    marker_size = 10;
    
    t1 = tau;
    t2 = tau+12;
    dt = t_vec(end) - t_vec(1);
    
    % Plot signal and samples
    figure
    set(gcf, 'Position', [100 100 560*3/5 420*3/5])
    stem(locs, amps, '^r', 'fill', 'MarkerSize', marker_size)
    hold on
    stem(seq_locs, seq_amps, 'pb', 'fill', 'MarkerSize', marker_size)
    hdl = legend('Original Diracs', 'Estimated Diracs', 'Location', 'SouthEast');
    set(hdl, 'FontSize', font_size)
    hdl = xlabel('Time (s)');
    set(hdl, 'FontSize', 16)
    axis([t1 t2 -.2 1.3])
    
    figure
    set(gcf, 'Position', [150 150 560*3/5 420*3/5])
    stem(t_n, y_n, '.r')
    hold on
    stem(t_n, e_n, '.k')
    hdl = legend('Noiseless samples', 'Noise');
    set(hdl, 'FontSize', font_size)
    hdl = xlabel('Time (s)');
    set(hdl, 'FontSize', 16)
    axis([t1 t2 -.2 1.3])
    
    t1 = tau;
    t2 = tau+8;
    dt = t_vec(end) - t_vec(1);
    
    % Plot retrieved locations scatter
    figure
    set(gcf, 'Position', [200 200 560*3/5 420*3/5])
    idx1_0     = min(win_idx);
    idx1_end   = max(win_idx);
    idx1_delta = idx1_end - idx1_0;
    orig_t_k = [locs(:).'; locs(:).'];
    orig_i   = repmat([idx1_0; idx1_end], 1, length(locs));
    h1 = plot(orig_i(:,1), orig_t_k(:,1), 'b');
    hold on
    for ith = 2:length(orig_i)
        plot(orig_i(:,ith), orig_t_k(:,ith), 'b');
    end
    h2 = scatter(win_idx, tk_vec, 12, 'r', 'filled');
    axis([t1/dt*idx1_delta t2/dt*idx1_delta t1 t2])
%     hdl = title(['Exp decaying recovery. Sliding window length = ' num2str(win_len1)]);
%     set(hdl, 'FontSize', 16)
    hdl = xlabel('n_i');
    set(hdl, 'FontSize', 16)
    hdl = ylabel('Time (s)');
    set(hdl, 'FontSize', 16)
    hdl = legend([h1, h2], 'True locations of Diracs', 'Detected locations', 'Location', 'SouthEast');
    set(hdl, 'FontSize', 14)
    
    figure
    set(gcf, 'Position', [200 200 560*3/5 420*3/5])
    stem(locs, max_detect*ones(length(locs),1), '.b', 'LineWidth', 1)
    hold on
    plot(seq_hist_t, seq_hist, 'r', 'LineWidth', 2)
    plot([seq_hist_t(1) seq_hist_t(end)], [threshold threshold], 'k')
    axis([t1 t2 -0.1*max_detect 1.4*max_detect])
    hdl = xlabel('Time (s)');
    set(hdl, 'FontSize', font_size)
    hdl = legend('True Diracs', 'Histogram', 'Threshold');
    set(hdl, 'FontSize', 14)
    
    % Plot histogram computed sequentially
    figure
    set(gcf, 'Position', [200 200 560*3/5 420*3/5])
    stem(locs, max_detect*ones(length(locs),1), '.b', 'LineWidth', 1)
    hold on
    plot(seq_hist_t, seq_hist, 'r', 'LineWidth', 2)
    stem(seq_hit_locs, max_detect*ones(length(seq_hit_locs),1), 'g', 'LineWidth', 1)
    if ~isempty(seq_fp_locs)
        stem(seq_fp_locs, max_detect*ones(length(seq_fp_locs),1), 'k', 'LineWidth', 1)
    end
    plot([seq_hist_t(1) seq_hist_t(end)], [threshold threshold], 'k')
    axis([t1 t2 -0.1*max_detect 1.4*max_detect])
    hdl = xlabel('Time (s)');
    set(hdl, 'FontSize', font_size)
    if ~isempty(seq_fp_locs)
        hdl = legend('Locations of true Dircas', 'Locations histogram', 'Hit locations', 'False positives', 'Histogram threshold');
    else
        hdl = legend('Locations of true Dircas', 'Locations histogram', 'Hit locations', 'Histogram threshold');
    end
    set(hdl, 'FontSize', 14)
%     hdl = title('Histogram computed sequentially');
%     set(hdl, 'FontSize', font_size)
    
    if compute_post_histogram
        % Plot histogram computed a posteriori
        figure
        set(gcf, 'Position', [150 150 1200 300])
        stem(locs, max_detect*ones(length(locs),1), '.b', 'LineWidth', 1)
        hold on
        plot(post_hist_t, post_hist, 'r', 'LineWidth', 2)
        stem(post_locs, max_detect*ones(length(post_locs),1), 'g', 'LineWidth', 1)
        plot([post_hist_t(1) post_hist_t(end)], [threshold threshold], 'k')
        axis([t1 t2 -0.1*max_detect 1.4*max_detect])
        hdl = xlabel('Time (s)');
        set(hdl, 'FontSize', font_size)
        hdl = legend('Locations of true Dircas', 'Locations histogram', 'Detected locations');
        set(hdl, 'FontSize', 14)
        hdl = title('A posteriori histogram');
        set(hdl, 'FontSize', font_size)
    end
    
end
