
clc, clear all, close all

SNR_vec    = [5, 10, 15, 20];
% SNR_vec = [10];
num_snr    = length(SNR_vec);
iterations = 100;
det_rates  = zeros(num_snr, iterations);
fps        = zeros(num_snr, iterations);
std_devs   = zeros(num_snr, iterations);
valid_it   = true(num_snr, iterations);

% Signal's characteristics
K   = 5;
tau = K^2 / 8;
N   = 2 * K^2;
T   = 1 / 16;

% Parameters
% P     = 2*K-1;
P = 22;
L = P+1;
TTs   = 64;
T_s   = T / TTs;
thres = 0.25;

load('thousand_diracs.mat', 'locs', 'amps')

n0 = floor((locs(1)-tau)/T);
nF = ceil((locs(end)+tau)/T+L);

t_vec = ( n0*T : T_s : nF*T )';
n_vec = ( n0 : nF )';
t_n   = n_vec*T;

% Generate the continuous-time signal
i_locs    = round((locs- t_vec(1))/T_s) + 1;
x         = zeros(size(t_vec));
x(i_locs) = amps;

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

for ith_snr = 1 : num_snr
    SNR = SNR_vec(ith_snr);
    
    disp(['SNR = ' num2str(SNR) ' dB'])
    
    for it = 1 : iterations
    
        disp([' - iteration ' num2str(it)])
    
        % Add noise
        P_y   = y_n' * y_n / length(y_n);
        e_n   = randn(size(y_n));
        P_e   = e_n' * e_n / length(e_n);
        A_e   = sqrt(10^(-SNR/10) * P_y / P_e);
        e_n   = A_e * e_n;
        sigma = sqrt((e_n' * e_n) / length(e_n));

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
            nA = n_vec(ith+1);
            nB = n_vec(ith) + N;

            true_locs = locs(locs>t0 & locs<tF);
            true_amps = amps(locs>t0 & locs<tF);

            % Obtain samples that correspond to this time interval
            yi_n = y_n(ith+1:ith+N) + e_n(ith+1:ith+N);


            % Apply mask and compute moments
            yi_n = yi_n .* mask;
            si_m = c_m_n * yi_n;

            % Oracle to get number of Diracs within the time interval or fixed K
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
            t_k(abs(a_k)<0.25) = [];
            a_k(abs(a_k)<0.25) = [];

            if ~isempty(t_k)
                tk_vec  = [tk_vec; t_k];
                ak_vec  = [ak_vec; a_k];
                win_idx = [win_idx; ones(length(t_k),1)*ith];
            end

            % Retrieve Dirac in ]t0,t0+T] from histogram
            [locats, amplis, histo, t_bins] = get_locations_from_histogram(tk_vec, ak_vec, ...
                                                                           t0-T/bins_per_T, ...
                                                                           t0+T, ...
                                                                           bins_per_T+2, ...
                                                                           threshold);
            seq_locs   = [seq_locs; locats];
            seq_amps   = [seq_amps; amplis];
            seq_hist   = [seq_hist histo(2:end-1)];
            seq_hist_t = [seq_hist_t t_bins(2:end-1)];
        end

        % Retrieve Diracs of last interval
        t_last = t0;
        for t0 = t0 : T : t_n(end)
            [locats, amplis, histo, t_bins] = get_locations_from_histogram(tk_vec, ak_vec, ...
                                                                           t0-T/bins_per_T, ...
                                                                           t0+T, ...
                                                                           bins_per_T+2, ...
                                                                           threshold);
            seq_locs   = [seq_locs; locats];
            seq_amps   = [seq_amps; amplis];
            seq_hist   = [seq_hist histo(2:end-1)];
            seq_hist_t = [seq_hist_t t_bins(2:end-1)];
        end

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

        % Check consistency of detected diracs
        if sum(hit_sp) ~= length(seq_hit_locs)
            valid_it(ith_snr,it) = false;
            continue;
        end
        
        % Accuracy of detected spikes
        det_rates(ith_snr,it) = sum(hit_sp) / length(hit_sp) * 100;
        fps(ith_snr,it)       = length(seq_fp_locs);
        std_devs(ith_snr,it)  = sqrt( mean( (locs(hit_sp) - seq_hit_locs).^2 ) );


    end
end

%%
save(['results/perfs_noisy_' num2str(iterations) 'its.mat'], 'SNR_vec', 'det_rates', 'fps', 'std_devs', 'valid_it')

