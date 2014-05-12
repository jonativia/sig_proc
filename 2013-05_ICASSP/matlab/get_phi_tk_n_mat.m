function phi_mat = get_phi_tk_n_mat(phi, t_phi, t_k, n_vec, T, T_s)

phi_mat = zeros(length(n_vec), length(t_k));
for ith_n = 1:length(n_vec)
    n = n_vec(ith_n);
    [~,idx_p,idx_t] = intersect(round((t_phi+n*T)/T_s), round(t_k/T_s));
    if ~isempty(idx_t)
        phi_mat(ith_n,idx_t) = phi(idx_p);
    end
end

