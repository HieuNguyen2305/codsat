function [p_bs_opt, p_st_opt, R_real, Q_real, info] = sca_subproblem(P, CH, p_bs_k, p_st_k, lambda)
% SCA_SUBPROBLEM: Solves Eq (23) with NORMALIZED units.

    N = P.N; K = P.K; S = P.S; M = P.M; zeta = P.zeta; 

    % --- NORMALIZE CHANNELS & NOISE ---
    % Scale H2 by G0 so that P * H2 is comparable to Noise=1
    H2_bs = real(abs(CH.h_bs).^2) * P.G0;
    H2_st = real(abs(CH.h_st).^2) * P.G0;
    
    % Normalized Noise (approx 1)
    sig_BS = P.sigma_BS_norm;
    sig_ST = P.sigma_ST_norm;
    
    % Penalty Weights (Scaled to match Rate in Mbps, e.g., 10-100 range)
    W_penalty = 1e3; 

    % --- 1. LINEARIZE POWER CONSUMPTION ---
    % Gradients
    Grad_Q_BS = zeros(N,K,S);
    Q_BS_val_k = 0;
    for n=1:N
        p_sum_n = sum(sum(p_bs_k(n,:,:)));
        term_exp = exp(-zeta * p_sum_n);
        Grad_Q_BS(n,:,:) = 1 + P.psi_n * zeta * term_exp; 
        Q_BS_val_k = Q_BS_val_k + P.P0_n + p_sum_n + P.psi_n*(1 - term_exp);
    end
    
    Grad_Q_ST = zeros(M,K);
    Q_ST_val_k = 0;
    for m=1:M
        p_sum_m = sum(p_st_k(m,:));
        term_exp = exp(-zeta * p_sum_m);
        Grad_Q_ST(m,:) = 1 + P.psi_m * zeta * term_exp;
        Q_ST_val_k = Q_ST_val_k + P.P0_m + p_sum_m + P.psi_m*(1 - term_exp);
    end
    Q_total_k = Q_BS_val_k + Q_ST_val_k;

    % --- 2. PREPARE RATE POINTS ---
    [I_bs_val, I_st_val] = calc_interference(P, H2_bs, H2_st, p_bs_k, p_st_k, sig_BS, sig_ST);
    mu_BS_point = log(I_bs_val); % ln domain
    mu_ST_point = log(I_st_val); % ln domain

    % --- 3. CVX SOLVER ---
    cvx_begin quiet
        cvx_solver mosek
        cvx_precision default % 'low' is risky if unscaled, 'default' is better now
        
        variables p_bs(N,K,S) p_st(M,K)
        variables mu_BS(N,K,S) mu_ST(M,K)
        variables Omega_BS(N,K,S) Omega_ST(M,K)
        
        % Slack Variables
        variables sl_c1(N,S) sl_c2(N,K) sl_c3(K) sl_c4(K) sl_c5(K)
        variables sl_rate(K)

        p_bs >= 0; p_st >= 0;
        sl_c1 >= 0; sl_c2 >= 0; sl_c3 >= 0; sl_c4 >= 0; sl_c5 >= 0; sl_rate >= 0;

        % --- OBJECTIVE FUNCTION (Normalized to Mbps) ---
        Rate_Sum_Nats = 0;
        % Use NORMALIZED Bandwidth (MHz) -> Output is Mbps
        for n=1:N, for k=1:K, for s=1:S
            Rate_Sum_Nats = Rate_Sum_Nats + P.W_SC_Norm * (Omega_BS(n,k,s) - mu_BS(n,k,s));
        end; end; end
        for m=1:M, for k=1:K
            Rate_Sum_Nats = Rate_Sum_Nats + P.W_ST_Norm * (Omega_ST(m,k) - mu_ST(m,k));
        end; end
        
        % Total Rate in Mbps
        Rate_Total_Mbps = Rate_Sum_Nats; 

        % Linearized Power (Watts)
        Q_lin = Q_total_k + sum(vec(Grad_Q_BS .* (p_bs - p_bs_k))) + ...
                            sum(vec(Grad_Q_ST .* (p_st - p_st_k)));

        Slack_Penalty = W_penalty * (sum(sl_c1(:)) + sum(sl_c2(:)) + sum(sl_c3(:)) + ...
                                     sum(sl_c4(:)) + sum(sl_c5(:)) + sum(sl_rate(:)));
        
        % Objective: Mbps - (Mbps/Watt)*Watt - Penalty
        maximize( Rate_Total_Mbps - lambda * Q_lin - Slack_Penalty - P.prox_weight*(sum_square(p_bs(:)-p_bs_k(:)) + sum_square(p_st(:)-p_st_k(:))) )
        
        % --- CONSTRAINTS ---
        % Power Budgets
        for n=1:N, sum(sum(p_bs(n,:,:))) <= P.pBS_max; end
        for m=1:M, sum(p_st(m,:)) <= P.PST_max; end

        % Association Constraints (Relaxed)
        % C1
        for n=1:N, for s=1:S
            lhs = 0;
            for k=1:K
                p0 = p_bs_k(n,k,s);
                grad = zeta * exp(-zeta * p0);
                val  = 1 - exp(-zeta * p0);
                lhs = lhs + (val + grad * (p_bs(n,k,s) - p0));
            end
            lhs <= 1 + sl_c1(n,s);
        end; end
        % C2
        for n=1:N, for k=1:K
            lhs = 0;
            for s=1:S
                p0 = p_bs_k(n,k,s);
                grad = zeta * exp(-zeta * p0);
                val  = 1 - exp(-zeta * p0);
                lhs = lhs + (val + grad * (p_bs(n,k,s) - p0));
            end
            lhs <= P.S_bar + sl_c2(n,k);
        end; end
        % C3
        for k=1:K
            lhs = 0;
            for n=1:N
                p_sum_0 = sum(p_bs_k(n,k,:));
                grad = zeta * exp(-zeta * p_sum_0);
                val  = 1 - exp(-zeta * p_sum_0);
                lhs = lhs + (val + grad * (sum(p_bs(n,k,:)) - p_sum_0));
            end
            lhs <= 1 + sl_c3(k);
        end
        % C4
        for k=1:K
            lhs = 0;
            for m=1:M
                p0 = p_st_k(m,k);
                grad = zeta * exp(-zeta * p0);
                val  = 1 - exp(-zeta * p0);
                lhs = lhs + (val + grad * (p_st(m,k) - p0));
            end
            lhs <= 1 + sl_c4(k);
        end
        % C5
        for k=1:K
            lhs = 0;
            for n=1:N
                p_sum_0 = sum(p_bs_k(n,k,:));
                grad = zeta * exp(-zeta * p_sum_0);
                val  = 1 - exp(-zeta * p_sum_0);
                lhs = lhs + (val + grad * (sum(p_bs(n,k,:)) - p_sum_0));
            end
            for m=1:M
                p0 = p_st_k(m,k);
                grad = zeta * exp(-zeta * p0);
                val  = 1 - exp(-zeta * p0);
                lhs = lhs + (val + grad * (p_st(m,k) - p0));
            end
            lhs >= 1 - sl_c5(k);
        end

        % --- RATE CONSTRAINTS (Normalized) ---
        % BS
        for n=1:N, for k=1:K, for s=1:S
            Interf = 0;
            for l=1:N
                if l~=n, Interf = Interf + sum(p_bs(l,:,s)) * H2_bs(l,k,s); end
            end
            Interf_ST = 0;
            for m_idx=1:M, Interf_ST = Interf_ST + sum(p_st(m_idx,:)) * H2_st(m_idx,k); end
            Interf = Interf + P.alpha_BW * Interf_ST;
            
            exp(mu_BS_point(n,k,s)) * (mu_BS(n,k,s) - mu_BS_point(n,k,s) + 1) >= Interf + sig_BS;
            Signal = p_bs(n,k,s)*H2_bs(n,k,s);
            -rel_entr(1, Signal + Interf + sig_BS) >= Omega_BS(n,k,s);
        end; end; end
        
        % ST
        for m=1:M, for k=1:K
            Interf = 0;
            for nn=1:N, for ss=1:S, Interf = Interf + sum(p_bs(nn,:,ss)) * H2_bs(nn,k,ss); end; end
            for mm=1:M
                if mm~=m
                    Interf = Interf + sum(p_st(mm,:)) * H2_st(mm,k);
                else
                    p_other = sum(p_st(m,:)) - p_st(m,k);
                    Interf = Interf + p_other * H2_st(m,k);
                end
            end

            exp(mu_ST_point(m,k)) * (mu_ST(m,k) - mu_ST_point(m,k) + 1) >= Interf + sig_ST;
            Signal = p_st(m,k)*H2_st(m,k);
            -rel_entr(1, Signal + Interf + sig_ST) >= Omega_ST(m,k);
        end; end
        
        % QoS Threshold (Mbps)
        for k=1:K
            R_k_Mbps = 0;
            for n=1:N, for s=1:S, R_k_Mbps = R_k_Mbps + P.W_SC_Norm*(Omega_BS(n,k,s) - mu_BS(n,k,s)); end; end
            for m=1:M, R_k_Mbps = R_k_Mbps + P.W_ST_Norm*(Omega_ST(m,k) - mu_ST(m,k)); end
            
            R_k_Mbps >= P.R_threshold_Mbps - sl_rate(k);
        end

    cvx_end
    
    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
        p_bs_opt = full(real(p_bs)); 
        p_st_opt = full(real(p_st));
        info.status = 'Success';
    else
        p_bs_opt = p_bs_k; 
        p_st_opt = p_st_k;
        info.status = 'Failed';
    end
    
    [R_real, Q_real] = calc_real_metrics_smooth(P, H2_bs, H2_st, p_bs_opt, p_st_opt, sig_BS, sig_ST);
end

function [I_bs, I_st] = calc_interference(P, H2_bs, H2_st, p_bs, p_st, sig_BS, sig_ST)
    N=P.N; K=P.K; S=P.S; M=P.M;
    I_bs = zeros(N,K,S);
    for n=1:N, for k=1:K, for s=1:S
        val = 0;
        for l=1:N, if l~=n, val = val + sum(p_bs(l,:,s)) * H2_bs(l,k,s); end; end
        st_interf = 0;
        for m=1:M, st_interf = st_interf + sum(p_st(m,:)) * H2_st(m,k); end
        val = val + P.alpha_BW * st_interf;
        I_bs(n,k,s) = val + sig_BS;
    end; end; end
    
    I_st = zeros(M,K);
    for m=1:M, for k=1:K
        val = 0;
        for nn=1:N, for ss=1:S, val = val + sum(p_bs(nn,:,ss)) * H2_bs(nn,k,ss); end; end
        for mm=1:M
             if mm~=m, val = val + sum(p_st(mm,:)) * H2_st(mm,k);
             else, val = val + (sum(p_st(m,:)) - p_st(m,k)) * H2_st(m,k); end
        end
        I_st(m,k) = val + sig_ST;
    end; end
end

function [R, Q] = calc_real_metrics_smooth(P, H2_bs, H2_st, p_bs, p_st, sigma_BS, sigma_ST)
    R = 0;
    for n=1:P.N, for k=1:P.K, for s=1:P.S
        if p_bs(n,k,s) > 1e-15
            I = 0; 
            for l=1:P.N, if l~=n, I = I + sum(p_bs(l,:,s))*H2_bs(l,k,s); end; end
            st_int = 0;
            for m=1:P.M, st_int = st_int + sum(p_st(m,:))*H2_st(m,k); end
            I = I + P.alpha_BW * st_int;
            
            % Mbps Calculation (Using W_Norm)
            R = R + P.W_SC_Norm * log(1 + p_bs(n,k,s)*H2_bs(n,k,s)/(sigma_BS + I));
        end
    end; end; end
    
    for m=1:P.M, for k=1:P.K
        if p_st(m,k) > 1e-15
            I = 0; 
            for nn=1:P.N, for ss=1:P.S, I = I + sum(p_bs(nn,:,ss))*H2_bs(nn,k,ss); end; end
            for mm=1:P.M
                if mm~=m, I = I + sum(p_st(mm,:)) * H2_st(mm,k);
                else, I = I + (sum(p_st(m,:)) - p_st(m,k)) * H2_st(m,k); end
            end
            
            R = R + P.W_ST_Norm * log(1 + p_st(m,k)*H2_st(m,k)/(sigma_ST + I));
        end
    end; end
    % No need to multiply by RateScale anymore, W_Norm handles it.
    
    Q_BS = 0;
    for n=1:P.N
        p_sum = sum(sum(p_bs(n,:,:)));
        Q_BS = Q_BS + P.P0_n + p_sum + P.psi_n * (1 - exp(-P.zeta * p_sum));
    end
    Q_ST = 0;
    for m=1:P.M
        p_sum = sum(p_st(m,:));
        Q_ST = Q_ST + P.P0_m + p_sum + P.psi_m * (1 - exp(-P.zeta * p_sum));
    end
    Q = Q_BS + Q_ST;
end