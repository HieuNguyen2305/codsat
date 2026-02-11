function [p_bs_new, slack_val] = solve_bs_local(P, CH, p_bs_old, p_st_fixed, eta, lambda_st, lambda_serv)
% SOLVE_BS_LOCAL: Normalized + Service Reward
% Input: lambda_serv (Vector Kx1) - Giá thưởng cho việc kết nối User

    N = P.N; K = P.K; S = P.S;
    
    % Normalize Inputs
    H2_bs = real(abs(CH.h_bs).^2) * P.G0; 
    H2_st = real(abs(CH.h_st).^2) * P.G0; 
    sig_BS = P.sigma_BS_norm;
    W_SLACK = 1e3; 

    % Interference from ST
    I_from_ST = zeros(N,K,S);
    for k=1:K
        i_st_val = sum(p_st_fixed(:,k) .* H2_st(:,k));
        I_from_ST(:,k,:) = P.alpha_BW * i_st_val; 
    end

    % SCA Prep
    I_bs_val_old = zeros(N,K,S);
    for n=1:N, for k=1:K, for s=1:S
        inter_intra = 0;
        for l=1:N, if l~=n, inter_intra = inter_intra + p_bs_old(l,k,s) * H2_bs(l,k,s); end; end
        I_bs_val_old(n,k,s) = inter_intra + I_from_ST(n,k,s);
    end; end; end
    mu_old = log(I_bs_val_old + sig_BS);

    % Gradients
    grad_Q_bs = zeros(N,1);
    for n=1:N
        p_sum_n = sum(sum(p_bs_old(n,:,:)));
        grad_Q_bs(n) = 1 + P.zeta * P.psi_n * exp(-P.zeta * p_sum_n);
    end

    val_f_bs  = 1 - exp(-P.zeta * p_bs_old); 
    grad_f_bs = P.zeta * exp(-P.zeta * p_bs_old);
    H_BS2ST_norm = 1e-13 * P.G0; 

    % --- CVX ---
    cvx_begin
        cvx_solver mosek
        cvx_precision default
        cvx_quiet(true);
        
        variable p_bs(N,K,S) nonnegative
        variable omega_bs(N,K,S) nonnegative 
        variable mu_bs(N,K,S)                
        
        variable s_c1(N,S) nonnegative
        variable s_c2(N,K) nonnegative
        variable s_c3(K)   nonnegative
        variable s_rate(N,K,S) nonnegative

        % --- OBJECTIVE ---
        
        % 1. Rate (Mbps)
        obj_Rate = sum(sum(sum(omega_bs - mu_bs))) * P.W_SC_Norm; 
        
        % 2. Energy
        obj_Energy = 0;
        for n=1:N
            obj_Energy = obj_Energy + grad_Q_bs(n) * sum(sum(p_bs(n,:,:)));
        end
        
        % 3. Interf Cost
        obj_Interf_Cost = lambda_st * sum(sum(sum(p_bs))) * H_BS2ST_norm;
        
        % 4. NEW: Service Reward (Dual Decomposition for C5)
        % Khuyến khích BS kết nối với User k bằng cách thưởng lambda_serv(k)
        obj_Reward = 0;
        for k=1:K
            term_k = 0;
            for n=1:N
                % Linearize: Reward * grad_assoc * p_bs
                % Alpha_nk = 1 - exp(-zeta * sum(p))
                p_sum_old_nk = sum(p_bs_old(n,k,:));
                grad_assoc = P.zeta * exp(-P.zeta * p_sum_old_nk);
                
                term_k = term_k + grad_assoc * sum(p_bs(n,k,:));
            end
            obj_Reward = obj_Reward + lambda_serv(k) * term_k;
        end
        
        obj_Prox = P.prox_weight * sum_square(p_bs(:) - p_bs_old(:));
        obj_Penalty = W_SLACK * (sum(s_c1(:)) + sum(s_c2(:)) + sum(s_c3(:)) + sum(s_rate(:)));

        % Maximize with Reward
        maximize( obj_Rate - eta * obj_Energy - obj_Interf_Cost + obj_Reward - obj_Prox - obj_Penalty )

        subject to
            for n=1:N, sum(sum(p_bs(n,:,:))) <= P.pBS_max; end

            % Rate SCA
            for n=1:N, for k=1:K, for s=1:S
                Signal = p_bs(n,k,s) * H2_bs(n,k,s);
                Interf = 0;
                for l=1:N, if l~=n, Interf = Interf + p_bs(l,k,s) * H2_bs(l,k,s); end; end
                Interf_Total = Interf + I_from_ST(n,k,s) + sig_BS;

                log(Signal + Interf_Total) >= omega_bs(n,k,s);
                I0 = I_bs_val_old(n,k,s) + sig_BS;
                log(I0) + (Interf_Total - I0)/I0 <= mu_bs(n,k,s) + s_rate(n,k,s);
            end; end; end

            % Constraints C1, C2, C3 (Relaxed)
            for n=1:N, for s=1:S
                lhs = 0; for k=1:K, lhs = lhs + val_f_bs(n,k,s) + grad_f_bs(n,k,s)*(p_bs(n,k,s)-p_bs_old(n,k,s)); end; lhs <= 1+s_c1(n,s);
            end; end
            for n=1:N, for k=1:K
                lhs = 0; for s=1:S, lhs = lhs + val_f_bs(n,k,s) + grad_f_bs(n,k,s)*(p_bs(n,k,s)-p_bs_old(n,k,s)); end; lhs <= P.S_bar+s_c2(n,k);
            end; end
            for k=1:K
                lhs = 0; for n=1:N
                    p_sum_old = sum(p_bs_old(n,k,:));
                    grad_sum = P.zeta * exp(-P.zeta * p_sum_old);
                    val_sum  = 1 - exp(-P.zeta * p_sum_old);
                    lhs = lhs + val_sum + grad_sum*(sum(p_bs(n,k,:))-p_sum_old);
                end; lhs <= 1+s_c3(k);
            end
    cvx_end

    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
        p_bs_new = full(p_bs);
        slack_val = sum(s_c1(:)) + sum(s_c2(:)) + sum(s_c3(:));
    else
        p_bs_new = p_bs_old;
        slack_val = 1e6; 
    end
end