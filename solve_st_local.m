function [p_st_new, slack_val] = solve_st_local(P, CH, p_st_old, p_bs_fixed, eta, lambda_bs, lambda_serv)
% SOLVE_ST_LOCAL: Normalized + Service Reward

    M = P.M; K = P.K;
    
    H2_st = real(abs(CH.h_st).^2) * P.G0;
    H2_bs = real(abs(CH.h_bs).^2) * P.G0;
    sig_ST = P.sigma_ST_norm;
    W_SLACK = 1e3; 

    I_from_BS = zeros(M,K);
    for k=1:K
        val_bs = 0;
        for n=1:P.N, for s=1:P.S
            val_bs = val_bs + p_bs_fixed(n,k,s) * H2_bs(n,k,s);
        end; end
        I_from_BS(:,k) = val_bs; 
    end

    I_st_val_old = zeros(M,K);
    for m=1:M, for k=1:K
        i_intra = 0;
        for mm=1:M
            if mm ~= m, i_intra = i_intra + p_st_old(mm,k) * H2_st(mm,k);
            else, p_other = sum(p_st_old(m,:)) - p_st_old(m,k); i_intra = i_intra + p_other * H2_st(m,k); end
        end
        I_st_val_old(m,k) = i_intra + I_from_BS(m,k);
    end; end
    
    mu_old = log(I_st_val_old + sig_ST);

    grad_Q_st = zeros(M,1);
    for m=1:M
        p_sum_m = sum(p_st_old(m,:));
        grad_Q_st(m) = 1 + P.zeta * P.psi_m * exp(-P.zeta * p_sum_m);
    end

    val_f_st  = 1 - exp(-P.zeta * p_st_old); 
    grad_f_st = P.zeta * exp(-P.zeta * p_st_old);
    H_ST2BS_norm = 1e-13 * P.G0;

    % --- CVX ---
    cvx_begin
        cvx_solver mosek
        cvx_precision default
        cvx_quiet(true);
        
        variable p_st(M,K) nonnegative
        variable omega_st(M,K) nonnegative
        variable mu_st(M,K)
        
        variable s_c4(K) nonnegative
        variable s_rate(M,K) nonnegative

        % --- OBJECTIVE ---
        obj_Rate = sum(sum(omega_st - mu_st)) * P.W_ST_Norm; 
        
        obj_Energy = 0;
        for m=1:M
            obj_Energy = obj_Energy + grad_Q_st(m) * sum(p_st(m,:));
        end
        
        scaling_factor = P.W_SC_Hz / P.W_ST_Hz;
        obj_Interf_Cost = lambda_bs * scaling_factor * sum(sum(p_st)) * H_ST2BS_norm;
        
        % 4. NEW: Service Reward (Dual Decomposition)
        obj_Reward = 0;
        for k=1:K
            term_k = 0;
            for m=1:M
                % Grad of Beta_mk = grad_f_st * p_st
                term_k = term_k + grad_f_st(m,k) * p_st(m,k);
            end
            obj_Reward = obj_Reward + lambda_serv(k) * term_k;
        end
        
        obj_Prox = P.prox_weight * sum_square(p_st(:) - p_st_old(:));
        obj_Penalty = W_SLACK * (sum(s_c4) + sum(s_rate(:)));

        maximize( obj_Rate - eta * obj_Energy - obj_Interf_Cost + obj_Reward - obj_Prox - obj_Penalty )

        subject to
            for m=1:M, sum(p_st(m,:)) <= P.PST_max; end
            
            for m=1:M, for k=1:K
                Signal = p_st(m,k) * H2_st(m,k);
                Interf_Intra = 0;
                for mm=1:M
                    if mm ~= m, Interf_Intra = Interf_Intra + p_st(mm,k) * H2_st(mm,k);
                    else, p_other = sum(p_st(m,:)) - p_st(m,k); Interf_Intra = Interf_Intra + p_other * H2_st(m,k); end
                end
                Interf_Total = Interf_Intra + I_from_BS(m,k) + sig_ST;

                log(Signal + Interf_Total) >= omega_st(m,k);
                I0 = I_st_val_old(m,k) + sig_ST;
                log(I0) + (Interf_Total - I0)/I0 <= mu_st(m,k) + s_rate(m,k);
            end; end
            
            % C4
            for k=1:K
                lhs = 0; for m=1:M, lhs = lhs + val_f_st(m,k) + grad_f_st(m,k)*(p_st(m,k)-p_st_old(m,k)); end; lhs <= 1 + s_c4(k);
            end
    cvx_end
    
    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
        p_st_new = full(p_st);
        slack_val = sum(s_c4) + sum(s_rate(:));
    else
        p_st_new = p_st_old;
        slack_val = 1e6;
    end
end