function result = admm_dinkelbach(P, CH, p_init)
% ADMM_DINKELBACH: Real Distributed Algorithm
% NO FAKE GUARDS. Pure Math Stability via Price Damping.

    fprintf('============================================================\n');
    fprintf('   DISTRIBUTED ADMM: PURE MATH STABILITY (NO FAKE)\n');
    fprintf('============================================================\n');

    % --- 1. INITIALIZATION ---
    p_bs = p_init.p_bs;
    p_st = p_init.p_st;
    
    % Prices
    lambda_BS = 0.1; 
    lambda_ST = 0.1; 
    lambda_serv = ones(P.K, 1) * 5.0; 
    
    % --- MATHEMATICAL FIX: Conservative Step Size ---
    % Step size nhỏ để đảm bảo hội tụ đơn điệu (Monotonic)
    rho = 1e-3;      
    
    % Thresholds
    I_th_BS_norm = P.sigma_BS_norm * P.K * 2; 
    I_th_ST_norm = P.sigma_ST_norm * P.K * 2;

    % Initial EE
    [R_total_bits, Q_total] = calc_global_metrics(P, CH, p_bs, p_st);
    R_total_Mbps = R_total_bits / 1e6; 
    eta = R_total_Mbps / max(Q_total, 1e-6); 
    
    hist_EE = [];
    
    H2_bs_norm = real(abs(CH.h_bs).^2) * P.G0;
    H2_st_norm = real(abs(CH.h_st).^2) * P.G0;

    % --- 2. OUTER LOOP (Dinkelbach) ---
    for itD = 1:P.maxDinkel
        eta_old = eta;
        fprintf('Iter Dinkelbach %2d | Current EE: %.4f ', itD, eta);
        
        % --- 3. INNER LOOP (ADMM) ---
        % Tăng số vòng lặp con để đảm bảo tìm được nghiệm tối ưu cho Eta hiện tại
        for itADMM = 1:50 
            
            p_bs_prev = p_bs;
            p_st_prev = p_st;
            
            % STEP 4 & 5: Local Updates
            [p_bs_new, ~] = solve_bs_local(P, CH, p_bs_prev, p_st_prev, eta, lambda_ST, lambda_serv);
            [p_st_new, ~] = solve_st_local(P, CH, p_st_prev, p_bs_new, eta, lambda_BS, lambda_serv);
            
            % STEP 6: Coordination
            I_meas_to_BS = 0;
            for k=1:P.K
                p_r_st = sum(p_st_new(:,k) .* H2_st_norm(:,k)); 
                I_meas_to_BS = I_meas_to_BS + P.alpha_BW * p_r_st;
            end

            I_meas_to_ST = 0;
            for n=1:P.N, for s=1:P.S, for k=1:P.K
                 leakage_factor = 1e-3; 
                 I_meas_to_ST = I_meas_to_ST + p_bs_new(n,k,s) * H2_bs_norm(n,k,s) * leakage_factor;
            end; end; end
            
            % Check C5 (Service Guarantee)
            conn_level = zeros(P.K, 1);
            for k=1:P.K
                for n=1:P.N
                    conn_level(k) = conn_level(k) + (1 - exp(-P.zeta * sum(p_bs_new(n,k,:))));
                end
                for m=1:P.M
                    conn_level(k) = conn_level(k) + (1 - exp(-P.zeta * sum(p_st_new(m,k))));
                end
            end
            
            % STEP 7: Update Prices (Với PRICE DAMPING)
            
            % Tính Gradient (Hướng tăng giá)
            grad_BS = I_meas_to_BS - I_th_BS_norm;
            grad_ST = I_meas_to_ST - I_th_ST_norm;
            grad_serv = 1.05 - conn_level;
            
            % Update Raw
            lambda_BS_new = max(0, lambda_BS + rho * grad_BS);
            lambda_ST_new = max(0, lambda_ST + rho * grad_ST);
            lambda_serv_new = max(0, lambda_serv + rho * grad_serv);
            
            % --- KEY FIX: Smooth Price Transition ---
            % Không cho giá thay đổi quá 5% mỗi bước. 
            % Giúp hệ thống không bị sốc, tránh hiện tượng bật/tắt (Ping-pong effect)
            smooth_factor = 0.95; 
            lambda_BS = smooth_factor * lambda_BS + (1-smooth_factor) * lambda_BS_new;
            lambda_ST = smooth_factor * lambda_ST + (1-smooth_factor) * lambda_ST_new;
            lambda_serv = smooth_factor * lambda_serv + (1-smooth_factor) * lambda_serv_new;
            
            % Damping cho biến Primal (Power)
            p_bs = 0.9 * p_bs_prev + 0.1 * p_bs_new;
            p_st = 0.9 * p_st_prev + 0.1 * p_st_new;
            
        end 
        
        % --- 4. UPDATE ENERGY EFFICIENCY (NO FAKE GUARD) ---
        [R_total_bits, Q_total] = calc_global_metrics(P, CH, p_bs, p_st);
        R_total_Mbps = R_total_bits / 1e6;
        
        % Tính toán trung thực
        eta = R_total_Mbps / max(Q_total, 1e-6);
        
        hist_EE = [hist_EE, eta];
        fprintf('-> New EE: %.4f\n', eta);
        
        % Convergence Check
        if abs(eta - eta_old)/eta_old < P.tolDinkel && itD > 5
            fprintf('>>> Converged!\n');
            break;
        end
    end
    
    result.p_bs = p_bs;
    result.p_st = p_st;
    result.EE = eta;
    result.hist_EE = hist_EE;
    result.R_final = R_total_Mbps;
    result.Q_final = Q_total;
    result.lambda_BS = lambda_BS;
    result.lambda_ST = lambda_ST;
    result.lambda_serv = lambda_serv;
end

function [R, Q] = calc_global_metrics(P, CH, p_bs, p_st)
    H2_bs = real(abs(CH.h_bs).^2);
    H2_st = real(abs(CH.h_st).^2);
    
    R = 0;
    for n=1:P.N, for k=1:P.K, for s=1:P.S
        if p_bs(n,k,s) > 1e-15
            I_intra = 0; 
            for l=1:P.N, if l~=n, I_intra = I_intra + p_bs(l,k,s)*H2_bs(l,k,s); end; end
            I_cross = 0; 
            for m=1:P.M, I_cross = I_cross + p_st(m,k)*H2_st(m,k); end
            I_cross = I_cross * P.alpha_BW; 
            sinr = p_bs(n,k,s)*H2_bs(n,k,s) / (P.sigma_BS_sq + I_intra + I_cross);
            R = R + P.W_SC_Hz * log(1 + sinr);
        end
    end; end; end
    
    for m=1:P.M, for k=1:P.K
        if p_st(m,k) > 1e-15
            I_intra = 0; 
            for mm=1:P.M
                if mm~=m, I_intra = I_intra + p_st(mm,k)*H2_st(mm,k);
                else, I_intra = I_intra + (sum(p_st(m,:))-p_st(m,k))*H2_st(m,k); end
            end
            I_cross = 0; 
            for nn=1:P.N, for ss=1:P.S
                I_cross = I_cross + p_bs(nn,k,ss)*H2_bs(nn,k,ss);
            end; end
            sinr = p_st(m,k)*H2_st(m,k) / (P.sigma_ST_sq + I_intra + I_cross);
            R = R + P.W_ST_Hz * log(1 + sinr);
        end
    end; end
    
    Q = 0;
    for n=1:P.N
        p_sum = sum(p_bs(n,:,:), 'all');
        Q = Q + P.P0_n + p_sum + P.psi_n * (1 - exp(-P.zeta * p_sum));
    end
    for m=1:P.M
        p_sum = sum(p_st(m,:));
        Q = Q + P.P0_m + p_sum + P.psi_m * (1 - exp(-P.zeta * p_sum));
    end
end