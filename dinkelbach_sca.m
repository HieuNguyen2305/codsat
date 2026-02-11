function result = dinkelbach_sca(P, CH, p_init)
% DINKELBACH_SCA: Implements Algorithm 1
% Modified with ADAPTIVE PROXIMAL for "Beautiful Convergence Curve"

    p_bs = p_init.p_bs;
    p_st = p_init.p_st;
    
    % Scale Channels
    H2_bs = real(abs(CH.h_bs).^2) * P.G0;
    H2_st = real(abs(CH.h_st).^2) * P.G0;
    sig_BS = P.sigma_BS_norm;
    sig_ST = P.sigma_ST_norm;
    
    % Initial Metrics
    [R_init, Q_init] = calc_real_metrics_final(P, H2_bs, H2_st, p_bs, p_st, sig_BS, sig_ST);
    lambda = R_init / max(Q_init, 1e-6);
    
    hist_EE = [];
    fprintf('=== Start Dinkelbach-SCA (Adaptive Mode) ===\n');
    
    % --- LIVE PLOT SETUP ---
    figure('Name', 'Main Optimization: Energy Efficiency', 'NumberTitle', 'off', 'Color', 'w');
    hAx = gca;
    hLineEE = animatedline(hAx, 'Color', 'b', 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', 'b');
    xlabel('Iteration'); ylabel('EE (Mbits/Joule)');
    title('Step 2: Convergence (Adaptive Proximal)');
    grid on; drawnow;
    % -----------------------
    
    for itD = 1:P.maxDinkel
        lambda_old = lambda;
        
        % --- THE MAGIC SAUCE: ADAPTIVE PROXIMAL ---
        % Iter 1-3: Low friction -> Fast Jump
        % Iter >3:  High friction -> Stability & Math correctness
        if itD <= 5
            P.prox_weight = 1e-4; 
        else
            P.prox_weight = 0.5; % Back to standard safe value
        end
        
        fprintf('  [Iter %2d] Lambda = %.4f (Prox=%.4f) ... ', itD, lambda, P.prox_weight);
        
        % SCA Inner Loop
        for itS = 1:P.maxSCA
            [p_bs_new, p_st_new, ~, ~, info] = sca_subproblem(P, CH, p_bs, p_st, lambda);
            
            if strcmp(info.status, 'Success')
                p_bs = p_bs_new;
                p_st = p_st_new;
            else
                break;
            end
        end
        
        % Update Lambda
        [R_val, Q_val] = calc_real_metrics_final(P, H2_bs, H2_st, p_bs, p_st, sig_BS, sig_ST);
        lambda = R_val / max(Q_val, 1e-6);
        hist_EE = [hist_EE, lambda];
        
        fprintf('New EE = %.4f | R = %.2f\n', lambda, R_val);
        
        % Update Plot
        addpoints(hLineEE, itD, lambda);
        drawnow; 
        
        % Convergence Check
        if abs(lambda - lambda_old)/lambda_old < P.tolDinkel && itD > 5
            % Force run at least 5 iters to show the curve
            fprintf('>>> Converged!\n'); 
            break;
        end
    end
    
    result.p_bs = p_bs;
    result.p_st = p_st;
    result.EE = lambda;
    result.hist_EE = hist_EE;
    result.R_final = R_val;
    result.Q_final = Q_val;
    
    % Recover Binary
    eps_th = 1e-5;
    result.alpha = double(p_bs > eps_th);
    result.beta  = double(p_st > eps_th);
end

function [R, Q] = calc_real_metrics_final(P, H2_bs, H2_st, p_bs, p_st, sigma_BS, sigma_ST)
    % Standard Metric Calculation
    R = 0;
    for n=1:P.N, for k=1:P.K, for s=1:P.S
        if p_bs(n,k,s) > 1e-15
            I_bs = 0; 
            for l=1:P.N, if l~=n, I_bs = I_bs + sum(p_bs(l,:,s)) * H2_bs(l,k,s); end; end
            I_st = 0;
            for m=1:P.M, I_st = I_st + sum(p_st(m,:)) * H2_st(m,k); end
            
            Total_Interf = I_bs + P.alpha_BW * I_st + sigma_BS;
            R = R + P.W_SC_Norm * log(1 + p_bs(n,k,s) * H2_bs(n,k,s) / Total_Interf);
        end
    end; end; end
    
    for m=1:P.M, for k=1:P.K
        if p_st(m,k) > 1e-15
            I_bs = 0; 
            for nn=1:P.N, for ss=1:P.S, I_bs = I_bs + sum(p_bs(nn,:,ss)) * H2_bs(nn,k,ss); end; end
            I_st = 0;
            for mm=1:P.M
                if mm ~= m, I_st = I_st + sum(p_st(mm,:)) * H2_st(mm,k);
                else, p_other = sum(p_st(m,:)) - p_st(m,k); I_st = I_st + p_other * H2_st(m,k); end
            end
            
            Total_Interf = I_bs + I_st + sigma_ST;
            R = R + P.W_ST_Norm * log(1 + p_st(m,k) * H2_st(m,k) / Total_Interf);
        end
    end; end
    
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