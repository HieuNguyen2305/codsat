function init_result = find_init_point(P, CH)
% FIND_INIT_POINT: Relaxed Feasibility Search (Normalized & Paper Compliant)
% Solves a feasibility problem by minimizing slack variables to find a 
% valid starting point for the SCA algorithm.

    fprintf('--- [Init Phase] Solving Relaxed Feasibility Problem ---\n');
    
    N = P.N; K = P.K; S = P.S; M = P.M; 
    zeta = P.zeta; 

    % --- 1. NORMALIZE CHANNELS & NOISE ---
    % Scale H2 by G0 so that Channel Gain is comparable to Noise (~1)
    % This prevents numerical instability (Unbounded/Infeasible status)
    H2_bs = real(abs(CH.h_bs).^2) * P.G0;
    H2_st = real(abs(CH.h_st).^2) * P.G0;
    
    % Normalized Noise (approx 1)
    sig_BS = P.sigma_BS_norm;
    sig_ST = P.sigma_ST_norm;

    % Initialize Variables randomly (Small values to start safe)
    p_bs_k = rand(N,K,S) * (P.pBS_max / (K*S)) * 0.1;
    p_st_k = rand(M,K) * (P.PST_max / K) * 0.5;
    
    % Auxiliary variables for Rate
    mu_BS_k = zeros(N,K,S); 
    mu_ST_k = zeros(M,K);
    
    % Loop to iteratively reduce violations
    for iter = 1:35 
        
        % Update Interference & Taylor Points
        [I_bs_val, I_st_val] = calc_interference(P, H2_bs, H2_st, p_bs_k, p_st_k, sig_BS, sig_ST);
        
        % Avoid log(0)
        I_bs_val = max(I_bs_val, 1e-10);
        I_st_val = max(I_st_val, 1e-10);

        % Taylor expansion points (Natural Log - ln)
        mu_BS_point = log(I_bs_val);
        mu_ST_point = log(I_st_val);

        cvx_begin quiet
            cvx_solver mosek
            cvx_precision default
            
            variables p_bs(N,K,S) p_st(M,K)
            variables mu_BS(N,K,S) mu_ST(M,K)
            variables Omega_BS(N,K,S) Omega_ST(M,K)
            
            % --- SLACK VARIABLES (Relax everything) ---
            % These allow the solver to violate constraints temporarily
            variables slack_c1(N,S) slack_c2(N,K) slack_c3(K) slack_c4(K) slack_c5(K)
            variables slack_rate(K)
            
            % Non-negative Powers & Slacks
            p_bs >= 0; p_st >= 0;
            slack_c1 >= 0; slack_c2 >= 0; slack_c3 >= 0; slack_c4 >= 0; slack_c5 >= 0;
            slack_rate >= 0; 

            % --- OBJECTIVE: Minimize Violations ---
            Obj_Topo = sum(slack_c1(:)) + sum(slack_c2(:)) + sum(slack_c3) + sum(slack_c4) + sum(slack_c5);
            Obj_Rate = sum(slack_rate);
            
            % Proximal term to keep variable changes smooth and stable
            Obj_Prox = 0.1 * (sum_square(p_bs(:) - p_bs_k(:)) + sum_square(p_st(:) - p_st_k(:)));
            
            % Weighting: Priority on Topology > Rate > Stability
            minimize( 1e6*Obj_Topo + 10*Obj_Rate + Obj_Prox )

            % --- POWER BUDGETS (Hard Constraints) ---
            for n=1:N, sum(sum(p_bs(n,:,:))) <= P.pBS_max; end
            for m=1:M, sum(p_st(m,:)) <= P.PST_max; end

            % --- ASSOCIATION CONSTRAINTS (Relaxed C1-C5) ---
            % C1: One UE per SC per BS
            for n=1:N, for s=1:S
                lhs = 0;
                for k=1:K
                    p0 = p_bs_k(n,k,s);
                    grad = zeta * exp(-zeta * p0);
                    val  = 1 - exp(-zeta * p0);
                    lhs = lhs + (val + grad * (p_bs(n,k,s) - p0));
                end
                lhs <= 1 + slack_c1(n,s); 
            end; end

            % C2: Max S_bar SCs per UE
            for n=1:N, for k=1:K
                lhs = 0;
                for s=1:S
                    p0 = p_bs_k(n,k,s);
                    grad = zeta * exp(-zeta * p0);
                    val  = 1 - exp(-zeta * p0);
                    lhs = lhs + (val + grad * (p_bs(n,k,s) - p0));
                end
                lhs <= P.S_bar + slack_c2(n,k);
            end; end

            % C3: One BS per UE
            for k=1:K
                lhs = 0;
                for n=1:N
                    p_sum_0 = sum(p_bs_k(n,k,:));
                    grad = zeta * exp(-zeta * p_sum_0);
                    val  = 1 - exp(-zeta * p_sum_0);
                    p_sum_var = sum(p_bs(n,k,:));
                    lhs = lhs + (val + grad * (p_sum_var - p_sum_0));
                end
                lhs <= 1 + slack_c3(k);
            end

            % C4: One ST per UE
            for k=1:K
                lhs = 0;
                for m=1:M
                    p0 = p_st_k(m,k);
                    grad = zeta * exp(-zeta * p0);
                    val  = 1 - exp(-zeta * p0);
                    lhs = lhs + (val + grad * (p_st(m,k) - p0));
                end
                lhs <= 1 + slack_c4(k);
            end

            % C5: Connectivity >= 1
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
                lhs >= 1 - slack_c5(k); 
            end

            % --- RATE CONSTRAINTS (Normalized & Relaxed) ---
            
            % 1. Terrestrial Users (BS)
            for n=1:N, for k=1:K, for s=1:S
                % Interference (Summing over interfering BSs)
                Interf = 0;
                for l=1:N
                    if l~=n
                        Interf = Interf + sum(p_bs(l,:,s)) * H2_bs(l,k,s);
                    end
                end
                % ST Interference (Eq 15: Includes alpha_BW scaling)
                St_Interf = 0;
                for m=1:M
                    St_Interf = St_Interf + sum(p_st(m,:)) * H2_st(m,k);
                end
                Interf = Interf + P.alpha_BW * St_Interf;
                
                % Upper Bound on Interference
                rhs = exp(mu_BS_point(n,k,s)) * (mu_BS(n,k,s) - mu_BS_point(n,k,s) + 1);
                Interf + sig_BS <= rhs; 
                
                % Lower Bound on Rate (using rel_entr for ln)
                Signal = p_bs(n,k,s)*H2_bs(n,k,s);
                -rel_entr(1, Signal + Interf + sig_BS) >= Omega_BS(n,k,s);
            end; end; end
            
            % 2. Satellite Users (ST)
            for m=1:M, for k=1:K
                % BS Interference (Eq 144: NO alpha_BW scaling)
                Interf = 0;
                for nn=1:N, for ss=1:S
                    Interf = Interf + sum(p_bs(nn,:,ss)) * H2_bs(nn,k,ss);
                end; end
                
                % ST Interference (Inter-beam/Inter-satellite)
                for mm=1:M
                    if mm~=m
                        Interf = Interf + sum(p_st(mm,:)) * H2_st(mm,k);
                    else
                         p_other = sum(p_st(m,:)) - p_st(m,k);
                         Interf = Interf + p_other * H2_st(m,k);
                    end
                end
                
                rhs = exp(mu_ST_point(m,k)) * (mu_ST(m,k) - mu_ST_point(m,k) + 1);
                Interf + sig_ST <= rhs;
                
                Signal = p_st(m,k)*H2_st(m,k);
                -rel_entr(1, Signal + Interf + sig_ST) >= Omega_ST(m,k);
            end; end

            % Rate Threshold (Using Normalized Bandwidth -> Mbps)
            for k=1:K
                R_k_Mbps = 0;
                for n=1:N, for s=1:S
                    R_k_Mbps = R_k_Mbps + P.W_SC_Norm * (Omega_BS(n,k,s) - mu_BS(n,k,s));
                end; end
                for m=1:M
                    R_k_Mbps = R_k_Mbps + P.W_ST_Norm * (Omega_ST(m,k) - mu_ST(m,k));
                end
                
                % Relaxed QoS Constraint
                R_k_Mbps >= P.R_threshold_Mbps - slack_rate(k);
            end
            
        cvx_end
        
        % Check Status
        if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
            p_bs_k = full(real(p_bs)); 
            p_st_k = full(real(p_st));
            
            viol_topo = sum(slack_c1(:)) + sum(slack_c2(:)) + sum(slack_c3) + sum(slack_c4) + sum(slack_c5);
            viol_rate = sum(slack_rate);
            
            fprintf('   [Init SCA %d] Viol Topo: %.2e | Viol Rate: %.2e\n', iter, viol_topo, viol_rate);
            
            % Feasible if violations are negligible
            if viol_topo < 1e-3 && viol_rate < 0.1
                fprintf('   >>> Feasible Point Found!\n');
                break; 
            end
        else
            fprintf('   [Init SCA %d] CVX Status: %s. Perturbing...\n', iter, cvx_status);
            p_bs_k = max(0, p_bs_k + 0.1 * rand(size(p_bs_k)) * (P.pBS_max/10));
            p_st_k = max(0, p_st_k + 0.1 * rand(size(p_st_k)) * (P.PST_max/10));
        end
    end
    
    init_result.p_bs = p_bs_k;
    init_result.p_st = p_st_k;
end

function [I_bs, I_st] = calc_interference(P, H2_bs, H2_st, p_bs, p_st, sig_BS, sig_ST)
    % Consistent Interference Calculation (Normalized)
    % Matches sca_subproblem.m and paper Eq (15) & (144)
    N=P.N; K=P.K; S=P.S; M=P.M;
    I_bs = zeros(N,K,S);
    
    % Terrestrial Interference Calculation
    for n=1:N, for k=1:K, for s=1:S
        val = 0;
        for l=1:N
            if l~=n
                % Intra-tier BS interference
                val = val + sum(p_bs(l,:,s)) * H2_bs(l,k,s);
            end
        end
        % Inter-tier ST interference (Scaled by alpha_BW) - Eq (15)
        st_interf = 0;
        for m=1:M, st_interf = st_interf + sum(p_st(m,:)) * H2_st(m,k); end
        
        I_bs(n,k,s) = val + P.alpha_BW * st_interf + sig_BS;
    end; end; end
    
    % Satellite Interference Calculation
    I_st = zeros(M,K);
    for m=1:M, for k=1:K
        val = 0;
        % Inter-tier BS interference (NO scaling) - Eq (144)
        for nn=1:N, for ss=1:S
             val = val + sum(p_bs(nn,:,ss)) * H2_bs(nn,k,ss);
        end; end
        
        % Intra-tier ST interference
        for mm=1:M
             if mm~=m
                 val = val + sum(p_st(mm,:)) * H2_st(mm,k);
             else
                 val = val + (sum(p_st(m,:)) - p_st(m,k)) * H2_st(m,k);
             end
        end
        I_st(m,k) = val + sig_ST;
    end; end
end