function CH = gen_channels(P)
% GEN_CHANNELS: Generates channel coefficients based on Section II.A
% Inputs: P (System Parameters)
% Outputs: CH struct containing h_bs (Terrestrial) and h_st (Satellite)

    N = P.N; K = P.K; S = P.S; M = P.M;
    
    % Initialize Matrices
    CH.h_bs = zeros(N, K, S);
    CH.h_st = zeros(M, K);

    %% 1. Terrestrial Channel (Eq 1 - Rician Fading) 
    kappa_TN = 5; % Rician Factor
    
    for n = 1:N
        for k = 1:K
            % Distance Calculation
            d_nk = norm(P.bsPos(n, :) - P.userPos(k, :));
            d_nk = max(d_nk, 10); % Avoid singularity
            
            % Pathloss (Standard UMi model often used with Eq 1)
            PL_dB = 28 + 22*log10(d_nk) + 20*log10(P.fc/1e9);
            PL_lin = 10^(-PL_dB/10) * P.G_BS * P.G_UE;
            
            for s = 1:S
                % Small-scale fading
                h_los = exp(1i * 2 * pi * rand); 
                h_nlos = (randn + 1i * randn) / sqrt(2);
                
                % Eq (1): Rician Fading
                CH.h_bs(n,k,s) = sqrt(PL_lin) * (...
                    sqrt(kappa_TN / (kappa_TN + 1)) * h_los + ...
                    sqrt(1 / (kappa_TN + 1)) * h_nlos);
            end
        end
    end

    %% 2. Satellite Channel (Eq 4, 5, 6, 7)
    
    for m = 1:M
        for k = 1:K
            % Satellite Coordinates (Correct indexing for P.M > 1)
            sat_pos = P.stPos(m, :); 
            ue_pos = P.userPos(k, :);
            
            % 3D Distance
            d_horiz = norm(sat_pos(1:2) - ue_pos(1:2));
            d_vert = abs(sat_pos(3) - ue_pos(3));
            d_mk = sqrt(d_horiz^2 + d_vert^2);
            
            % Off-axis Angle (Theta) for Spot Beam
            % Angle between nadir (vertical) and user line-of-sight
            theta_mk = atan(d_horiz / d_vert);
            
            % Beam Gain - Eq (5) & (6) 
            if abs(theta_mk) < 1e-6
                G_m = P.G_ST_max;
            else
                % Eq (6): Lambda/u factor
                u = 2.07123 * sin(theta_mk) / sin(P.theta_3dB);
                
                if abs(u) < 1e-6
                     gain_factor = 1;
                else
                     % Eq (5): Bessel function pattern
                     term = (besselj(1,u)/(2*u)) + (36*besselj(3,u)/(u^3));
                     gain_factor = term^2;
                end
                G_m = P.G_ST_max * gain_factor;
            end
            
            % Path Loss (FSPL) - derived from Eq (4) 
            % FSPL = (c / 4*pi*f*d)^2
            FSPL_dB = 20*log10(d_mk) + 20*log10(P.fc) - 147.55;
            PL_st_lin = 10^(-FSPL_dB/10) * G_m * P.G_UE;
            
            % Channel Coefficient - Eq (7) [cite: 60]
            % Using random phase for snapshot simulation instead of explicit Doppler time-shift
            CH.h_st(m,k) = sqrt(PL_st_lin) * exp(1i * 2 * pi * rand); 
        end
    end
end