function P = params()
% PARAMS: High Density Scenario (K=12 Users)
% Mục đích: Test khả năng chịu tải và quản lý nhiễu của thuật toán

    rng(99); % Đổi seed một chút để tạo vị trí mới ngẫu nhiên

    %% 1. Topology & Spectral
    P.K = 6;   % <--- TĂNG TỪ 6 LÊN 12 USER
    P.N = 3;    
    P.M = 2;    
    P.S = 8;    
    
    P.fc = 3.4e9;           
    P.c  = 3e8;             
    P.W_SC_Hz = 180e3;          
    P.W_ST_Hz = 20e6;           
    P.alpha_BW = P.W_SC_Hz / P.W_ST_Hz; 

    %% 2. Power Constraints (Watts)
    P.pBS_max = 40;  
    P.PST_max = 30;             
    
    % Noise Calculation
    N0_dBm_Hz = -174;
    NF_BS = 10; NF_ST = 10; 
    P.sigma_BS_sq = 10^((N0_dBm_Hz + 10*log10(P.W_SC_Hz) + NF_BS - 30)/10); 
    P.sigma_ST_sq = 10^((N0_dBm_Hz + 10*log10(P.W_ST_Hz) + NF_ST - 30)/10); 

    %% 3. Normalization
    Ref_Noise = min(P.sigma_BS_sq, P.sigma_ST_sq);
    P.G0 = 1 / Ref_Noise; 
    P.sigma_BS_norm = P.sigma_BS_sq * P.G0;
    P.sigma_ST_norm = P.sigma_ST_sq * P.G0;
    
    P.RateUnit = 1e6; 
    P.W_SC_Norm = P.W_SC_Hz / P.RateUnit;
    P.W_ST_Norm = P.W_ST_Hz / P.RateUnit;

    %% 4. Antenna & Pathloss
    P.G_BS = 10^(15/10);        
    P.G_UE = 1;                 
    P.G_ST_max = 10^(35/10);    
    P.theta_3dB = deg2rad(0.4); 

    %% 5. Position (Mixed Near & Far Scaled)
    P.bsPos = [0, 0, 30; 500, 0, 30; 250, 433, 30]; 
    P.stPos = [250, 250, 600e3; 250, 250 + 100e3, 600e3]; 
    
    % Phân bố lại User: 50% Gần (Ưu tiên BS), 50% Xa (Ưu tiên ST)
    K_near = floor(P.K / 2); % 6 User
    K_far = P.K - K_near;    % 6 User
    
    % Nhóm 1: Gần BS (0 - 500m)
    pos_near = 500 * rand(K_near, 2); 
    
    % Nhóm 2: Xa BS (1000 - 1500m) -> Cell Edge / Satellite Coverage
    pos_far  = 1000 + 500 * rand(K_far, 2); 
    
    P.userPos = [pos_near; pos_far];
    P.userPos(:,3) = 1.5; 

    %% 6. Energy Consumption Model
    P.P0_n = 40;   
    P.psi_n = 5;   
    P.P0_m = 50;   
    P.psi_m = 10;  

    %% 7. Algorithm Hyperparameters
    P.zeta = 10;            
    P.S_bar = 8;            
    P.R_threshold_Mbps = 0.5; % Giữ nguyên QoS 0.5 Mbps
    
    % Cấu hình vẽ đồ thị mượt (Smooth Plot)
    P.maxDinkel = 40;       
    P.tolDinkel = 1e-3;     
    P.maxSCA = 10;           
    P.tolSCA = 1e-3;
    P.prox_weight = 1e-3;   
end