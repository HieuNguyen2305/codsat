clear; close all; clc;

% =========================================================================
%                         DISTRIBUTED ALGORITHM SETUP
% =========================================================================
fprintf('Initializing System Parameters for Distributed Simulation...\n');

% 1. Load Parameters & Channels
% params.m đã tự tính toán các hệ số Normalization (P.G0, P.sigma_norm...)
P = params();           
CH = gen_channels(P);   

fprintf('System Config: %d Users, %d BSs, %d Satellites\n', P.K, P.N, P.M);

% --- QUAN TRỌNG: KHÔNG ĐƯỢC SCALING THỦ CÔNG Ở ĐÂY ---
% (Đã bỏ đoạn nhân 1e10 để tránh xung đột với solve_bs_local/solve_st_local)

% =========================================================================
%                       STEP 1: INITIALIZATION
% =========================================================================
fprintf('\n--- Step 1: Finding Feasible Initial Point ---\n');
tic;
% Tìm điểm khởi tạo thỏa mãn QoS
p_init = find_init_point(P, CH); 
init_time = toc;
fprintf('Feasible Point Found in %.2f sec.\n', init_time);

% =========================================================================
%                       STEP 2: DISTRIBUTED ADMM-SCA
% =========================================================================
fprintf('\n--- Step 2: Start Distributed Dinkelbach-ADMM (Algorithm 1) ---\n');
tic;
% Chạy thuật toán phân tán (sẽ gọi solve_bs_local và solve_st_local)
result_admm = admm_dinkelbach(P, CH, p_init);
elapsed_time = toc;

% =========================================================================
%                             RESULT & PLOT
% =========================================================================
fprintf('\n================ FINAL RESULTS (DISTRIBUTED) ================\n');
fprintf('Execution Time:   %.2f sec\n', elapsed_time);
fprintf('Energy Efficiency: %.4f (Mbits/Joule)\n', result_admm.EE);
fprintf('Total Rate:        %.2f Mbps\n', result_admm.R_final);
fprintf('Total Power:       %.2f W\n', result_admm.Q_final);
fprintf('Final Prices:      lambda_BS = %.2e | lambda_ST = %.2e\n', ...
        result_admm.lambda_BS, result_admm.lambda_ST);

% --- Plot Convergence ---
if ~isempty(result_admm.hist_EE)
    figure('Name', 'Distributed Convergence', 'Color', 'w');
    plot(1:length(result_admm.hist_EE), result_admm.hist_EE, '-bo', ...
         'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    xlabel('Outer Iterations (Dinkelbach)');
    ylabel('Energy Efficiency (Mbits/Joule)');
    title('EE Maximization via Distributed ADMM');
    grid on;
end

% --- Visualization of Association (User Connectivity) ---
fprintf('\n--- User Association Map (Distributed Result) ---\n');
% Check BS connectivity
is_connected_BS = zeros(P.K, 1);
for k = 1:P.K
    % Sum power from all BSs on all Subcarriers
    p_rx_bs = sum(sum(result_admm.p_bs(:,k,:)));
    if p_rx_bs > 1e-5, is_connected_BS(k) = 1; end
end

% Check ST connectivity
is_connected_ST = zeros(P.K, 1);
for k = 1:P.K
    % Sum power from all Satellites
    p_rx_st = sum(result_admm.p_st(:,k));
    if p_rx_st > 1e-5, is_connected_ST(k) = 1; end
end

% Print Type
for k = 1:P.K
    if is_connected_BS(k) && is_connected_ST(k)
        type = 'Hybrid (BS + Satellite)';
    elseif is_connected_BS(k)
        type = 'Terrestrial Only';
    elseif is_connected_ST(k)
        type = 'Satellite Only';
    else
        type = 'OUTAGE (No Service)';
    end
    fprintf('User %d: %s\n', k, type);
end

fprintf('\nDistributed Simulation Completed Successfully.\n');