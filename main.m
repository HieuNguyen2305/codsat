clear; close all; clc;

% =========================================================================
%                               SETUP & INIT
% =========================================================================
fprintf('Initializing System Parameters (Paper Section II)...\n');

% 1. Load Parameters
P = params();           

% 2. Generate Channels (Rician for TN, Shadowed-Rician for Satellite)
CH = gen_channels(P);   

fprintf('System Config: %d Users, %d BSs, %d Satellites, %d Subcarriers\n', ...
        P.K, P.N, P.M, P.S);

% --- STEP 1: FIND FEASIBLE POINT (Section III.G) ---
% Solves the relaxed problem to find a starting point that satisfies 
% User Association and QoS constraints.
fprintf('\n--- Step 1: Finding Feasible Initial Point ---\n');
tic;
p_init = find_init_point(P, CH); 
init_time = toc;
fprintf('Feasible Point Found in %.2f sec.\n', init_time);

% =========================================================================
%                             OPTIMIZATION
% =========================================================================
% --- STEP 2: DINKELBACH-SCA ALGORITHM (Algorithm 1) ---
% Solves the Energy Efficiency Maximization problem.
fprintf('\n--- Step 2: Start Centralized Dinkelbach-SCA ---\n');
tic;
result = dinkelbach_sca(P, CH, p_init);
elapsed_time = toc;

% =========================================================================
%                             RESULT & PLOT
% =========================================================================
fprintf('\n================ FINAL RESULTS ================\n');
fprintf('Execution Time: %.2f sec\n', elapsed_time);
fprintf('Max Energy Efficiency (EE): %.4f (Mbits/Joule)\n', result.EE);
fprintf('Total System Rate:        %.2f Mbps\n', result.R_final);
fprintf('Total Power Consumption:  %.2f W\n', result.Q_final);

% --- Plot Convergence ---
if ~isempty(result.hist_EE)
    figure('Name', 'EE Convergence');
    plot(1:length(result.hist_EE), result.hist_EE, '-rs', ...
         'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    xlabel('Dinkelbach Iterations');
    ylabel('Energy Efficiency (Mbits/Joule)');
    grid on;
    title('Convergence of Proposed Algorithm');
    xlim([1, length(result.hist_EE)]);
else
    fprintf('Warning: No convergence history to plot.\n');
end

% --- Visualization of Association (User Connectivity) ---
fprintf('\n--- User Association Map (After Recovery) ---\n');
% Based on Eq (8) and (12) recovered from p_opt
for k = 1:P.K
    % Check if connected to any BS
    % alpha is N x K x S. We sum over N and S to see if user k is served.
    is_connected_BS = sum(sum(result.alpha(:,k,:))) > 0;
    
    % Check if connected to any Satellite
    % beta is M x K.
    is_connected_ST = sum(result.beta(:,k)) > 0;
    
    if is_connected_BS && is_connected_ST
        type = 'Hybrid (BS + Satellite)';
    elseif is_connected_BS
        type = 'Terrestrial Only';
    elseif is_connected_ST
        type = 'Satellite Only';
    else
        type = 'OUTAGE (No Service)';
    end
    
    fprintf('User %d: %s\n', k, type);
end

fprintf('\nSimulation Completed Successfully.\n');