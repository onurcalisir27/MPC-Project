classdef InverseDynamicsISMMPCController < handle

    properties
        inverse_dynamics  % Inverse dynamics controller
        mpc               % MPC controller
        ism               % ISM controller
        nu                % Number of joints (control inputs)
        initialized       % Flag to track initialization
        
        % For performance analysis
        u_mpc_history     % History of MPC control inputs
        u_ism_history     % History of ISM control inputs
        sigma_history     % History of sliding variables
        eta_history       % History of uncertainties
        time_history      % History of time steps
    end
    
    methods
        function this = InverseDynamicsISMMPCController(inverse_dynamics, mpc, ism)
            % Constructor
            % inverse_dynamics: Inverse dynamics controller for feedback linearization
            % mpc: Model Predictive Controller for optimal control
            % ism: Integral Sliding Mode controller for uncertainty compensation
            
            this.inverse_dynamics = inverse_dynamics;
            this.mpc = mpc;
            this.ism = ism;
            this.nu = numel(inverse_dynamics.robot.Link);
            this.initialized = false;
            
            % Initialize history arrays
            this.u_mpc_history = zeros(this.nu, 0);
            this.u_ism_history = zeros(this.nu, 0);
            this.sigma_history = zeros(this.nu, 0);
            this.eta_history = zeros(this.nu, 0);
            this.time_history = [];
        end
        
        function [tau, v, info] = computeControl(this, t, x, xref)
            % Compute control input using the hierarchical structure
            % t: Current time
            % x: Current state [q; q_dot]
            % xref: Reference state [qref; qref_dot]
            
            % Extract positions and velocities
            q = x(1:this.nu);
            q_dot = x(this.nu+1:end);
            
            % Initialize controller if this is the first call
            if ~this.initialized || t == 0
                this.initialized = true;
                % Initialize ISM controller with initial state
                this.ism.initialize(t, x);
            end
            
            % Step 1: Compute the MPC control input (outer loop)
            [u_mpc, mpc_info] = this.mpc.computeControl(x, xref);
            
            % Ensure it's a column vector
            u_mpc = reshape(u_mpc, [this.nu, 1]);
            
            % Step 2: Compute the ISM control input (middle loop)
            % This rejects the matched uncertainties not compensated by inverse dynamics
            try
                [u_ism, sigma] = this.ism.computeControl(x, u_mpc, t);
                
                % Ensure proper dimensions
                u_ism = reshape(u_ism, [this.nu, 1]);
                sigma = reshape(sigma, [this.nu, 1]);
            catch e
                fprintf('ISM computation at t=%f failed: %s\n', t, e.message);
                u_ism = zeros(this.nu, 1);
                sigma = zeros(this.nu, 1);
            end
            
            % Step 3: Combine MPC and ISM control components
            v = u_mpc + u_ism;
            
            % Step 4: Apply inverse dynamics to compute torques (inner loop)
            [tau, eta] = this.inverse_dynamics.computeControl(q, q_dot, v);
            
            % Store data for analysis
            if isempty(this.time_history)
                % First assignment
                this.u_mpc_history = u_mpc;
                this.u_ism_history = u_ism;
                this.sigma_history = sigma;
                this.eta_history = eta;
                this.time_history = t;
            else
                % Append to existing arrays
                this.u_mpc_history(:, end+1) = u_mpc;
                this.u_ism_history(:, end+1) = u_ism;
                this.sigma_history(:, end+1) = sigma;
                this.eta_history(:, end+1) = eta;
                this.time_history(end+1) = t;
            end
            
            % Return additional information
            info.u_mpc = u_mpc;
            info.u_ism = u_ism;
            info.sigma = sigma;
            info.eta = eta;
            info.mpc_info = mpc_info;
            info.v = v;
        end
        
        function visualizePerformance(this)
            % Visualize controller performance
            
            % Skip if no data is available
            if isempty(this.time_history)
                warning('No data available for visualization');
                return;
            end
            
            % 1. Plot control inputs (MPC, ISM, and combined)
            figure('Name', 'MPC-ISM Control Components');
            for i = 1:this.nu
                subplot(this.nu, 1, i);
                plot(this.time_history, this.u_mpc_history(i, :), 'b-', 'LineWidth', 2);
                hold on;
                plot(this.time_history, this.u_ism_history(i, :), 'r-', 'LineWidth', 2);
                plot(this.time_history, this.u_mpc_history(i, :) + this.u_ism_history(i, :), 'k--', 'LineWidth', 1.5);
                grid on;
                title(['Control Input for Joint ' num2str(i)]);
                xlabel('Time (s)');
                ylabel('Acceleration (rad/s^2)');
                legend('MPC', 'ISM', 'Combined');
            end
            
            % 2. Plot uncertainties vs ISM compensation
            figure('Name', 'Uncertainty Compensation');
            for i = 1:this.nu
                subplot(this.nu, 1, i);
                plot(this.time_history, this.eta_history(i, :), 'b-', 'LineWidth', 2);
                hold on;
                plot(this.time_history, this.u_ism_history(i, :), 'r-', 'LineWidth', 2);
                grid on;
                title(['Uncertainty vs ISM for Joint ' num2str(i)]);
                xlabel('Time (s)');
                ylabel('Value (rad/s^2)');
                legend('Uncertainty', 'ISM Compensation');
            end
            
            % 3. Plot sliding variables
            figure('Name', 'Sliding Variables');
            for i = 1:this.nu
                subplot(this.nu, 1, i);
                plot(this.time_history, this.sigma_history(i, :), 'LineWidth', 2);
                grid on;
                title(['Sliding Variable \sigma for Joint ' num2str(i)]);
                xlabel('Time (s)');
                ylabel('Value');
                
                % Add horizontal line at zero
                hold on;
                plot([this.time_history(1), this.time_history(end)], [0, 0], 'k--');
                hold off;
            end
        end
        
        function data = getPerformanceData(this)
            % Return performance data for analysis
            data.time = this.time_history;
            data.u_mpc = this.u_mpc_history;
            data.u_ism = this.u_ism_history;
            data.sigma = this.sigma_history;
            data.eta = this.eta_history;
        end
        
        function resetController(this)
            % Reset the controller
            this.initialized = false;
            
            % Clear history
            this.u_mpc_history = zeros(this.nu, 0);
            this.u_ism_history = zeros(this.nu, 0);
            this.sigma_history = zeros(this.nu, 0);
            this.eta_history = zeros(this.nu, 0);
            this.time_history = [];
            
            % Reset ISM controller
            this.ism.resetIntegralTerm();
        end
    end
end