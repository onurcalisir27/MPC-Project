classdef RobotSimulator < handle
    % Simulator for robot manipulator control
    % Tests both MPC and MPC-ISM control schemes
    
    properties
        robot              % Robot model
        reference_traj     % Reference trajectory
        controllers        % Map of controllers to test
        dt                 % Simulation time step
        tf                 % Final simulation time
        results            % Simulation results
    end
    
    methods
        function this = RobotSimulator(robot, dt, tf)
            % Constructor
            this.robot = robot;
            this.dt = dt;
            this.tf = tf;
            this.controllers = containers.Map();
            this.results = containers.Map();
        end
        
        function setReferenceTrajectory(this, q_ref, dq_ref, t)
            % Set reference trajectory
            this.reference_traj.q = q_ref;   % Joint positions
            this.reference_traj.dq = dq_ref; % Joint velocities
            this.reference_traj.t = t;       % Time vector
        end
        
        function addController(this, name, controller)
            % Add a controller to test
            this.controllers(name) = controller;
        end
        
        function results = runSimulations(this)
            % Run simulations for all controllers
            
            % Loop through all controllers
            controller_names = this.controllers.keys;
            for i = 1:length(controller_names)
                name = controller_names{i};
                controller = this.controllers(name);
                
                % Run simulation with current controller
                disp(['Running simulation with controller: ' name]);
                this.results(name) = this.simulateControl(controller, name);
                
                % Reset controller if needed
                if isprop(controller, 'resetController') && ismethod(controller, 'resetController')
                    controller.resetController();
                end
            end
            
            % Return results
            results = this.results;
        end
        
        function result = simulateControl(this, controller, name)
            % Simulate robot control with a given controller
            
            % Time settings
            t_span = 0:this.dt:this.tf;
            n_steps = length(t_span);
            
            % Number of joints
            nu = length(this.robot.Link);
            
            % Initialize result structure
            result = struct();
            result.t = t_span;
            result.q = zeros(nu, n_steps);
            result.dq = zeros(nu, n_steps);
            result.q_ref = zeros(nu, n_steps);
            result.dq_ref = zeros(nu, n_steps);
            result.tau = zeros(nu, n_steps);
            result.u = zeros(nu, n_steps);
            result.error = zeros(nu, n_steps);
            result.controller_data = cell(1, n_steps);
            
            % Initial state
            q0 = zeros(nu, 1);
            dq0 = zeros(nu, 1);
            x = [q0; dq0];
            
            % Inverse dynamics controller (assuming it's part of the controller or separate)
            if isprop(controller, 'inverse_dynamics')
                inv_dyn = controller.inverse_dynamics;
            elseif isfield(controller, 'inverse_dynamics')
                inv_dyn = controller.inverse_dynamics;
            else
                inv_dyn = [];
            end
            
            % Simulation loop
            for k = 1:n_steps
                t = t_span(k);
                
                % Get reference at current time
                [q_ref, dq_ref] = this.interpolateReference(t);
                x_ref = [q_ref; dq_ref];
                
                % Store reference
                result.q_ref(:, k) = q_ref;
                result.dq_ref(:, k) = dq_ref;
                
                % Store current state
                result.q(:, k) = x(1:nu);
                result.dq(:, k) = x(nu+1:end);
                
                % Compute tracking error
                result.error(:, k) = q_ref - x(1:nu);
                
                % Compute control input
                if strcmp(name, 'MPC-ISM')
                    % For MPC-ISM controller
                    [tau, v, info] = controller.computeControl(t, x, x_ref);
                    
                    % Make sure v is a column vector
                    if size(v, 2) > 1
                        % If v is a matrix, take the first column
                        v = v(:,1);
                    end
                    
                    % Ensure tau is a column vector
                    if size(tau, 2) > 1
                        tau = tau(:,1);
                    end
                    
                    result.u(:, k) = v;
                    result.tau(:, k) = tau;
                    result.controller_data{k} = info;
                else
                    % For standard MPC or other controllers
                    if ismethod(controller, 'computeControl')
                        [u, info] = controller.computeControl(x, x_ref);
                        
                        % Make sure u is a column vector
                        if size(u, 2) > 1
                            u = u(:,1);
                        end
                        
                        result.u(:, k) = u;
                        
                        % If we have inverse dynamics, compute torques
                        if ~isempty(inv_dyn)
                            [tau, ~] = inv_dyn.computeControl(x(1:nu), x(nu+1:end), u);
                            
                            % Ensure tau is a column vector
                            if size(tau, 2) > 1
                                tau = tau(:,1);
                            end
                            
                            result.tau(:, k) = tau;
                        else
                            % Without inverse dynamics, assume u is torque
                            result.tau(:, k) = u;
                        end
                        
                        result.controller_data{k} = info;
                    else
                        error('Controller must have computeControl method');
                    end
                end
                
                % Simulate robot dynamics for one time step
                if ~isempty(inv_dyn)
                    % If we have inverse dynamics, use it to simulate
                    [t_out, x_out] = ode45(@(t, x) inv_dyn.simulateDynamics(t, x, result.tau(:, k)), ...
                        [t, t+this.dt], x);
                else
                    % Simplified double integrator model
                    [t_out, x_out] = ode45(@(t, x) this.simpleDynamics(t, x, result.u(:, k)), ...
                        [t, t+this.dt], x);
                end
                
                % Update state
                x = x_out(end, :)';
            end
            
            % Compute performance metrics
            result.metrics = this.computePerformanceMetrics(result);
        end
        
        function [q_ref, dq_ref] = interpolateReference(this, t)
            % Interpolate reference trajectory at time t
            
            % Find time indices
            t_ref = this.reference_traj.t;
            
            % Clamp time to valid range
            t = min(max(t, t_ref(1)), t_ref(end));
            
            % Interpolate
            q_ref = interp1(t_ref, this.reference_traj.q', t)';
            dq_ref = interp1(t_ref, this.reference_traj.dq', t)';
        end
        
        function dx = simpleDynamics(this, t, x, u)
            % Simple double integrator dynamics for testing
            % dx = [q_dot; q_ddot], q_ddot = u
            
            nu = length(u);
            q_dot = x(1:nu);
            
            % Double integrator dynamics
            dx = [x(nu+1:end); u];
        end
        
        function metrics = computePerformanceMetrics(this, result)
            % Compute performance metrics for analysis
            
            % Extract data
            q = result.q;
            q_ref = result.q_ref;
            u = result.u;
            error = result.error;
            
            % Number of joints
            nu = size(q, 1);
            
            % RMS tracking error for each joint
            rms_error = zeros(nu, 1);
            for i = 1:nu
                rms_error(i) = sqrt(mean(error(i, :).^2));
            end
            
            % Control effort (RMS of control inputs)
            control_effort = zeros(nu, 1);
            for i = 1:nu
                control_effort(i) = sqrt(mean(u(i, :).^2));
            end
            
            % Store metrics
            metrics.rms_error = rms_error;
            metrics.control_effort = control_effort;
            metrics.steady_state_error = mean(abs(error(:, end-10:end)), 2);
        end
        
        function visualizeResults(this)
            % Visualize and compare simulation results
            
            % Get controller names
            controller_names = this.results.keys;
            n_controllers = length(controller_names);
            
            % Colors for different controllers
            colors = lines(n_controllers);
            
            % 1. Joint positions and tracking errors
            figure('Name', 'Joint Positions and Tracking Errors');
            
            % Get number of joints from first result
            first_result = this.results(controller_names{1});
            nu = size(first_result.q, 1);
            
            % Plot for each joint
            for i = 1:nu
                % Position plot
                subplot(nu, 2, 2*i-1);
                hold on;
                
                % Reference
                plot(first_result.t, first_result.q_ref(i, :), 'k--', 'LineWidth', 1.5);
                
                % Actual positions for each controller
                for j = 1:n_controllers
                    result = this.results(controller_names{j});
                    plot(result.t, result.q(i, :), 'Color', colors(j, :), 'LineWidth', 1.5);
                end
                
                grid on;
                title(['Joint ' num2str(i) ' Position']);
                xlabel('Time (s)');
                ylabel('Position (rad)');
                if i == 1
                    legend(['Reference', controller_names], 'Location', 'best');
                end
                
                % Error plot
                subplot(nu, 2, 2*i);
                hold on;
                
                % Tracking errors for each controller
                for j = 1:n_controllers
                    result = this.results(controller_names{j});
                    plot(result.t, result.error(i, :), 'Color', colors(j, :), 'LineWidth', 1.5);
                end
                
                grid on;
                title(['Joint ' num2str(i) ' Tracking Error']);
                xlabel('Time (s)');
                ylabel('Error (rad)');
                if i == 1
                    legend(controller_names, 'Location', 'best');
                end
            end
            
            % 2. Control inputs
            figure('Name', 'Control Inputs');
            
            % Plot for each joint
            for i = 1:nu
                subplot(nu, 1, i);
                hold on;
                
                % Control inputs for each controller
                for j = 1:n_controllers
                    result = this.results(controller_names{j});
                    plot(result.t, result.u(i, :), 'Color', colors(j, :), 'LineWidth', 1.5);
                end
                
                grid on;
                title(['Joint ' num2str(i) ' Control Input']);
                xlabel('Time (s)');
                ylabel('Control (rad/s^2)');
                if i == 1
                    legend(controller_names, 'Location', 'best');
                end
            end
        end
    end
end