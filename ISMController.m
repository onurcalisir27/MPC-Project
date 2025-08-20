classdef ISMController < handle

    
    properties
        c               % Sliding variable parameters (row vector [c_1, c_2, c_3])
        Umax            % Maximum control values (row vector [Umax_1, Umax_2, Umax_3])
        mu              % Filter time constants (row vector [mu_1, mu_2, mu_3])
        nu              % Number of joints (control inputs)
        t0              % Initial time
        x0              % Initial state
        prev_time       % Previous time step
        prev_u_ism      % Previous ISM control (filtered)
        prev_u_ism_raw  % Previous raw discontinuous control
        integral_term   % Integral term for sliding variable
        sliding_var     % Current sliding variable values
    end
    
    methods
        function this = ISMController(c, Umax, mu)
            % Constructor
            % c: Sliding variable parameter (vector for each joint)
            % Umax: Maximum control values (vector for each joint)
            % mu: Filter time constants (vector for each joint)
            
            this.c = c(:)';         % Ensure row vector
            this.Umax = Umax(:)';   % Ensure row vector
            this.mu = mu(:)';       % Ensure row vector
            this.nu = length(c);    % Number of joints
            
            % Initialize storage
            this.prev_u_ism = zeros(this.nu, 1);
            this.prev_u_ism_raw = zeros(this.nu, 1);
            this.integral_term = zeros(this.nu, 1);
            this.sliding_var = zeros(this.nu, 1);
            this.prev_time = 0;
            this.x0 = zeros(2*this.nu, 1); % Initialize x0 to prevent errors
        end
        
        function initialize(this, t0, x0)
            % Initialize controller with initial time and state
            % as per equation (15) in the paper
            this.t0 = t0;
            this.x0 = x0;
            this.prev_time = t0;
            this.integral_term = zeros(this.nu, 1);
            this.prev_u_ism = zeros(this.nu, 1);
            this.prev_u_ism_raw = zeros(this.nu, 1);
        end
        
        function [u_ism, sigma] = computeControl(this, x, u_mpc, t)
            % Compute ISM control input based on equations (13)-(15) in the paper
            % x: Current state [q; q_dot]
            % u_mpc: Control input from MPC
            % t: Current time
            
            % Extract positions and velocities
            q = x(1:this.nu);
            q_dot = x(this.nu+1:end);
            
            % Compute time step
            dt = t - this.prev_time;
            if dt <= 0
                dt = 1e-6; % Avoid division by zero
            end
            
            % Initialize output
            u_ism = zeros(this.nu, 1);
            sigma = zeros(this.nu, 1);
            
            % Ensure u_mpc is a column vector with the right dimensions
            u_mpc = reshape(u_mpc, [this.nu, 1]);
            
            % Compute ISM control for each joint
            for i = 1:this.nu
                % Position error relative to initial state
                position_error = q(i) - this.x0(i);
                
                % Update integral term based on equation (15) in the paper
                if t > this.t0
                    % Approximate the integral term using trapezoidal rule
                    v_term = u_mpc(i) + this.prev_u_ism(i);
                    this.integral_term(i) = this.integral_term(i) + ...
                        (q_dot(i) + v_term - this.prev_u_ism_raw(i)) * dt;
                end
                
                % Compute sliding variable according to equation (15)
                sigma(i) = this.c(i) * position_error + q_dot(i) - this.integral_term(i);
                
                % Store sliding variable
                this.sliding_var(i) = sigma(i);
                
                % Compute raw discontinuous control using equation (13)
                u_ism_raw_i = -this.Umax(i) * sign(sigma(i));
                
                % Apply filtering to reduce chattering using equation (14)
                if t == this.t0
                    % At initial time, use raw control
                    u_ism(i) = u_ism_raw_i;
                else
                    % Filter the discontinuous control
                    u_ism(i) = this.prev_u_ism(i) + ...
                        (dt/this.mu(i)) * (u_ism_raw_i - this.prev_u_ism(i));
                end
                
                % Store raw discontinuous control for next iteration
                this.prev_u_ism_raw(i) = u_ism_raw_i;
            end
            
            % Update storage for next iteration
            this.prev_time = t;
            this.prev_u_ism = u_ism;
            
            % Apply final saturation for safety
            u_ism = min(max(u_ism, -this.Umax'), this.Umax');
        end
        
        function resetIntegralTerm(this)
            % Reset the integral term (useful for restarting control)
            this.integral_term = zeros(this.nu, 1);
        end
        
        function sigma = getSlidingVariable(this)
            % Get current sliding variable values
            sigma = this.sliding_var;
        end
    end
end