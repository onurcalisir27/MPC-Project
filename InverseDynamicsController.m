classdef InverseDynamicsController < handle

    properties
        robot           % Reference to the robot model
        uncertainty_mag % Magnitude of uncertainty for each joint
    end
    
    methods
        function this = InverseDynamicsController(robot, uncertainty_mag)
            % Constructor
            this.robot = robot;
            
            % Default uncertainty magnitude if not provided (from the paper)
            if nargin < 2
                this.uncertainty_mag = [20, 30, 80]; % From the paper
            else
                this.uncertainty_mag = uncertainty_mag;
            end
        end
        
        function [tau, eta] = computeControl(this, q, q_dot, v)
            % Compute inverse dynamics control tau = M(q)v + n(q,q_dot)
            % and return the uncertainty term eta
            
            % Get inertia matrix
            M = this.computeInertiaMatrix(q);
            
            % Get nonlinear terms
            n = this.computeNonlinearTerms(q, q_dot);
            
            % Compute control torque using inverse dynamics
            tau = M * v + n;
            
            % Generate uncertainty (for simulation purposes)
            eta = this.generateUncertainty(q, q_dot);
        end
        
        function M = computeInertiaMatrix(this, q)
            % Simplified computation of inertia matrix
            % For a 3-link robot in 2D (planar robot)
            
            % Extract link lengths
            L = this.robot.Link;
            
            % For simplicity, use a diagonal inertia matrix with position-dependent terms
            M = diag([1 + 0.1*sin(q(1))^2, ...
                      1 + 0.2*L(2)*sin(q(2))^2, ...
                      1 + 0.1*L(3)*sin(q(3))^2]);
            
            % Add coupling between joints (off-diagonal terms)
            M(1,2) = 0.05 * sin(q(1)+q(2));
            M(2,1) = M(1,2);
            M(2,3) = 0.05 * sin(q(2)+q(3));
            M(3,2) = M(2,3);
        end
        
        function n = computeNonlinearTerms(this, q, q_dot)
            % Computation of nonlinear terms including:
            % Coriolis, centrifugal, friction, and gravity
            
            % Extract link lengths
            L = this.robot.Link;
            
            % Gravity terms
            g = 9.81; % Gravity constant
            gravity = [0; ...
                       g * L(2) * cos(q(2)); ...
                       g * L(3) * cos(q(2) + q(3))];
            
            % Coriolis and centrifugal terms
            coriolis = [0.1 * q_dot(1) * q_dot(2) * sin(q(2)); ...
                        -0.1 * q_dot(1)^2 * sin(q(2)) + 0.1 * q_dot(2) * q_dot(3) * sin(q(3)); ...
                        -0.1 * q_dot(2)^2 * sin(q(3))];
            
            % Friction (viscous and static)
            friction = 0.2 * q_dot + 0.1 * sign(q_dot);
            
            % Combine all nonlinear terms
            n = gravity + coriolis + friction;
        end
        
        function eta = generateUncertainty(this, q, q_dot)
            % Generate uncertainty term for simulation as described in the paper
            % The uncertainty bounds are [20, 30, 80] rad/s^2 for the three joints
            
            % Number of joints
            n_joints = length(q);
            
            % Generate uncertainty within bounds
            eta = zeros(n_joints, 1);
            for i = 1:n_joints
                % Maximum uncertainty amplitude from the paper
                max_eta = this.uncertainty_mag(i);
                
                % Create a combination of deterministic and random components
                % to stay within the bounds specified in the paper
                
                % Deterministic component based on position
                pos_comp = 0.5 * max_eta * sin(q(i));
                
                % Deterministic component based on velocity
                vel_comp = 0.2 * max_eta * sign(q_dot(i)) * min(abs(q_dot(i)), 1);
                
                % Random component (bounded)
                rand_comp = 0.3 * max_eta * (2*rand() - 1);
                
                % Combined uncertainty
                eta(i) = pos_comp + vel_comp + rand_comp;
                
                % Ensure bounds are respected
                eta(i) = min(max(eta(i), -max_eta), max_eta);
            end
        end
        
        function [x_dot] = simulateDynamics(this, t, x, u)
            % Simulate the robot dynamics for testing
            % x = [q; q_dot], u is joint torque
            
            % Extract position and velocity
            n_joints = length(u);
            q = x(1:n_joints);
            q_dot = x(n_joints+1:end);
            
            % Compute inverse dynamics for the current state
            M = this.computeInertiaMatrix(q);
            n = this.computeNonlinearTerms(q, q_dot);
            
            % Add uncertainty
            eta = this.generateUncertainty(q, q_dot);
            
            % Compute acceleration: q_ddot = M⁻¹(τ - n + η)
            q_ddot = M \ (u - n + eta);
            
            % Return the state derivative
            x_dot = [q_dot; q_ddot];
        end
        
        function [x_dot] = simulateLinearizedDynamics(this, t, x, v)
            % Simulate the linearized system with uncertainty
            % x = [q; q_dot], v is desired acceleration
            
            % Extract position and velocity
            n_joints = length(v);
            q = x(1:n_joints);
            q_dot = x(n_joints+1:end);
            
            % Generate uncertainty
            eta = this.generateUncertainty(q, q_dot);
            
            % The linearized system is a double integrator with uncertainty
            % q̈ = v - η
            q_ddot = v - eta;
            
            % Return the state derivative
            x_dot = [q_dot; q_ddot];
        end
    end
end