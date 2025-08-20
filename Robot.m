classdef Robot < handle
    % Robot - A simplified 3-link robot manipulator for MPC control
    
    properties
        Link          % Link lengths; a vector [L1, L2, L3]
        Mass          % Link masses; a vector [m1, m2, m3]
    end
    
    properties (Access = private)
        JointAngle    % Angles of each joint; a vector [q1, q2, q3]
    end
    
    methods
        %% Constructor
        function this = Robot(link, initial_joint_angles)
            if nargin < 1
                error('Link lengths must be specified.');
            end
            
            this.Link = link; % Set link lengths
            this.Mass = [1, 1, 1]; % Unit mass for simplification
            
            if nargin < 2
                % Default joint angles to zeros if not provided
                this.JointAngle = zeros(1, length(link));
            else
                % Set initial joint angles
                this.JointAngle = initial_joint_angles;
            end
        end
        
        % Setter for Joint Angles
        function setJointAngle(this, value)
            this.JointAngle = value;
        end
        
        % Getter for Joint Angles
        function value = getJointAngle(this)
            value = this.JointAngle;
        end
        
        %% Forward Kinematics - Simplified for 3 links
        function FK = forward_kinematics(this, dh_table)
            FK = eye(4);
            for i = 1:length(dh_table(:,1))
                H = FK_helper(this, dh_table(i, :));
                FK = FK*H;
            end
            FK = simplify(FK);
        end
        
        % FK Helper Function
        function H = FK_helper(this, link_params)
            a = link_params(1);
            alpha = link_params(2);
            d = link_params(3);
            theta = link_params(4);
            H = [cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta);
                 sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
                          0,             sin(alpha),             cos(alpha),            d;
                          0,                      0,                      0,            1];
        end
        
        % Numeric FK Calculation
        function T = Numeric_FK(this, FK, q_des)
            % Define symbolic variables for joint angles and link lengths
            syms theta1 theta2 theta3 L1 L2 L3
            
            % Extract link lengths
            links = this.Link;
            
            % Substitute the joint angles (q_des) into the symbolic FK expression
            P = subs(FK, [theta1, theta2, theta3], q_des);
            
            % Substitute the link lengths into the FK expression
            P = subs(P, [L1, L2, L3], links);
            
            % Simplify and evaluate the result
            T = double(vpa(P)); % Use vpa to handle symbolic numbers before conversion
        end

        %% Position Calculation Methods - Modified for 3 links
        % Position of the end of Link 1
        function posA = calcPosA(this)
            L = this.Link;
            x = 0;
            y = 0;
            z = L(1);
            posA = [x; y; z];
        end
        
        % Position of the end of Link 2
        function posB = calcPosB(this)
            L = this.Link;
            q = this.JointAngle;
            
            x = L(2)*cos(q(1))*cos(q(2));
            y = L(2)*cos(q(2))*sin(q(1)); 
            z = L(1) + L(2)*sin(q(2));   
            posB = [x; y; z];
        end
        
        % Position of the end of Link 3 (end effector)
        function posC = calcPosC(this)
            L = this.Link;
            q = this.JointAngle;
            
            x = cos(q(1))*(L(3) * cos(q(2) + q(3)) + L(2) * cos(q(2)));
            y = sin(q(1))*(L(3) * cos(q(2) + q(3)) + L(2) * cos(q(2)));
            z = L(1) + L(3) * sin(q(2) + q(3)) + L(2) * sin(q(2));
            posC = [x; y; z];
        end
        
        %% Jacobian Calculation - Modified for 3 links
        function J = calcJacobian(this, q)
            % Define small perturbation for numeric differentiation
            delta = 1e-6;
            
            % Preallocate Jacobian matrix
            J = zeros(3, length(q));
            
            % Get current end-effector position
            posC = this.calcPosC(); % Current end-effector position
            
            % Loop through each joint to compute partial derivatives
            for i = 1:length(q)
                % Create perturbed joint angle vector
                q_perturbed = q;
                q_perturbed(i) = q_perturbed(i) + delta;
                
                % Update joint angles
                this.setJointAngle(q_perturbed);
                
                % Compute perturbed end-effector position
                pos_perturbed = this.calcPosC();
                
                % Numeric differentiation (finite difference)
                J(:, i) = (pos_perturbed - posC) / delta;
            end
            
            % Reset joint angles to the original values
            this.setJointAngle(q);
        end
        
        %% Inverse Kinematics - Modified for 3 links
        function Thetas = IK(this, dx, dy, dz, phi)
            % Extract link lengths
            L = this.Link;
            
            % Calculate theta1
            theta1 = atan2(dy, dx);
            
            % Intermediate terms
            d = dx - L(3) * cos(theta1) * cos(phi);
            e = dy - L(3) * sin(theta1) * cos(phi);
            f = dz - L(1) - L(3) * sin(phi);
            
            % Calculate theta3
            cos_theta3 = (d^2 + e^2 + f^2 - L(2)^2 - L(3)^2) / (2 * L(2) * L(3));
            cos_theta3 = min(max(cos_theta3, -1), 1);
            theta3 = acos(cos_theta3);
            
            % Intermediate terms for theta2
            a = L(3) * sin(theta3);
            b = L(2) + L(3) * cos(theta3);
            c = dz - L(1) - L(3) * sin(phi);
            r = sqrt(a^2 + b^2);
            
            % Calculate theta2
            sqrt_term = max(0, r^2 - c^2); % Ensure no complex numbers
            theta2 = atan2(c, sqrt(sqrt_term)) - atan2(a, b);
            
            % Return joint angles
            Thetas = [theta1, theta2, theta3];
        end
        
        %% Trajectory Planning
        function [q_trajectory, q_dot_trajectory, q_ddot_trajectory, coeffs] = CubicPolynomialInterpolation(this, q_initial, q_final, t, n)          
            T = linspace(0, t, n); 
            num_joints = length(q_initial); 
            q_trajectory = zeros(num_joints, n); 
            q_dot_trajectory = zeros(num_joints, n); 
            q_ddot_trajectory = zeros(num_joints, n); 

            % Coefficients matrix A for cubic polynomial
            A = [0,  0,  0, 1;        % Position at t = 0
                 t^3, t^2, t, 1;      % Position at t = t
                 0,  0,  1, 0;        % Velocity at t = 0
                 3*t^2, 2*t, 1, 0];   % Velocity at t = t
        
            for i = 1:num_joints
                % Boundary conditions: positions and velocities
                B = [q_initial(i); q_final(i); 0; 0]; % Zero velocities at start and end
        
                % Solve for polynomial coefficients
                coeffs = A \ B;
        
                % Generate the trajectory for joint i
                q_trajectory(i, :) = coeffs(1) * T.^3 + coeffs(2) * T.^2 + ...
                                     coeffs(3) * T + coeffs(4);
        
                % Generate the velocities for joint i (first derivative)
                q_dot_trajectory(i, 2:end-1) = 3 * coeffs(1) * T(2:end-1).^2 + ...
                                               2 * coeffs(2) * T(2:end-1) + ...
                                               coeffs(3);
        
                % Generate the accelerations for joint i (second derivative)
                q_ddot_trajectory(i, 2:end-1) = 6 * coeffs(1) * T(2:end-1) + ...
                                               4 * coeffs(2);
                                         
                % Ensure zero velocity at start and end
                q_dot_trajectory(i, 1) = 0;
                q_dot_trajectory(i, end) = 0;
                q_ddot_trajectory(i, 1) = 0;
                q_ddot_trajectory(i, end) = 0;
            end
        end
        
        %% Animation Function - Simplified for 3 links
        function animate(this, q_trajectory, targets)
            % Set up 3D figure
            figure;
            hold on;
            
            % Setup ground plane
            fill3([-5, 5, 5, -5], [-5, -5, 5, 5], [0, 0, 0, 0], 'c', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); 
            
            % Set up 3D visualization
            axis equal;
            grid on;
            xlim([-3, 3]); ylim([-3, 3]); zlim([-1, 3]);
            xlabel('X'); ylabel('Y'); zlabel('Z');
            view(45, 30);
            title('Robot Simulation')
            
            % Number of frames in the trajectory
            num_frames = size(q_trajectory, 2);
            
            % Animation loop
            for i = 1:num_frames
                % Update joint angles
                this.setJointAngle(q_trajectory(:, i)');
                
                % Compute positions of the links
                posA = this.calcPosA();
                posB = this.calcPosB();
                posC = this.calcPosC();
                
                % Clear and re-plot the robot
                cla;
                
                % Re-plot ground plane
                fill3([-5, 5, 5, -5], [-5, -5, 5, 5], [0, 0, 0, 0], 'c', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                
                % Plot links
                plot3([0, posA(1)], [0, posA(2)], [0, posA(3)], 'k', 'LineWidth', 3); % Link 1
                plot3([posA(1), posB(1)], [posA(2), posB(2)], [posA(3), posB(3)], 'k', 'LineWidth', 3); % Link 2
                plot3([posB(1), posC(1)], [posB(2), posC(2)], [posB(3), posC(3)], 'k', 'LineWidth', 3); % Link 3
                
                % Plot joints
                plot3(0, 0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); % Base
                plot3(posA(1), posA(2), posA(3), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); % Joint A
                plot3(posB(1), posB(2), posB(3), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); % Joint B
                plot3(posC(1), posC(2), posC(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % End effector
                
                % Plot target points
                for j = 1:size(targets, 1)
                    plot3(targets(j, 1), targets(j, 2), targets(j, 3), 'g*', 'MarkerSize', 10);
                    text(targets(j, 1), targets(j, 2), targets(j, 3) + 0.1, sprintf('T%d', j), 'FontSize', 10);
                end
                
                % Update display
                drawnow;
                pause(0.01);
            end
        end
    end
end