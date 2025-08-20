
%% Setup
clear all;
close all;
clc;

% Robot parameters
links = [1, 1, 0.5]; % 3 links with lengths L1=1, L2=1, L3=0.5
robot = Robot(links, [0, 0, 0]);

% System dimensions
nu = length(links); % Number of joints
nx = 2*nu;          % Full state dimension (positions and velocities)

% Constraints (from Table I in the paper)
state_constraints = [pi/2, pi/2, pi/2; 2, 3.5, 6.3]; % [position limits; velocity limits]
control_constraints = [145, 250, 350]; % acceleration limits

% Time settings
dt = 0.005;     % Simulation sampling time (from the paper: 0.001s)
T = 0.02;       % MPC sampling time (from the paper: 0.02s)
tf = 10;        % Total simulation time
t = 0:dt:tf;    % Time vector
n_steps = length(t);

%% Reference Trajectory Generation
% Generate desired trajectory: move from initial to target position
q_initial = [0, 0, 0]';
q_target = [pi/4, pi/3, pi/2]'; % Target from the paper

% Generate smooth trajectory using cubic polynomial interpolation
[q_trajectory, q_dot_trajectory, q_ddot_trajectory, ~] = ...
    robot.CubicPolynomialInterpolation(q_initial, q_target, tf, n_steps);

% Add terminal phase to hold position
terminal_phase_time = 5; % seconds
terminal_steps = terminal_phase_time / dt;
t_extended = [t, t(end) + dt:dt:t(end) + terminal_phase_time];

% Hold final position with zero velocity
q_trajectory_extended = [q_trajectory, repmat(q_target, 1, terminal_steps)];
q_dot_trajectory_extended = [q_dot_trajectory, zeros(nu, terminal_steps)];

% Update the time and trajectory variables
t = t_extended;
n_steps = length(t);
q_trajectory = q_trajectory_extended;
q_dot_trajectory = q_dot_trajectory_extended;

%% setup

% 1. Inverse Dynamics Controller
uncertainty_mag = [20, 30, 80]; % From the paper
inv_dyn = InverseDynamicsController(robot, uncertainty_mag);

% 2. MPC Controller (parameter values from the paper)
% Cost matrices (Q = diag(100, 100) and R = 0.1 for each joint as per paper)
Q = zeros(nx, nx);
R = zeros(nu, nu);
for i = 1:nu
    Q(2*i-1:2*i, 2*i-1:2*i) = diag([100, 100]); % Position and velocity weights
    R(i, i) = 0.1; % Control weight
end

% MPC parameters
N = 10;                           % Prediction horizon from the paper
umax = control_constraints';      % Control constraints
xmax = reshape(state_constraints, [], 1); % State constraints

% Create MPC controller
mpc_controller = MPController(Q, R, N, T, umax, xmax);

% 3. ISM Controller Parameters (values from the paper)
c = [10, 10, 10];              % Sliding variable parameters
Umax = [20, 35, 85];          % Maximum control values slightly less than uncertainty_mag
mu = [0.05, 0.05, 0.05];      % Filter time constants

% Create ISM controller
ism_controller = ISMController(c, Umax, mu);

% 4. Combined MPC-ISM Controller
mpc_ism_controller = InverseDynamicsISMMPCController(inv_dyn, mpc_controller, ism_controller);

%% Simulation Setup
% Create simulator
simulator = RobotSimulator(robot, dt, tf);

% Set reference trajectory
simulator.setReferenceTrajectory(q_trajectory, q_dot_trajectory, t);

% Add controllers to test
simulator.addController('MPC', mpc_controller);
simulator.addController('MPC-ISM', mpc_ism_controller);

%% run sim
disp('Starting simulations...');
results = simulator.runSimulations();
disp('Simulations completed!');

%% Visualize Results (only mpc ism)
simulator.visualizeResults();
mpc_ism_controller.visualizePerformance();

%% Print Performance Metrics
disp('Performance Metrics');
controller_names = results.keys;
for i = 1:length(controller_names)
    name = controller_names{i};
    disp(['Controller: ' name]);
    
    % RMS Tracking Error (eRMS from Table II)
    disp('RMS Tracking Error (rad):');
    for j = 1:nu
        fprintf('  Joint %d: %.4f\n', j, results(name).metrics.rms_error(j));
    end
    
    % Control Effort (Ec from Table II)
    disp('Control Effort (rad/s^2):');
    for j = 1:nu
        fprintf('  Joint %d: %.4f\n', j, results(name).metrics.control_effort(j));
    end
    
    % Additional metric: Steady-state Error
    disp('Steady-state Error (rad):');
    for j = 1:nu
        fprintf('  Joint %d: %.6f\n', j, results(name).metrics.steady_state_error(j));
    end
    
end

%% 3D Animation
% Choose the best controller for animation
best_controller = 'MPC-ISM';
best_result = results(best_controller);

% Calculate target position in Cartesian coordinates
robot.setJointAngle(q_target');
target_pos = robot.calcPosC();

% Animate the robot following the trajectory
disp('Starting 3D animation...');
figure;
robot.animate(best_result.q, target_pos');
title(['Robot Animation: ' best_controller ' Controller']);

%% Additional Comparison Plots
figure('Name', 'Joint Position Tracking Comparison');
subplot(3, 1, 1);
hold on;
plot(results('MPC').t, results('MPC').q(1, :), 'b-', 'LineWidth', 1);
plot(results('MPC-ISM').t, results('MPC-ISM').q(1, :), 'r-', 'LineWidth', 1);
plot(results('MPC').t, results('MPC').q_ref(1, :), 'k--', 'LineWidth', 1.5);
grid on;
title('Joint 1 Position Tracking');
xlabel('Time (s)');
ylabel('Position (rad)');
legend('MPC', 'MPC-ISM', 'Reference');

subplot(3, 1, 2);
hold on;
plot(results('MPC').t, results('MPC').q(2, :), 'b-', 'LineWidth', 1);
plot(results('MPC-ISM').t, results('MPC-ISM').q(2, :), 'r-', 'LineWidth', 1);
plot(results('MPC').t, results('MPC').q_ref(2, :), 'k--', 'LineWidth', 1.5);
grid on;
title('Joint 2 Position Tracking');
xlabel('Time (s)');
ylabel('Position (rad)');

subplot(3, 1, 3);
hold on;
plot(results('MPC').t, results('MPC').q(3, :), 'b-', 'LineWidth', 1);
plot(results('MPC-ISM').t, results('MPC-ISM').q(3, :), 'r-', 'LineWidth', 1);
plot(results('MPC').t, results('MPC').q_ref(3, :), 'k--', 'LineWidth', 1.5);
grid on;
title('Joint 3 Position Tracking');
xlabel('Time (s)');
ylabel('Position (rad)');

%% Compare Tracking Errors (as in Fig. 4 and Fig. 6 in the paper)
figure('Name', 'Tracking Error Comparison');
subplot(3, 1, 1);
hold on;
plot(results('MPC').t, abs(results('MPC').error(1, :)), 'b-', 'LineWidth', 1);
plot(results('MPC-ISM').t, abs(results('MPC-ISM').error(1, :)), 'r-', 'LineWidth', 1);
grid on;
title('Joint 1 Absolute Tracking Error');
xlabel('Time (s)');
ylabel('Error (rad)');
legend('MPC', 'MPC-ISM');

subplot(3, 1, 2);
hold on;
plot(results('MPC').t, abs(results('MPC').error(2, :)), 'b-', 'LineWidth', 1);
plot(results('MPC-ISM').t, abs(results('MPC-ISM').error(2, :)), 'r-', 'LineWidth', 1);
grid on;
title('Joint 2 Absolute Tracking Error');
xlabel('Time (s)');
ylabel('Error (rad)');

subplot(3, 1, 3);
hold on;
plot(results('MPC').t, abs(results('MPC').error(3, :)), 'b-', 'LineWidth', 1);
plot(results('MPC-ISM').t, abs(results('MPC-ISM').error(3, :)), 'r-', 'LineWidth', 1);
grid on;
title('Joint 3 Absolute Tracking Error');
xlabel('Time (s)');
ylabel('Error (rad)');

disp('All visualizations completed!');