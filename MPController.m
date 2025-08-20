classdef MPController < handle

    
    properties
        Q           % State error cost weight
        R           % Control cost weight
        P           % Terminal cost weight
        N           % Prediction horizon
        dt          % Sampling time
        umax        % Maximum control input
        xmax        % Maximum state values
        nx          % State dimension
        nu          % Control dimension
        A           % Discrete system matrix
        B           % Discrete input matrix
        Xf          % Terminal constraint set (optional)
        KLQ         % LQR gain for terminal cost
        execution_time % For computational performance tracking
    end
    
    methods
        function this = MPController(Q, R, N, dt, umax, xmax)
            % Constructor
            this.Q = Q;
            this.R = R;
            this.N = N;
            this.dt = dt;
            this.umax = umax;
            this.xmax = xmax;
            
            % Extract dimensions
            this.nx = size(Q, 1);
            this.nu = length(umax);
            
            % Create system matrices for double integrator
            A_c = [zeros(this.nu), eye(this.nu); 
                   zeros(this.nu), zeros(this.nu)];
            B_c = [zeros(this.nu); eye(this.nu)];
            
            % Discretize system
            [this.A, this.B] = this.discretizeSystem(A_c, B_c, dt);
            
            % Compute LQR gain and terminal cost matrix
            [this.KLQ, this.P] = this.calculateTerminalCost();
            
            % Initialize performance tracking
            this.execution_time = [];
        end
        
        function [u, info] = computeControl(this, x, xref)
            % Compute MPC control input and return diagnostic info
            tic; % Start timing
            
            % Check dimensions
            if length(x) ~= this.nx || length(xref) ~= this.nx
                error('State or reference dimensions do not match expected dimensions');
            end
            
            % Compute tracking error (x - xref)
            error_state = x - xref;
            
            % Build prediction matrices for the entire horizon
            [Phi, Gamma] = this.buildPredictionMatrices();
            
            % Build cost function matrices for QP
            [H, f] = this.buildCostFunction(Phi, Gamma, x, xref);
            
            % Build constraint matrices for QP
            [Aineq, bineq, Aeq, beq] = this.buildConstraints(Phi, Gamma, x);
            
            % Solve the QP problem
            options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex');
            try
                [U, cost, exitflag, output] = quadprog(H, f, Aineq, bineq, Aeq, beq, [], [], [], options);
                
                % Check if solution was found
                if exitflag < 0
                    warning('QP solver failed with exit flag: %d', exitflag);
                    % Fallback to LQR controller if optimization fails
                    u = -this.KLQ * error_state;
                    u = min(max(u, -this.umax), this.umax); % Clip control input
                else
                    % Extract the first control input - receding horizon control
                    u = U(1:this.nu);
                end
                
            catch e
                warning('QP solver error: %s', e.message);
                % Fallback to LQR controller
                u = -this.KLQ * error_state;
                u = min(max(u, -this.umax), this.umax); % Clip control input
            end
            
            % Record execution time
            this.execution_time(end+1) = toc;
            
            % Prepare diagnostic information
            info.execution_time = this.execution_time(end);
            info.mean_execution_time = mean(this.execution_time);
            if exist('exitflag', 'var')
                info.exitflag = exitflag;
                info.cost = cost;
            else
                info.exitflag = -999; % Error code
                info.cost = NaN;
            end
        end
        
        function [A, B] = discretizeSystem(this, A_c, B_c, dt)
            % Discretize continuous-time system using zero-order hold
            nx = size(A_c, 1);
            nu = size(B_c, 2);
            
            % Using matrix exponential method
            M = expm([A_c, B_c; zeros(nu, nx), zeros(nu, nu)] * dt);
            A = M(1:nx, 1:nx);
            B = M(1:nx, nx+1:end);
        end
        
        function [Phi, Gamma] = buildPredictionMatrices(this)
            % Build state and input prediction matrices for the entire horizon
            A = this.A;
            B = this.B;
            N = this.N;
            nx = this.nx;
            nu = this.nu;
            
            % Initialize matrices
            Phi = zeros(nx*N, nx);    % Maps initial state to future states
            Gamma = zeros(nx*N, nu*N); % Maps control inputs to future states
            
            % Fill in the matrices block by block
            for i = 1:N
                % State prediction from initial state
                Phi((i-1)*nx+1:i*nx, :) = A^i;
                
                % Control effect on future states
                for j = 1:i
                    row_block = (i-1)*nx+1:i*nx;
                    col_block = (j-1)*nu+1:j*nu;
                    Gamma(row_block, col_block) = A^(i-j) * B;
                end
            end
        end
        
        function [H, f] = buildCostFunction(this, Phi, Gamma, x, xref)
            % Build cost function matrices for QP: min 0.5*U'*H*U + f'*U
            N = this.N;
            nx = this.nx;
            nu = this.nu;
            
            % Compute error = x - xref
            error_state = x - xref;
            
            % Extend reference over horizon
            x_ref_N = repmat(xref, N, 1);
            
            % Build block diagonal Q matrix
            Q_blk = kron(eye(N-1), this.Q);
            Q_blk = blkdiag(Q_blk, this.P);  % Terminal cost
            
            % Build block diagonal R matrix
            R_blk = kron(eye(N), this.R);
            
            % Compute quadratic term H
            H = Gamma' * Q_blk * Gamma + R_blk;
            
            % Make sure H is symmetric (handle numerical issues)
            H = (H + H')/2;
            
            % Compute linear term f
            f = Gamma' * Q_blk * (Phi * x - x_ref_N);
        end
        
        function [Aineq, bineq, Aeq, beq] = buildConstraints(this, Phi, Gamma, x)
            % Build constraint matrices for QP
            N = this.N;
            nx = this.nx;
            nu = this.nu;
            
            % Input constraints: -umax <= u <= umax
            I_Nu = eye(nu*N);
            Aineq_u = [I_Nu; -I_Nu];
            bineq_u = [repmat(this.umax, N, 1); repmat(this.umax, N, 1)];
            
            % State constraints: -xmax <= x <= xmax
            % We apply these to the predicted states X = Phi*x + Gamma*U
            predicted_x0_effect = Phi * x;  % Effect of initial state
            
            Aineq_x = [Gamma; -Gamma];
            bineq_x = [repmat(this.xmax, N, 1) - predicted_x0_effect; 
                       repmat(this.xmax, N, 1) + predicted_x0_effect];
                   
            % Combine all inequality constraints
            Aineq = [Aineq_u; Aineq_x];
            bineq = [bineq_u; bineq_x];
            
            % No equality constraints in this case
            Aeq = [];
            beq = [];
        end
        
        function [KLQ, P] = calculateTerminalCost(this)
            % Calculate LQR gain and terminal cost matrix P
            % Solves discrete-time Riccati equation for infinite horizon LQR
            
            % LQR design for terminal cost
            [KLQ, P] = dlqr(this.A, this.B, this.Q, this.R);
            
            % Note: dlqr returns K for u = -K*x, so negate it
            KLQ = -KLQ;
            
            % Verify that P satisfies the Lyapunov equation:
            % (A-B*K)'*P*(A-B*K) - P = -(Q + K'*R*K)
            % This is necessary for the terminal cost to ensure stability
            closed_loop = this.A + this.B*KLQ;
            residual = closed_loop'*P*closed_loop - P + this.Q + KLQ'*this.R*KLQ;
            
            % Ensure P is positive definite (for numerical stability)
            P = (P + P')/2;
            
            % Check numerical issues in Riccati solution
            if norm(residual, 'fro') > 1e-10
                warning('Terminal cost calculation may have numerical issues');
            end
        end
    end
end