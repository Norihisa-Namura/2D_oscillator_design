classdef stuart_landau
    % This class contains information on Stuart-Landau (SL) oscillators.
    %
    %
    % properties (principals only):
    %
    % - alpha: Parameter "alpha"
    %
    % - beta: Parameter "beta"
    %
    % - dt: Time interval
    %
    % - order: Maximum order of polynomials
    %
    % - lambda_tol: Tolerance value for linear stability
    %
    % - gamma: Regularization parameter for coefficients of polynomials
    %
    %
    % functions: 
    %
    % - func: Dynamics of SL oscillator
    %
    % - func_jacobi: Jacobian matrix of SL oscillator

    properties
        alpha = 2; % 2
        beta = 1; % 1
        dt = 2*pi/2000; % time interval
        order = 3; % polynomial order
        lambda_tol = -1; % tolerance value
        gamma = 1e-0; % regularization parameter
        initial = [1;0];
        x_tick = [-1,0,1];
        y_tick = [-1,0,1];
        x_lim = [-1.5,1.5];
        y_lim = [-1.5,1.5];
        name = "sl";
        y_basis = 0;
    end

    methods
        function f = func(obj,x)
            % input:
            % 
            % - x: State (matrix: 2*length)
            %
            %
            % output:
            %
            % - f: Dynamics at the state (matrix: 2*length)

            f1 = x(1,:) - obj.alpha * x(2,:) - (x(1,:) - obj.beta * x(2,:)) .* (x(1,:).^2 + x(2,:).^2);
            f2 = obj.alpha * x(1,:) + x(2,:) - (obj.beta * x(1,:) + x(2,:)) .* (x(1,:).^2 + x(2,:).^2);
            f = [f1;f2];
        end
        
        function J = func_jacobi(obj,x1_lc,x2_lc)
            % input:
            % 
            % - x1_lc/x2_lc: 1st/2nd component of state (scalar)
            %
            %
            % output:
            %
            % - J: Jacobian matrix (matrix: 2*2)

            x12 = x1_lc.^2 + x2_lc.^2;
            J = [1 - x12 - 2*x1_lc*(x1_lc - obj.beta*x2_lc),-obj.alpha + obj.beta*x12 - 2*x2_lc*(x1_lc - obj.beta*x2_lc);...
               obj.alpha - obj.beta*x12 - 2*x1_lc*(obj.beta*x1_lc + x2_lc),1 - x12 - 2*x2_lc*(obj.beta*x1_lc + x2_lc)];
        end
    end
end