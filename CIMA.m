classdef CIMA
    % This class contains information on Chlorite-iodide-malonic acid (CIMA) oscillators.
    %
    %
    % properties (principals only):
    %
    % - a: Parameter "a"
    %
    % - b: Parameter "b"
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
    % - func: Dynamics of CIMA oscillator
    %
    % - func_jacobi: Jacobian matrix of CIMA oscillator

    properties
        a = 10;
        b = 2;
        dt = 0.001; % time interval
        order = 5; % polynomial order
        lambda_tol = -1; % tolerance value
        gamma = 1e-5; % regularization parameter
        initial = [2;4];
        x_tick = [0,2,4];
        y_tick = [2,5,8];
        x_lim = [0,4];
        y_lim = [2,8];
        name = "CIMA";
        y_basis = 4;
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

            f1 = obj.a - x(1,:) - 4*x(1,:).*x(2,:) ./ (1 + x(1,:).^2);
            f2 = obj.b * x(1,:) .* (1 - x(2,:) ./ (1 + x(1,:).^2));
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

            tmp = 1 + x1_lc^2;
            J = [-1 - 4*x2_lc*(1 - x1_lc^2) / tmp^2,-4*x1_lc / tmp;...
               obj.b * (1 + (-x2_lc + x1_lc^2*x2_lc) / tmp^2),-obj.b*x1_lc / tmp];
        end
    end
end