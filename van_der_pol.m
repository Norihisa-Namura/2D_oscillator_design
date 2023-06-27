classdef van_der_pol
    % This class contains information on van der Pol (vdP) oscillators.
    %
    %
    % properties (principals only):
    %
    % - nu: Parameter "nu"
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
    % - func: Dynamics of vdP oscillator
    %
    % - func_jacobi: Jacobian matrix of vdP oscillator

    properties
        nu = 3; % 3
        dt = 5e-3; % time interval
        order = 3; % polynomial order
        lambda_tol = -0.5; % tolerance value
        gamma = 1e-2; % regularization parameter
        initial = [2;0];
        x_tick = [-2,0,2];
        y_tick = [-6,-3,0,3,6];
        x_lim = [-3,3];
        y_lim = [-6,6];
        name = "vdp";
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

            f1 = x(2,:);
            f2 = obj.nu * (1 - x(1,:).^2) .* x(2,:) - x(1,:);
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

            J = [0,1;-1 + obj.nu*(-2*x1_lc)*x2_lc,obj.nu*(1 - x1_lc^2)];
        end
    end
end