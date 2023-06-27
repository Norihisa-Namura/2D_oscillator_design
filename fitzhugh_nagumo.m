classdef fitzhugh_nagumo
    % This class contains information on FitzHugh-Nagumo (FHN) oscillators.
    %
    %
    % properties (principals only):
    %
    % - a: Parameter "a"
    %
    % - b: Parameter "b"
    %
    % - c: Parameter "c"
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
    % - func: Dynamics of FHN oscillator
    %
    % - func_jacobi: Jacobian matrix of FHN oscillator

    properties
        a = -0.1;
        b = 0.5;
        c = 0.01;
        dt = 0.1; % time interval
        order = 7; % polynomial order
        lambda_tol = -0.5; % tolerance value
        gamma = 1e-7; % regularization parameter
        initial = [0.95;0.05];
        x_tick = [-0.5,0,0.5,1];
        y_tick = [-0.1,0,0.1,0.2,0.3];
        x_lim = [-0.5,1];
        y_lim = [-0.1,0.3];
        name = "fhn";
        y_basis = 0;
    end
    
    methods
        function dxdt = func(obj,x)
            % input:
            % 
            % - x: State (matrix: 2*length)
            %
            %
            % output:
            %
            % - f: Dynamics at the state (matrix: 2*length)

            dx1dt = x(1,:) .* (x(1,:) - obj.a) .* (1 - x(1,:)) - x(2,:);
            dx2dt = (x(1,:) - obj.b * x(2,:)) * obj.c;
            dxdt = [dx1dt;dx2dt];
        end
        
        function J = func_jacobi(obj,x1_lc,x2_lc)
            % input:
            % 
            % - x1_lc: 1st component of state (scalar)
            %
            % - x2_lc (not necessary but needed to align function format): 2nd component of state (scalar)
            %
            %
            % output:
            %
            % - J: Jacobian matrix (matrix: 2*2)

            J = [-3*x1_lc.^2 + 2*(obj.a+1)*x1_lc - obj.a,-1;1*obj.c,-obj.b*obj.c];
        end
    end
end