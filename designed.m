classdef designed
    % This class contains information on designed oscillators.
    %
    %
    % properties (principals only):
    %
    % - dt: Time interval
    %
    % - order: Maximum order of polynomials
    %
    % - zeta1/zeta2: Polynomial coefficients of 1st/2nd component of vector field
    %
    %
    % functions: 
    %
    % - func: Dynamics of designed oscillator
    %
    % - func_jacobi: Jacobian matrix of designed oscillator

    properties
        dt;
        order;
        zeta1;
        zeta2;
        initial;
        x_tick;
        y_tick;
        x_lim;
        y_lim;
        name;
        y_basis;
    end
    methods
        function obj = designed(dt,order,xi,cla)
            obj.dt = dt;
            obj.order = order;
            obj.zeta1 = xi(1:end/2);
            obj.zeta2 = xi(end/2+1:end);
            obj.initial = cla.initial;
            obj.x_tick = cla.x_tick;
            obj.y_tick = cla.y_tick;
            obj.x_lim = cla.x_lim;
            obj.y_lim = cla.y_lim;
            obj.name = cla.name;
            obj.y_basis = cla.y_basis;
        end
        function F = func(obj,x)
            % input:
            % 
            % - x: State (matrix: 2*length)
            %
            %
            % output:
            %
            % - f: Dynamics at the state (matrix: 2*length)
            
            U = funcs.polynomial(x,obj.order);
            
            F = [U'*obj.zeta1,U'*obj.zeta2]';
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

            dUdx = funcs.polynomial_diff([x1_lc;x2_lc],obj.order);

            J = [dUdx(:,1)'*obj.zeta1,dUdx(:,2)'*obj.zeta1;dUdx(:,1)'*obj.zeta2,dUdx(:,2)'*obj.zeta2];
        end
    end
end