classdef star
    % This class contains information on star-shaped oscillators.
    %
    %
    % properties (principals only):
    %
    % - T: Period of oscillation
    %
    % - M: Number of time steps in one oscillation
    %
    % - dt: Time interval
    %
    % - p: Periodic trajectory
    %
    % - Z: Phase sensitivity function (PSF)
    %
    % - dpdt: Velocity on the periodic trajectory
    %
    % - dZdt: Time derivative of the PSF
    %
    % - Gamma: Phase coupling function (PCF) analytically obtained in the papar 
    %          (Depends on the input to the oscillator)
    %
    % - stable: Stable fixed points of the PCF
    %
    % - unstable: Unstable fixed points of the PCF
    %
    % - order: Maximum order of polynomials
    %
    % - lambda_tol: Tolerance value for linear stability
    %
    % - gamma: Regularization parameter for coefficients of polynomials

    properties
        T;
        M = 1000;
        t;
        dt;
        p;
        Z;
        dpdt;
        dZdt;
        Gamma;
        stable;
        unstable;
        order = 10;
        lambda_tol = -1;
        gamma = 1e-3;
        initial;
        x_tick = [-2,-1,0,1,2];
        y_tick = [-2,-1,0,1,2];
        x_lim = [-2,2];
        y_lim = [-2,2];
        name = "star";
        r1 = sqrt(2);
        r2;
        m = 4;
        evol_basin = 12;
        y_basis = 0;
    end
    methods
        function obj = star()
            obj.r2 = sqrt(obj.r1^2 - 1) / obj.m;
            obj.T = 2*pi;
            obj.dt = obj.T/obj.M;
            obj.t = linspace(0,obj.T-obj.dt,obj.M);
            obj.p = [obj.r1*cos(obj.t) + obj.r2*sin(obj.m*obj.t); obj.r1*sin(obj.t) + obj.r2*cos(obj.m*obj.t)];
            obj.dpdt = [-obj.r1*sin(obj.t) + obj.m*obj.r2*cos(obj.m*obj.t); obj.r1*cos(obj.t) - obj.m*obj.r2*sin(obj.m*obj.t)];
            obj.Z = [-obj.r1*sin(obj.t) - obj.m*obj.r2*cos(obj.m*obj.t); obj.r1*cos(obj.t) + obj.m*obj.r2*sin(obj.m*obj.t)];
            obj.dZdt = [-obj.r1*cos(obj.t) + obj.m^2*obj.r2*sin(obj.m*obj.t); -obj.r1*sin(obj.t) + obj.m^2*obj.r2*cos(obj.m*obj.t)];
            obj.Gamma = sqrt(2)*cos(obj.t.');
            obj.initial = obj.p(:,1);
            obj.stable = [pi/2];
            obj.unstable = [3*pi/2];
        end
    end
end