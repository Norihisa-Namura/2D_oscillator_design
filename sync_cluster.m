classdef sync_cluster
    % This class contains information on oscillators with high-harmonic PSF.
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
    % - Gamma_5: Phase coupling function (PCF) for 5 clusters
    %
    % - Gamma_3: Phase coupling function (PCF) for 3 clusters
    %
    % - stable_5: Stable fixed points of the PCF for 5 clusters
    %
    % - stable_3: Stable fixed points of the PCF for 3 clusters
    %
    % - unstable_5: Unstable fixed points of the PCF for 5 clusters
    %
    % - unstable_3: Unstable fixed points of the PCF for 3 clusters
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
        Gamma_5;
        Gamma_3;
        stable_5;
        stable_3;
        unstable_5;
        unstable_3;
        order = 7;
        lambda_tol = -1;
        gamma = 1e-2;
        initial = [1;0];
        x_tick = [-1,0,1];
        y_tick = [-1,0,1];
        x_lim = [-1.5,1.5];
        y_lim = [-1.5,1.5];
        name = "sync_cluster";
        evol_basin = 12;
        y_basis = 0;
    end
    methods
        function obj = sync_cluster()
            obj.T = 2*pi;
            obj.dt = obj.T/obj.M;
            obj.t = linspace(0,obj.T-obj.dt,obj.M);
            obj.p = [cos(obj.t);sin(obj.t)];
            obj.Z = [-sin(5*obj.t);2*cos(obj.t)-2*cos(3*obj.t)+cos(5*obj.t)];
            obj.dpdt = [-sin(obj.t);cos(obj.t)];
            obj.dZdt = [-5*cos(5*obj.t);-2*sin(obj.t)+6*sin(3*obj.t)-5*sin(5*obj.t)];
            obj.Gamma_5 = cos(5*obj.t.')/2;
            obj.stable_5 = [1,5,9,13,17] * pi/10;
            obj.unstable_5 = [3,7,11,15,19] * pi/10;
            obj.Gamma_3 = -cos(3*obj.t.');
            obj.stable_3 = [3,7,11] * pi/6;
            obj.unstable_3 = [1,5,9,] * pi/6;
        end
    end
end