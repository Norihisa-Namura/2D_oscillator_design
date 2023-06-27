function [lambda_true,Z,dZdt,v1_,u0_,u1_] = floquet(T,dt,init,cla)
    % This function calculates the PSF of two-dimensional limit-cycle oscillators.
    %
    %
    % input: 
    %
    % - T: Period of oscillation
    %
    % - dt: Time interval
    %
    % - init: Initial state for time evolution (column vector)
    %
    % - cla: Class of a limit-cycle oscillator
    %
    % output: 
    %
    % - lambda_true: Second Floquet exponent
    %
    % - Z: PSF
    %
    % - dZdt: Time derivative of PSF

    Tnum = round(T/dt);
    omega = 2*pi/T;

    [X0_, u0_] = Calc_X0_u0(init, dt, Tnum,cla);

    disp("calc PSF")
    [v0_, v0_dif] = Calc_v0(X0_, dt, Tnum,cla); 
    Z = omega*v0_;
    dZdt = omega*v0_dif;

    [lambda_true, u1_] = Calc_u1(X0_, dt, u0_, v0_, T, Tnum,cla);

    v1_ = Calc_v1(X0_, dt, u0_, v0_, u1_, Tnum, lambda_true,cla);

    % check biorthogonality
    prod1 = sum(u0_.*v0_,1);
    prod2 = sum(u1_.*v0_,1);
    prod3 = sum(u0_.*v1_,1);
    prod4 = sum(u1_.*v1_,1);
    figure()
    plot(prod1,'o')
    hold on
    plot(prod2,'o')
    legend("$u_{0} \cdot v_{0}$", "$u_{1} \cdot v_{0}$")
    figure()
    plot(prod3,'o')
    hold on
    plot(prod4,'o')
    legend("$u_{0} \cdot v_{1}$", "$u_{1} \cdot v_{1}$")
    hold off
end

%% functions
% Time evolution for Jacobian matrix 
function y = Jacobi_evol_U(X, dt, X_half_next, X_next, y, lam, cla)
    % 4-th order Runge-Kutta method for Jacobian matrix

    SIZE = size(X,1);
    k1 = dt * (cla.func_jacobi(X(1,:),X(2,:))*y - lam *y);
    k2 = dt * (cla.func_jacobi(X_half_next(1,:),X_half_next(2,:)) - lam*eye(SIZE))*(y+k1/2);
    k3 = dt * (cla.func_jacobi(X_half_next(1,:),X_half_next(2,:)) - lam*eye(SIZE))*(y+k2/2);
    k4 = dt * (cla.func_jacobi(X_next(1,:),X_next(2,:)) - lam*eye(SIZE))*(y+k3);
    y = y + (k1+2*k2+2*k3+k4)/6 ;
end
function y = Jacobi_evol_V(X, dt, X_half_next, X_next, y, lam,cla)
    % 4-th order Runge-Kutta method for transposed Jacobian matrix

    SIZE = size(X,1);
    k1 = dt * (cla.func_jacobi(X(1,:),X(2,:)).'*y - lam *y);
    k2 = dt * (cla.func_jacobi(X_half_next(1,:),X_half_next(2,:)) - lam*eye(SIZE)).'*(y+k1/2);
    k3 = dt * (cla.func_jacobi(X_half_next(1,:),X_half_next(2,:)) - lam*eye(SIZE)).'*(y+k2/2);
    k4 = dt * (cla.func_jacobi(X_next(1,:),X_next(2,:)) - lam*eye(SIZE)).'*(y+k3);
    y = y + (k1+2*k2+2*k3+k4)/6 ;
end

% Calculation for Floquet vectors
% u0
function [X0_, u0_] = Calc_X0_u0(X, dt, Tnum,cla)
    X0_ = zeros(size(X, 1),Tnum);
    u0_ = zeros(size(X, 1),Tnum);
    for tt = 1:Tnum
        X0_(:,tt) = X;
        f = cla.func(X);
        u0_(:,tt) = f;
        X = funcs.runge_kutta_4(X,dt,cla);
    end
end

% v0
function [v0_, v0_dif] = Calc_v0(X0_, dt, Tnum,cla)
    v0_ = zeros(size(X0_, 1),Tnum);
    v0_dif = zeros(size(X0_, 1),Tnum);
    v0 = ones(size(X0_, 1),1); 

    for rep = 1:3
        for tt = 1:Tnum
            [v0,v0dif] = evol_v0(v0,X0_,Tnum,tt,dt,cla);
        end
    end

    % logging
    for tt = 1:Tnum
        v0_(:,Tnum-tt+1) = v0;
        v0_dif(:,Tnum-tt+1) = v0dif;
        [v0,v0dif] = evol_v0(v0,X0_,Tnum,tt,dt,cla);
    end
end
function [v0,v0dif] = evol_v0(v0,X0_,Tnum,tt,dt,cla)
    X = X0_(:,Tnum-tt+1);
    h = -1/2*dt;
    X_half_next = funcs.runge_kutta_4(X,h,cla);
    X_next = X0_(:,mod((Tnum-tt-1),Tnum)+1); %X_next = X(t-dt)
    v0 = Jacobi_evol_V(X, dt, X_half_next, X_next, v0, 0,cla);
    f = cla.func(X_next);
    prob = v0.'* f; % production <v0, u0>
    v0 = v0 / prob;
    
    v0dif = -cla.func_jacobi(X_next(1),X_next(2))'*v0;
end

% u1
function [lambda1, u1_] = Calc_u1(X0_, dt, u0_, v0_, T, Tnum,cla)
    u1_ = zeros(size(X0_, 1),Tnum);
    u1 = ones(size(X0_, 1),1); 
    u1 = u1/norm(u1);
    for rep = 1:5 % for convergence
        for tt = 1:Tnum
            u1 = evol_u1(u1,X0_,u0_,v0_,Tnum,tt,dt,0,cla);
        end
        norm1 = norm(u1);
        u1 = u1/norm1;
    end
    lambda1 = log(norm1)/T;

    % logging
    for tt = 1:Tnum
        u1_(:,tt) = u1;
        u1 = evol_u1(u1,X0_,u0_,v0_,Tnum,tt,dt,lambda1,cla);
    end
end
function u1 = evol_u1(u1,X0_,u0_,v0_,Tnum,tt,dt,lambda1,cla)
    prod = v0_(:,tt).'* u1; % production <v0, u0>
    u1 = u1 - prod*u0_(:,tt);
    X = X0_(:,tt);
    h = 1/2*dt;
    X_half_next = funcs.runge_kutta_4(X,h,cla);
    X_next = X0_(:,mod(tt,Tnum)+1); %X_next = X(t+dt)
    
    u1 = Jacobi_evol_U(X, dt, X_half_next, X_next, u1, lambda1,cla);
end

% v1
function [v1_] = Calc_v1(X0_, dt, u0_, v0_, u1_, Tnum, lambda1,cla)
% calculation of amplitude sensitivity function

    v1_ = zeros(size(X0_, 1),Tnum);
    v1 = ones(size(X0_, 1),1); % initial value of v1
    v1 = v1/norm(v1);

    for rep = 1:3 % for convergence
        for tt = 1:Tnum
            v1 = evol_v1(v1,X0_,u0_,v0_,Tnum,tt,dt,lambda1,cla);
        end
        v1 = v1 / (v1.'*u1_(:,Tnum));
    end

    % logging
    for tt = 1:Tnum
        v1_(:,Tnum-tt+1) = v1;
        v1 = evol_v1(v1,X0_,u0_,v0_,Tnum,tt,dt,lambda1,cla);
    end
end
function v1 = evol_v1(v1,X0_,u0_,v0_,Tnum,tt,dt,lambda1,cla)
% evolution of v1

    prod = u0_(:,Tnum-tt+1).'* v1; % production <v0, u0>
    v1 = v1 - prod*v0_(:,Tnum-tt+1);
    X = X0_(:,Tnum-tt+1);
    h = -1/2*dt;
    X_half_next = funcs.runge_kutta_4(X,h,cla);
    X_next = X0_(:,mod((Tnum-tt-1),Tnum)+1); %X_next = X(t-dt)
    
    v1 = Jacobi_evol_V(X, dt, X_half_next, X_next, v1, lambda1,cla);
end