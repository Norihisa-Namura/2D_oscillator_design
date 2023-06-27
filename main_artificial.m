% The code is for designing the vector field of artificial two-dimensional limit-cycle oscillators 
% with prescribed trajectories and phase-response characteristics and
% also for entrainmenting limit-cycke oscillators to periodic inputs
%
% If you use this code, please cite the following paper:
%
%   XXXXXXXXXXXXXXXXXXXXXXX
%
%
% Norihisa Namura (2023)
% 
% namura.n.aa@m.titech.ac.jp

%% Setting
clear
close all
mystyle
rng('default')

save_ind = 1;

%% Select an oscillator
cla = star();
%cla = sync_cluster();

%% Load data
%load("data/results_" + cla.name);

%% Parameters and functions
T = cla.T;
dt = cla.dt;
M = cla.M;

p = cla.p;
Z = cla.Z;
dpdt = cla.dpdt;
dZdt = cla.dZdt;

omega = 2*pi/T;
time_lc = cla.t.';

order = cla.order;
P = sum(1:order+1);
lambda_tol = cla.lambda_tol;

%% Normalization
p_norm = mean(vecnorm(dpdt,2,1));
pZ_rate = mean(vecnorm(dZdt,2,1)) / p_norm;

Z_norm = Z / pZ_rate;
dZdt_norm = dZdt / pZ_rate;

%% Design the vector field
U_period = funcs.polynomial(p,order);
dUdx = funcs.polynomial_diff(p,order);
[A,b,Q,A_ineq,b_ineq,sigma_A,sigma_b] = funcs.mat_2D(U_period,dUdx,dpdt,Z_norm,dZdt_norm,M,P,lambda_tol);
gamma = cla.gamma;
[xi] = funcs.optimizer(A,b,Q,A_ineq,b_ineq,gamma,sigma_A,sigma_b);

%% Designed oscillator
initial_tmp = cla.initial;
rep = designed(dt,order,xi(1:2*P),cla);
[T_rep,omega_rep,~,time_lc_rep] = funcs.period(dt,rep);
disp(abs(T_rep - T) / T)

if omega_rep ~= 0
    p_rep = funcs.phase_map(T_rep,dt,initial_tmp,rep);
    dpdt_rep = rep.func(p_rep);
    [lambda_rep,Z_rep] = floquet(T_rep,dt,initial_tmp,rep);
end

%% Output data for drawing vector field by python
utils.output_data(p_rep,cla,rep)

%% Show results 
fig = utils.show_results(time_lc,time_lc_rep,p,dpdt,p_rep,dpdt_rep,Z,Z_rep,cla);
utils.save_figure(save_ind,fig,"results",rep);

%% Entrainment
disp("entrainment")

if cla.name == "sync_cluster"
    theta = 2*pi*rand(1,100);  
    x = [cos(theta);sin(theta)];
    x_plot = zeros(size(x,1),size(x,2),4);
    count = 1;
    x_plot(:,:,count) = x;

    epsilon = 5e-2;

    num_five = 16;
    num_three = 8;
    num = num_five + num_three;
    
    for i = 1:round(num*2*pi/dt)
        t = i*dt;
        tt = [t,t+dt,t+2*dt];
        
        if i < round(num_five*2*pi/dt)
            Omega = 5*omega;
            u = [zeros(size(tt));cos(Omega*tt)];
            x = funcs.runge_kutta_4_input(x,epsilon*u,dt,rep);
        else
            Omega = 3*omega;
            u = [zeros(size(tt));cos(Omega*tt)];
            x = funcs.runge_kutta_4_input(x,epsilon*u,dt,rep);
        end
        
        if ismember(i,round([3*2*pi,num_five*2*pi,num*2*pi]/dt))
            count = count + 1;
            x_plot(:,:,count) = x;
        end
    end
elseif cla.name == "star"
    theta = 2*pi*rand(1,100);  
    x = [cla.r1*cos(theta) + cla.r2*sin(cla.m*theta); cla.r1*sin(theta) + cla.r2*cos(cla.m*theta)]; % initial state
    x_plot = zeros(size(x,1),size(x,2),4);
    count = 1;
    x_plot(:,:,count) = x;

    epsilon = 1e-2;

    for i = 1:round(100*2*pi/dt)
        t = i*dt;
        tt = [t,t+dt,t+2*dt];
        
        u = [-sin(omega*tt);cos(omega*tt)];
        x = funcs.runge_kutta_4_input(x,epsilon*u,dt,rep);
        
        if ismember(i,round([20*2*pi,40*2*pi,100*2*pi]/dt))
            count = count + 1;
            x_plot(:,:,count) = x;
        end
    end
end

%% Show oscillators
fig = utils.show_oscillators(p_rep,x_plot,cla);
utils.save_figure(save_ind,fig,"oscillators",rep);

%% Show PCF
if cla.name == "star"
    fig = utils.show_pcf(time_lc,cla.Gamma,cla);
    utils.save_figure(save_ind,fig,"pcf",rep);
elseif cla.name == "sync_cluster"
    fig = figure();
    fig.Position(3:4) = [700,300];

    pos = [0.1,0.25,0.35,0.7];
    subplot('Position',pos)
    utils.show_pcf(time_lc,cla.Gamma_5,cla,cla.stable_5,"$\Gamma(\phi;5\Omega)$",cla.unstable_5);
    title("(a)", 'Units', 'normalized', 'Position', [0.5, -0.36, 0],'Interpreter','none');
    box on

    pos = [0.6,0.25,0.35,0.7];
    subplot('Position',pos)
    utils.show_pcf(time_lc,cla.Gamma_3,cla,cla.stable_3,"$\Gamma(\phi;3\Omega)$",cla.unstable_3);
    title("(b)", 'Units', 'normalized', 'Position', [0.5, -0.36, 0],'Interpreter','none');
    box on
    utils.save_figure(save_ind,fig,"pcf",rep);
end

%% Save data for MATLAB
if save_ind == 1
    close all
    save("data/results_" + cla.name);
end
