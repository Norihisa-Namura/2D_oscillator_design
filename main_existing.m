% The code is for designing the vector field of existing two-dimensional limit-cycle oscillators 
% with prescribed trajectories and phase-response characteristics
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

save_ind = 0;

%% Select an oscillator
%class_name = @stuart_landau;
class_name = @van_der_pol;
%class_name = @fitzhugh_nagumo;
%class_name = @CIMA;
cla = class_name();

%% Load data
%load("data/results_" + cla.name); 

%% Parameters and Functions
dt = cla.dt;
[T,omega,initial_tmp,time_lc] = funcs.period(dt,cla);
M = length(time_lc);
order = cla.order;
P = sum(1:order+1);
lambda_tol = cla.lambda_tol;

if cla.name == "sl"
    T = 2*pi;
    omega = 1;
    lambda_true = -2;
    p = [cos(time_lc)';sin(time_lc)'];
    dpdt = [-sin(time_lc)';cos(time_lc)'];
    Z = [-sin(time_lc)' - cla.beta*cos(time_lc)'; cos(time_lc)' - cla.beta*sin(time_lc)'];
    dZdt = [-cos(time_lc)' + cla.beta*sin(time_lc)'; -sin(time_lc)' - cla.beta*cos(time_lc)'];
else
    p = funcs.phase_map(T,dt,initial_tmp,cla);
    dpdt = cla.func(p);
    
    [lambda_true,Z,dZdt] = floquet(T,dt,initial_tmp,cla);
end

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

%% Reconstruction
rep = designed(dt,order,xi,cla);

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

%% Save data for MATLAB
if save_ind == 1
    close all
    save("data/results_" + cla.name);
end
