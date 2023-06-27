classdef funcs
    methods
    end
    methods (Static)
        function U = polynomial(x,order)
            n = 1:order+1;
            l = length(x(1,:));
            
            U = zeros(sum(n),l);
            for i = 1:order+1
                for j = 1:n(i)
                    U(sum(n(1:i)) -j+1,:) = x(1,:).^(j-1) .* x(2,:).^(n(i)-j); % powered n(i)-1
                end
            end
        end
        
        function [dUdx] = polynomial_diff(x,order)
            % # of parameter * dimension * # of data (3rd order tensor)
            n = 1:order+1;
            num = length(x(1,:));
            dUdx = zeros(sum(n),2,num);
            for i = 1:order+1
                for j = 1:n(i)
                    if j < 2
                        dUdx(sum(n(1:i)) -j+1,1,:) = 0;
                    else
                        dUdx(sum(n(1:i)) -j+1,1,:) = (j-1) * x(1,:).^(j-2) .* x(2,:).^(n(i)-j); % powered n(i)-1
                    end
                    if n(i) - j - 1 < 0
                        dUdx(sum(n(1:i)) -j+1,2,:) = 0;
                    else
                        dUdx(sum(n(1:i)) -j+1,2,:) = x(1,:).^(j-1) .* ((n(i) - j) * x(2,:).^(n(i)-j-1)); % powered n(i)-1
                    end
                end
            end
        end
        
        function [M] = issemiposdef(M)
            if issymmetric(M)
                d = eig(M);
                tf = all(d >= 0);
                disp(tf)
                if tf == 0
                    val = 10 * min(d);
                    disp(val)
                    M = M - val*eye(size(M));
                end
            else 
                disp("Matrix is not sym")
            end
        end

        function [T,omega,initial_tmp,time_lc] = period(dt,cla)
            % calculate period from data
            %dt = 0.001;
            M = 1e8; % max iteration
            xx = zeros(length(cla.initial),3);
            xx(:,2) = cla.initial;
            xx(:,3) = funcs.runge_kutta_4(xx(:,2),dt,cla);
            count = 0;
            icount = 0;
            for i = 1:M-1
                xx(:,1:2) = xx(:,2:3);
                xx(:,3) = funcs.runge_kutta_4(xx(:,2),dt,cla);
                cond_var1 = (xx(2,1) - cla.y_basis) * (xx(2,2) - cla.y_basis);
                if i > icount + 30 && cond_var1 < 0
                    count = count + 1;
                    icount = i;
                    disp(xx(:,3).')
                    if count == 28
                        i_start = i;
                    elseif count == 30 
                        T = (i - i_start) * dt;
                    	break;
                    end
                else
                    T = 0;
                end
            end
            if T == 0
                T = Inf;
            end
            omega = 2*pi/T;
            initial_tmp = xx(:,3);
            if omega ~= 0
                time_lc = (0:dt:T-dt).';
            else
                time_lc = 0;
            end
        end
        
        function [x_lc,theta_lc] = phase_map(T,dt,initial,cla)
            % phase map on limit cycle
            omega = 2*pi/T;
            n = round(T/dt); % max iteration
            x_lc = zeros(length(initial),n);
            theta_lc = zeros(n,1);
            
            x_lc(:,1) = initial;
            
            for i = 1:n-1
                x_lc(:,i+1) = funcs.runge_kutta_4(x_lc(:,i),dt,cla);
                theta_lc(i+1) = omega*dt*i;
            end
        end
        
        function [xi] = optimizer(A,b,Q,A_ineq,b_ineq,gamma,sigma_A,sigma_b)
            % proposed convex optimization (quadratic programming)
            AA = funcs.issemiposdef(A.'*A);
            funcs.issemiposdef(AA + gamma*Q);
            [xi,resnorm] = quadprog(AA + gamma*Q,-A.'*b,A_ineq,b_ineq);
            
            disp("residual")
            disp(resnorm - gamma*norm(xi)^2/2 + b.'*b/2)
            
            xi = xi * sigma_b ./ sigma_A';
        end
        
        function [A,b,Q,A_ineq,b_ineq,sigma_A,sigma_b] = mat_2D(Ux,dUdx,dpdt,Z,dZdt,M,P,lambda_tol)
            % calculate matrices for optimization
            A1 = horzcat(Ux.',zeros(M,P));
            A2 = horzcat(zeros(M,P),Ux.');

            % jacobi matrix
            dUdx1 = squeeze(dUdx(:,1,:));
            dUdx2 = squeeze(dUdx(:,2,:));
            A3 = horzcat((Z(1,:).'*ones(1,P)).*dUdx1.',(Z(2,:).'*ones(1,P)).*dUdx1.');
            A4 = horzcat((Z(1,:).'*ones(1,P)).*dUdx2.',(Z(2,:).'*ones(1,P)).*dUdx2.');

            Q = eye(2*P);

            % lambda < 0
            A_ineq = [sum(dUdx1,2).',sum(dUdx2,2).']/M;
            
            b_ineq = lambda_tol;

            A = vertcat(A1,A2,A3,A4);
            b = vertcat(dpdt(1,:).',dpdt(2,:).',-dZdt(1,:).',-dZdt(2,:).'); 

            % standardization
            [~,mu_A,sigma_A] = zscore(A);
            [~,mu_b,sigma_b] = zscore(b);
            A = (A - mu_A) ./ sigma_A;
            b = (b - mu_b) ./ sigma_b;
            A_ineq = (A_ineq - mu_A) ./ sigma_A;
            b_ineq = (b_ineq - mu_b) ./ sigma_b;
        end
        
        function fx = runge_kutta_4(x,dt,cla)
            ka = cla.func(x);
            kb = cla.func(x + dt*ka/2);
            kc = cla.func(x + dt*kb/2);
            kd = cla.func(x + dt*kc);
            fx = x + dt*(ka + 2*kb + 2*kc + kd)/6;
        end
        
        function fx = runge_kutta_4_input(x,u,dt,cla)
            % u: N * 3
            ka = cla.func(x);
            ka = ka + u(:,1);
            kb = cla.func(x + dt*ka/2);
            kb = kb + u(:,2);
            kc = cla.func(x + dt*kb/2);
            kc = kc + u(:,2);
            kd = cla.func(x + dt*kc);
            kd = kd + u(:,3);
            fx = x + dt*(ka + 2*kb + 2*kc + kd)/6;
        end
    end
end