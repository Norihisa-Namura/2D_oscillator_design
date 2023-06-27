classdef utils
    methods
    end
    methods (Static)
        function [X,X1,X2,area_size] = mesh_grid(cla,n_delta)
            if ~exist("n_delta",'var')
                n_delta = 21;
            end
            x1 = linspace(cla.x_lim(1),cla.x_lim(2),n_delta);
            x2 = linspace(cla.y_lim(1),cla.y_lim(2),n_delta);
            [X1,X2] = meshgrid(x1,x2);
            area_size = size(X1);
            X = horzcat(X1(:),X2(:));
            X = X.';
        end

        function plot_lc(p,cla,p_rep)
            if exist("p_rep",'var')
                plot(p_rep(1,:),p_rep(2,:),'r')
                hold on
                plot(p(1,:),p(2,:),'k--')
            else
                plot(p(1,:),p(2,:),'k')
            end
            xlabel("$x_{1}$")
            ylabel("$x_{2}$")
            xlim(cla.x_lim)
            ylim(cla.y_lim)
            xticks(cla.x_tick)
            yticks(cla.y_tick)
            xtickangle(0)
            axis square
        end

        function plot_val(t,x,y_label,style)
            plot(t,x,style)
            xlabel("$t$")
            ylabel(y_label)
            xlim([0,t(end)])
            xticks([0,t(end)/3,2*t(end)/3,t(end)])
        end

        function plot_comp(t,t_rep,x,x_rep,y_label,fontsize)
            utils.plot_val(t_rep,x_rep,y_label,'r')
            hold on
            utils.plot_val(t,x,y_label,'k--')

            xlim([0,max(t(end),t_rep(end))])
            range = (max([x,x_rep]) - min([x,x_rep])) / 6;
            ylim([min([x,x_rep]) - range,max([x,x_rep]) + range])
            tmp = max(t(end),t_rep(end));
            xticks([0,tmp/3,2*tmp/3,tmp])
            xticklabels(["$0$","$\frac{1}{3}T$","$\frac{2}{3}T$","$T$"])
            xlabel("$t$",'FontSize',fontsize)
            ylabel(y_label,'FontSize',fontsize)
            xtickangle(0)
        end

        function fig = show_results(t,t_rep,p,dpdt,p_rep,dpdt_rep,Z,Z_rep,cla)
            fig = figure();
            fig.Position(3:4) = [700,800];
            str = ["(a)","(b)","(c)","(d)","(e)"];
            fontsize = 20;
            
            pos = [0.2,0.7,0.60,0.28];
            str_pos12 = pos(1:2) + [0.275,-0.065];
            str_pos = [str_pos12,0,0];
            subplot('Position',pos);
            utils.plot_lc(p,cla,p_rep);
            annotation('textbox',str_pos,'String',str(1),'EdgeColor','none','FontSize',fontsize)
            
            for i = 1:2
                pos = [0.47*(i-1)+0.1,0.4,0.35,0.18];
                str_pos12 = pos(1:2) + [0.152,-0.055];
                str_pos = [str_pos12,0,0];
                subplot('Position',pos);

                y_label = "$\dot{p}_{" + num2str(i) + "}$";
                utils.plot_comp(t,t_rep,dpdt(i,:),dpdt_rep(i,:),y_label,fontsize)
                if cla.name == "fhn" && i == 2
                    ax = gca;
                    ax.YAxis.Exponent = -3;
                end
                annotation('textbox',str_pos,'String',str(i+1),'EdgeColor','none','FontSize',fontsize)
            end
            
            for i = 1:2
                pos = [0.47*(i-1)+0.1,0.09,0.35,0.2];
                str_pos12 = pos(1:2) + [0.152,-0.055];
                str_pos = [str_pos12,0,0];
                subplot('Position',pos);

                y_label = "$\tilde{Z}_{" + num2str(i) + "}$";
                utils.plot_comp(t,t_rep,Z(i,:),Z_rep(i,:),y_label,fontsize);
                annotation('textbox',str_pos,'String',str(i+3),'EdgeColor','none','FontSize',fontsize)
            end
        end
        
        function plot_oscillators(count,p,x,cla)
            plot(p(1,:),p(2,:),'k','LineWidth',1.5)
            hold on
            if count == 3 || count == 4
                scatter(x(1,:),x(2,:),400,'r.')
            else
                scatter(x(1,:),x(2,:),200,'r.')
            end
            axis equal
            xlim(cla.x_lim)
            ylim(cla.y_lim)
            xticks(cla.x_tick)
            yticks(cla.y_tick)
            xlabel("$x_{1}$")
            ylabel("$x_{2}$")
            xtickangle(0)
        end

        function fig = show_oscillators(p_rep,x_plot,cla)
            fig = figure();
            fig.Position(3:4) = [700,700];
            str = ["(a)","(b)","(c)","(d)"];
            for count = 1:4
                pos = [rem(count+1,2)*0.45 + 0.1,1 - fix((count+1)/2)*1/2 + 0.12,0.35,0.35];
                str_pos12 = pos(1:2) + [0.15,-0.075];
                str_pos = [str_pos12,0,0];
                subplot('Position',pos);
                utils.plot_oscillators(count,p_rep,squeeze(x_plot(:,:,count)),cla)
                annotation('textbox',str_pos,'String',str(count),'EdgeColor','none','FontSize',20)
            end
        end

        function fig = show_pcf(phi,Gamma,cla,stable,label,unstable)
            ymax = max(Gamma);
            ymin = min(Gamma);
            
            if cla.name == "sync_cluster"
                fig = [];
            else
                fig = figure();
            end

            yline(0,'k--','LineWidth',2)
            hold on
            plot(phi,Gamma,'r')
            hold on
            if ~exist("stable",'var')
                for fp = cla.stable
                    scatter(fp,0,1000,'k.')
                    hold on
                    %plot(fp*ones(1,100),y,'--k')
                    xline(fp,'k--','LineWidth',2)
                    hold on
                end
            else
                for fp = stable
                    scatter(fp,0,1000,'k.')
                    hold on
                    %plot(fp*ones(1,100),y,'--k')
                    xline(fp,'k--','LineWidth',2)
                    hold on
                end
            end
            if ~exist("stable",'var')
                unstable = cla.unstable;
            end
            for fp = unstable
                plot(fp,0,'ok')
                hold on
                %plot(fp*ones(1,100),y,'--k')
                %hold on
            end

            box on

            if cla.name == "sync_cluster"
                xticks([0,pi/2,pi,3*pi/2,2*pi])
                xticklabels(["$0$","$\frac{1}{2}\pi$","$\pi$","$\frac{3}{2}\pi$","$2\pi$"])
                xlim([0,2*pi])
                ylim([ymin,ymax])
                xlabel("$\phi$")
                xtickangle(0)
                ylabel(label)
            else
                xticks([0,pi/2,pi,3*pi/2,2*pi])
                xticklabels(["$0$","$\frac{1}{2}\pi$","$\pi$","$\frac{3}{2}\pi$","$2\pi$"])
                xlim([0,2*pi])
                ylim([ymin,ymax])
                xlabel("$\phi$")
                ylabel("$\Gamma(\phi)$")
            end
        end
        
        function save_figure(save_ind,fig,save_name,cla)
            if save_ind == 1
                exportgraphics(fig,"figs/" + cla.name + "_" + save_name + ".pdf")
            end
        end

        function output_data(p_rep,cla,rep)
            [X,X1,X2,area_size] = utils.mesh_grid(cla);
            
            F_rep = rep.func(X);
            F1_rep = reshape(F_rep(1,:),area_size);
            F2_rep = reshape(F_rep(2,:),area_size);

            if ismethod(cla,"func")
                F_true = cla.func(X);
                F1_true = reshape(F_true(1,:),area_size);
                F2_true = reshape(F_true(2,:),area_size);
            else
                F1_true = zeros(size(F1_rep));
                F2_true = zeros(size(F2_rep));
            end
            
            utils.write_data(X1)
            utils.write_data(X2)
            utils.write_data(F1_true)
            utils.write_data(F2_true)
            utils.write_data(F1_rep)
            utils.write_data(F2_rep)
            utils.write_data(p_rep)
            
            x_tick = rep.x_tick;
            y_tick = rep.y_tick;
            utils.write_data(x_tick)
            utils.write_data(y_tick)
            writematrix(rep.name,"data/class_name.txt");
        end

        function write_data(data)
            savename = inputname(1);
            writematrix(data,"data/" + savename + ".dat");
        end
    end
end