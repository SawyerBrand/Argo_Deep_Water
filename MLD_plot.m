function MLD_plot(fn1,titl1,fn2,fn3,fn4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was made by Sawyer Brand, SIO/SOCCOM            %
% during research into the Argentine Basin mesoscale activities %
% and CO2/Heat Flux in that region                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MLD_plots_waterfall(fn1,titl1)
%MLD_plots_basic(fn1,titl1)
%MLD_plots_comparison(fn1,fn2,fn3,fn4)
    

    function MLD_plots_waterfall(fn,titl)
            f = load(fn);


            for i = 1:1:size(f.idat.rho,1)
                rho = f.idat.rho(i,:) - 1025;
                RhoII = rho(~isnan(rho));
                RhoTop = RhoII(:,1);

                newP = f.idat.P(i,:);
                newP = newP(find(rho >= (RhoTop + 0.3)));
                MLDval = newP(:,1);
                MLDdepth(i,:) = MLDval;

                H = f.idat.heat(i,:);
                H = H(~isnan(H));
                HEATI = find(rho >= (RhoTop + 0.3));
                HEATI = HEATI(:,1);

                HEATSum(i,:) = sum(H(:,1:HEATI));
            end

            for i = 1:1:size(f.idat.T,1)
                T = f.idat.T(i,:);
                TII = T(~isnan(T));
                TTop = TII(:,1);

                newP = f.idat.P(i,:);
                newP = newP(find(T <= (TTop - 0.4)));
                MLDval = newP(:,1);
                MLDdepthT(i,:) = MLDval;

                H = f.idat.heat(i,:);
                H = H(~isnan(H));
                HEATI = find(T <= (TTop - 0.4));
                HEATI = HEATI(:,1);

                HEATSumT(i,:) = sum(H(:,1:HEATI));;
            end

            figure
            
            pos1 = [0.1 0.55 0.85 0.42];
            %left bottom width height
            subplot('Position',pos1)
            
            hold on
            set(gca,'Fontsize',17)
            title(['Temperature vs Pressure ',titl])
            xlabel('In-Situ Temp (C) + 13')     %using in-situ for this with a delta T of 15 (eyeballed value)
            ylabel('In-situ Pressure')
            
            set(gca,'Ydir','reverse')       %Just because of depth
            
            for i = 1:size(f.idat.T,1)
                plot((f.idat.T(i,:)+0.75*(i)),f.idat.P(i,:),'b')    %adds 15 to the index each time, ends up being a good way to increment
            end
            hold off
            
            pos2 = [0.1 0.08 0.85 0.4];
            %left bottom width height
            subplot('Position',pos2)
            
            hold on
            grid on 

            set(gca,'ydir','reverse')
            title([titl,': Heat Mixed Layer Depth Analysis'],'FontSize',20)
            ylabel('Heat []','FontSize',18)
            xlabel('Profile Dates','FontSize',18)

            plot(f.idat.time,HEATSum,'r*-','LineWidth',1.5)
            plot(f.idat.time,HEATSumT,'b*-','LineWidth',1.5)

            xdates12881 = [datenum('102218','mmddyy'),datenum('110118','mmddyy'),...   %%%% I wanted monthly date ticks, but you could make a vector of
                datenum('111118','mmddyy'),datenum('112118','mmddyy'),...                  %%%% dates that correspond only to float profile days
                datenum('120118','mmddyy'),datenum('121118','mmddyy'),datenum('122118','mmddyy'),...
                datenum('123118','mmddyy'),datenum('011019','mmddyy'),...
                datenum('012019','mmddyy'),datenum('013019','mmddyy'),...
                datenum('020919','mmddyy'),datenum('021919','mmddyy'),...
                datenum('030119','mmddyy'),datenum('031219','mmddyy'),...
                datenum('032219','mmddyy'),datenum('040119','mmddyy'),...
                datenum('041119','mmddyy'),datenum('042119','mmddyy'),...
                datenum('050119','mmddyy'),datenum('051119','mmddyy'),...
                datenum('052119','mmddyy'),datenum('053119','mmddyy'),...
                datenum('061019','mmddyy'),datenum('062019','mmddyy'),...
                datenum('063019','mmddyy'),datenum('071019','mmddyy'),...
                datenum('072019','mmddyy')];

            xdates12700 = [datenum('102218','mmddyy'),datenum('110118','mmddyy'),...   %%%% I wanted monthly date ticks, but you could make a vector of
                datenum('111118','mmddyy'),datenum('112118','mmddyy'),...                  %%%% dates that correspond only to float profile days
                datenum('120118','mmddyy'),datenum('121118','mmddyy'),datenum('122118','mmddyy'),...
                datenum('123118','mmddyy'),datenum('011019','mmddyy'),...
                datenum('012019','mmddyy'),datenum('013019','mmddyy'),...
                datenum('020919','mmddyy'),datenum('021919','mmddyy'),...
                datenum('030119','mmddyy'),datenum('031219','mmddyy'),...
                datenum('032219','mmddyy'),datenum('040119','mmddyy'),...
                datenum('041119','mmddyy'),datenum('042119','mmddyy'),...
                datenum('050119','mmddyy'),datenum('051119','mmddyy'),...
                datenum('052119','mmddyy'),datenum('053119','mmddyy'),...
                datenum('061019','mmddyy'),datenum('062019','mmddyy'),...
                datenum('063019','mmddyy'),datenum('071019','mmddyy'),...
                datenum('072019','mmddyy')];

        %     if titl == 'Float 12700'
        %         plot((f.idat.time(17,:)),f.idat.P(17,:),'k-','LineWidth',2);
        %     end

            if titl == 'Float 12700'
                xtickangle(45);
                xticks(xdates12700);             %% this uses the list of dates above
                datetick('x','dd/mm/yy','keepticks') ;
            elseif titl == 'Float 12881'
                xtickangle(45);
                xticks(xdates12881);             %% this uses the list of dates above
                datetick('x','dd/mm/yy','keepticks') ;
            end

            set(gca,'ydir','reverse')
            ax = gca;
            ax.TickDir = 'out';
            ax.FontSize = 18;
            
            saveas(gcf,[titl,'_MLD_Plot_Waterfall.jpg'])
    end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% New Plot Starts %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function MLD_plots_basic(fn,titl)
            f = load(fn);


            for i = 1:1:size(f.idat.rho,1)
                rho = f.idat.rho(i,:) - 1025;
                RhoII = rho(~isnan(rho));
                RhoTop = RhoII(:,1);

                newP = f.idat.P(i,:);
                newP = newP(find(rho >= (RhoTop + 0.3)));
                MLDval = newP(:,1);
                MLDdepth(i,:) = MLDval;

                H = f.idat.heat(i,:);
                H = H(~isnan(H));
                HEATI = find(rho >= (RhoTop + 0.3));
                HEATI = HEATI(:,1);

                HEATSum(i,:) = sum(H(:,1:HEATI));
            end

            for i = 1:1:size(f.idat.T,1)
                T = f.idat.T(i,:);
                TII = T(~isnan(T));
                TTop = TII(:,1);

                newP = f.idat.P(i,:);
                newP = newP(find(T <= (TTop - 0.4)));
                MLDval = newP(:,1);
                MLDdepthT(i,:) = MLDval;

                H = f.idat.heat(i,:);
                H = H(~isnan(H));
                HEATI = find(T <= (TTop - 0.4));
                HEATI = HEATI(:,1);

                HEATSumT(i,:) = sum(H(:,1:HEATI));;
            end


            figure(1)
            hold on
            grid on

            set(gca,'ydir','reverse')
            title([titl,': Mixed Layer Depth Analysis'],'FontSize',20)
            ylabel('Pressure [dbars]','FontSize',18)
            xlabel('Profile Dates','FontSize',18)

            a1 = plot(f.idat.time,MLDdepth,'r*-','LineWidth',1.5);
            legend({'0.03 pden threshold'},'Location','northeast')
            
            a2 = plot(f.idat.time,MLDdepthT,'b*-','LineWidth',1.5);
            legend({'0.5 C threshold'},'Location','northeast')

            xdates12881 = [datenum('102218','mmddyy'),datenum('110118','mmddyy'),...   %%%% I wanted monthly date ticks, but you could make a vector of
                datenum('111118','mmddyy'),datenum('112118','mmddyy'),...                  %%%% dates that correspond only to float profile days
                datenum('120118','mmddyy'),datenum('121118','mmddyy'),datenum('122118','mmddyy'),...
                datenum('123118','mmddyy'),datenum('011019','mmddyy'),...
                datenum('012019','mmddyy'),datenum('013019','mmddyy'),...
                datenum('020919','mmddyy'),datenum('021919','mmddyy'),...
                datenum('030119','mmddyy'),datenum('031219','mmddyy'),...
                datenum('032219','mmddyy'),datenum('040119','mmddyy'),...
                datenum('041119','mmddyy'),datenum('042119','mmddyy'),...
                datenum('050119','mmddyy'),datenum('051119','mmddyy'),...
                datenum('052119','mmddyy'),datenum('053119','mmddyy'),...
                datenum('061019','mmddyy'),datenum('062019','mmddyy'),...
                datenum('063019','mmddyy'),datenum('071019','mmddyy'),...
                datenum('072019','mmddyy')];

            xdates12700 = [datenum('102218','mmddyy'),datenum('110118','mmddyy'),...   %%%% I wanted monthly date ticks, but you could make a vector of
                datenum('111118','mmddyy'),datenum('112118','mmddyy'),...                  %%%% dates that correspond only to float profile days
                datenum('120118','mmddyy'),datenum('121118','mmddyy'),datenum('122118','mmddyy'),...
                datenum('123118','mmddyy'),datenum('011019','mmddyy'),...
                datenum('012019','mmddyy'),datenum('013019','mmddyy'),...
                datenum('020919','mmddyy'),datenum('021919','mmddyy'),...
                datenum('030119','mmddyy'),datenum('031219','mmddyy'),...
                datenum('032219','mmddyy'),datenum('040119','mmddyy'),...
                datenum('041119','mmddyy'),datenum('042119','mmddyy'),...
                datenum('050119','mmddyy'),datenum('051119','mmddyy'),...
                datenum('052119','mmddyy'),datenum('053119','mmddyy'),...
                datenum('061019','mmddyy'),datenum('062019','mmddyy'),...
                datenum('063019','mmddyy'),datenum('071019','mmddyy'),...
                datenum('072019','mmddyy')];

        %     if titl == 'Float 12700'
        %         plot((f.idat.time(17,:)),f.idat.P(17,:),'k-','LineWidth',2);
        %     end

            if titl == 'Float 12700'
                xtickangle(45);
                xticks(xdates12700);             %% this uses the list of dates above
                datetick('x','dd/mm/yy','keepticks') ;
            elseif titl == 'Float 12881'
                xtickangle(45);
                xticks(xdates12881);             %% this uses the list of dates above
                datetick('x','dd/mm/yy','keepticks') ;
            end

            set(gca,'ydir','reverse')
            ax = gca;
            ax.TickDir = 'out';
            ax.FontSize = 18;

            legend({'0.03 pden threshold','0.5 C threshold'},'Location','northeast')
            hold off

            figure(2)

            hold on 
            grid on 

            set(gca,'ydir','reverse')
            title([titl,': Heat Mixed Layer Depth Analysis'],'FontSize',20)
            ylabel('Heat []','FontSize',18)
            xlabel('Profile Dates','FontSize',18)

            plot(f.idat.time,HEATSum,'r*-','LineWidth',1.5)
            plot(f.idat.time,HEATSumT,'b*-','LineWidth',1.5)

            xdates12881 = [datenum('102218','mmddyy'),datenum('110118','mmddyy'),...   %%%% I wanted monthly date ticks, but you could make a vector of
                datenum('111118','mmddyy'),datenum('112118','mmddyy'),...                  %%%% dates that correspond only to float profile days
                datenum('120118','mmddyy'),datenum('121118','mmddyy'),datenum('122118','mmddyy'),...
                datenum('123118','mmddyy'),datenum('011019','mmddyy'),...
                datenum('012019','mmddyy'),datenum('013019','mmddyy'),...
                datenum('020919','mmddyy'),datenum('021919','mmddyy'),...
                datenum('030119','mmddyy'),datenum('031219','mmddyy'),...
                datenum('032219','mmddyy'),datenum('040119','mmddyy'),...
                datenum('041119','mmddyy'),datenum('042119','mmddyy'),...
                datenum('050119','mmddyy'),datenum('051119','mmddyy'),...
                datenum('052119','mmddyy'),datenum('053119','mmddyy'),...
                datenum('061019','mmddyy'),datenum('062019','mmddyy'),...
                datenum('063019','mmddyy'),datenum('071019','mmddyy'),...
                datenum('072019','mmddyy')];

            xdates12700 = [datenum('102218','mmddyy'),datenum('110118','mmddyy'),...   %%%% I wanted monthly date ticks, but you could make a vector of
                datenum('111118','mmddyy'),datenum('112118','mmddyy'),...                  %%%% dates that correspond only to float profile days
                datenum('120118','mmddyy'),datenum('121118','mmddyy'),datenum('122118','mmddyy'),...
                datenum('123118','mmddyy'),datenum('011019','mmddyy'),...
                datenum('012019','mmddyy'),datenum('013019','mmddyy'),...
                datenum('020919','mmddyy'),datenum('021919','mmddyy'),...
                datenum('030119','mmddyy'),datenum('031219','mmddyy'),...
                datenum('032219','mmddyy'),datenum('040119','mmddyy'),...
                datenum('041119','mmddyy'),datenum('042119','mmddyy'),...
                datenum('050119','mmddyy'),datenum('051119','mmddyy'),...
                datenum('052119','mmddyy'),datenum('053119','mmddyy'),...
                datenum('061019','mmddyy'),datenum('062019','mmddyy'),...
                datenum('063019','mmddyy'),datenum('071019','mmddyy'),...
                datenum('072019','mmddyy')];

        %     if titl == 'Float 12700'
        %         plot((f.idat.time(17,:)),f.idat.P(17,:),'k-','LineWidth',2);
        %     end

            if titl == 'Float 12700'
                xtickangle(45);
                xticks(xdates12700);             %% this uses the list of dates above
                datetick('x','dd/mm/yy','keepticks') ;
            elseif titl == 'Float 12881'
                xtickangle(45);
                xticks(xdates12881);             %% this uses the list of dates above
                datetick('x','dd/mm/yy','keepticks') ;
            end

            set(gca,'ydir','reverse')
            ax = gca;
            ax.TickDir = 'out';
            ax.FontSize = 18;
    end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% New Plot Starts %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function MLD_plots_comparison(fn1,fn2,fn3,fn4)
        f1 = load(fn1);
        f2 = load(fn2);
        f3 = load(fn3);
        f4 = load(fn4);


        for i = 1:1:size(f1.idat.T,1)
            T = f1.idat.T(i,:);
            TII = T(~isnan(T));
            TTop = TII(:,1);

            newP = f1.idat.P(i,:);
            newP = newP(find(T <= (TTop - 0.4)));
            MLDval = newP(:,1);
            MLDdepthT1(i,:) = MLDval;

            H = f1.idat.heat(i,:);
            H = H(~isnan(H));
            HEATI = find(T <= (TTop - 0.4));
            HEATI = HEATI(:,1);

            HEATSumT1(i,:) = sum(H(:,1:HEATI));;
        end
        
        for i = 1:1:size(f2.idat.T,1)
            T = f2.idat.T(i,:);
            TII = T(~isnan(T));
            TTop = TII(:,1);

            newP = f2.idat.P(i,:);
            newP = newP(find(T <= (TTop - 0.4)));
            MLDval = newP(:,1);
            MLDdepthT2(i,:) = MLDval;

            H = f2.idat.heat(i,:);
            H = H(~isnan(H));
            HEATI = find(T <= (TTop - 0.4));
            HEATI = HEATI(:,1);

            HEATSumT2(i,:) = sum(H(:,1:HEATI));;
        end
        
        for i = 1:1:size(f3.idat.T,1)
            T = f3.idat.T(i,:);
            TII = T(~isnan(T));
            TTop = TII(:,1);

            newP = f3.idat.P(i,:);
            newP = newP(find(T <= (TTop - 0.4)));
            MLDval = newP(:,1);
            MLDdepthT3(i,:) = MLDval;

            H = f3.idat.heat(i,:);
            H = H(~isnan(H));
            HEATI = find(T <= (TTop - 0.4));
            HEATI = HEATI(:,1);

            HEATSumT3(i,:) = sum(H(:,1:HEATI));;
        end
        
        for i = 1:1:size(f4.idat.T,1)
            T = f4.idat.T(i,:);
            TII = T(~isnan(T));
            TTop = TII(:,1);

            newP = f4.idat.P(i,:);
            newP = newP(find(T <= (TTop - 0.4)));
            MLDval = newP(:,1);
            MLDdepthT4(i,:) = MLDval;

            H = f4.idat.heat(i,:);
            H = H(~isnan(H));
            HEATI = find(T <= (TTop - 0.4));
            HEATI = HEATI(:,1);

            HEATSumT4(i,:) = sum(H(:,1:HEATI));;
        end
        
        
        hold on 
        grid on
        
        set(gca,'ydir','reverse')
        title('Mixed Layer Depth Comparison of Argo in Argentine Basin','FontSize',20)
        ylabel('Pressure [dbars]')
        xlabel('Profiles')
        
        a1 = plot(f1.idat.time,MLDdepthT1,'g.-','LineWidth',1.5);
        a2 = plot(f2.idat.time,MLDdepthT2,'r.-','LineWidth',1.5);
        a3 = plot(f3.idat.time(101:122),MLDdepthT3(101:122),'b.-','LineWidth',1.5);
        a4 = plot(f4.idat.time,MLDdepthT4,'k.-','LineWidth',1.5);
        
        legend('Float 12700','Float 12881','Float 09646','Float 12778','Location','southwest');
        
        
        xdates12881 = [datenum('102218','mmddyy'),datenum('110118','mmddyy'),...   %%%% I wanted monthly date ticks, but you could make a vector of
                datenum('111118','mmddyy'),datenum('112118','mmddyy'),...                  %%%% dates that correspond only to float profile days
                datenum('120118','mmddyy'),datenum('121118','mmddyy'),datenum('122118','mmddyy'),...
                datenum('123118','mmddyy'),datenum('011019','mmddyy'),...
                datenum('012019','mmddyy'),datenum('013019','mmddyy'),...
                datenum('020919','mmddyy'),datenum('021919','mmddyy'),...
                datenum('030119','mmddyy'),datenum('031219','mmddyy'),...
                datenum('032219','mmddyy'),datenum('040119','mmddyy'),...
                datenum('041119','mmddyy'),datenum('042119','mmddyy'),...
                datenum('050119','mmddyy'),datenum('051119','mmddyy'),...
                datenum('052119','mmddyy'),datenum('053119','mmddyy'),...
                datenum('061019','mmddyy'),datenum('062019','mmddyy'),...
                datenum('063019','mmddyy'),datenum('071019','mmddyy'),...
                datenum('072019','mmddyy')];

            xdates12700 = [datenum('102218','mmddyy'),datenum('110118','mmddyy'),...   %%%% I wanted monthly date ticks, but you could make a vector of
                datenum('111118','mmddyy'),datenum('112118','mmddyy'),...                  %%%% dates that correspond only to float profile days
                datenum('120118','mmddyy'),datenum('121118','mmddyy'),datenum('122118','mmddyy'),...
                datenum('123118','mmddyy'),datenum('011019','mmddyy'),...
                datenum('012019','mmddyy'),datenum('013019','mmddyy'),...
                datenum('020919','mmddyy'),datenum('021919','mmddyy'),...
                datenum('030119','mmddyy'),datenum('031219','mmddyy'),...
                datenum('032219','mmddyy'),datenum('040119','mmddyy'),...
                datenum('041119','mmddyy'),datenum('042119','mmddyy'),...
                datenum('050119','mmddyy'),datenum('051119','mmddyy'),...
                datenum('052119','mmddyy'),datenum('053119','mmddyy'),...
                datenum('061019','mmddyy'),datenum('062019','mmddyy'),...
                datenum('063019','mmddyy'),datenum('071019','mmddyy'),...
                datenum('072019','mmddyy')];

        %     if titl == 'Float 12700'
        %         plot((f.idat.time(17,:)),f.idat.P(17,:),'k-','LineWidth',2);
        %     end

            if fn2 == '12881.mat'
                xtickangle(45);
                xticks(xdates12881);             %% this uses the list of dates above
                datetick('x','dd/mm/yy','keepticks') ;
            end

            set(gca,'ydir','reverse')
            ax = gca;
            ax.TickDir = 'out';
            ax.FontSize = 18;
            
            hold off
            
    end
end