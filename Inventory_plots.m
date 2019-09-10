function Inventory_plots(filename, titl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was made by Sawyer Brand, SIO/SOCCOM            %
% during research into the Argentine Basin mesoscale activities %
% and CO2/Heat Flux in that region                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%diffheat_tot_col(filename,titl)
DIC_Inventory_Plots(filename,titl)
    

%%
    
    function diffheat_tot_col(fn,titl)
        
        d = load(fn);
        %% plot most of column
        zmin = 0;
        zmax = 1800;
        dzlayer = 200;
        depths = fliplr(-[zmin:dzlayer:zmax-dzlayer; zmin+dzlayer:dzlayer:zmax]);
        nlayers = size(depths,2);

        heatLayers = nan(nlayers,1);
        cols = jet(nlayers);
        fig1 = figure(4);
        set(fig1,'units','normalized','outerposition',[0 0 1 1])
        hold on
        %lgdtxt = cell(1,nlayers);
        minval = inf;
        maxval = -inf;

        for i = 2:1:size(d.idat.layerHeat)
            heat(i,:) = 2*sum(d.idat.layerHeat(i,1:751)-d.idat.layerHeat(i-1,1:751));
            x(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
        end
        TOT = plot(x(:,1),heat,'color',cols(1,:),'LineWidth',1.5);

        yl = [minval - 0.1*(maxval-minval) maxval + 0.1*(maxval-minval)];

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
                datenum('063019','mmddyy'),...
                datenum('071019','mmddyy'),datenum('072019','mmddyy'),...
                datenum('073019','mmddyy')];
            
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
                datenum('063019','mmddyy'),...
                datenum('071019','mmddyy'),datenum('072019','mmddyy'),...
                datenum('073019','mmddyy')];

        if titl == 'Float 12700'
            xtickangle(45);
            xticks(xdates12700);             %% this uses the list of dates above
            datetick('x','dd/mm/yy','keepticks') ;
        else
            xtickangle(45);
            xticks(xdates12881);             %% this uses the list of dates above
            datetick('x','dd/mm/yy','keepticks') ;
        end  

        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',13)
        set(gcf,'color','w')

        grid on

        title([titl,' : Difference in Heat Integrated down to 1500m'],'fontsize',18)
        xlabel('time','fontsize',18)
        ylabel('depth-integrated heat content [J/m^2]','fontsize',18)
        legend([TOT],'Heat Integration');
        set(gca,'color','none')


        hold off
        
    end

%%
    function DIC_Inventory_Plots(fn,titl)
            d = load(fn);
            zmin = 0;
            zmax = 1800;
            dzlayer = 200;
            depths = fliplr(-[zmin:dzlayer:zmax-dzlayer; zmin+dzlayer:dzlayer:zmax]);
            nlayers = size(depths,2);

            heatLayers = nan(nlayers,1);
            cols = jet(nlayers);
            fig1 = figure(5);
            set(fig1,'units','normalized','outerposition',[0 0 1 1])
            hold on
            %lgdtxt = cell(1,nlayers);
            minval = inf;
            maxval = -inf;

            newDIC = d.idat.DIC.*d.idat.rho;
            
            for i = 1:1:size(newDIC,1)
                goodDIC = newDIC(i,10:750);
                DIC(i,:) = 2*sum(goodDIC);
                x(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
            end
            
            TOT = plot(d.idat.time(:,1),DIC,'LineWidth',1.5);
            
            if titl == 'Float 12700'
                plot((d.idat.time(17,1)),DIC(17,:),'k-')
            end

            yl = [minval - 0.1*(maxval-minval) maxval + 0.1*(maxval-minval)];

            xdates12881 = [datenum('102118','mmddyy'),datenum('103118','mmddyy'),...   %%%% I wanted monthly date ticks, but you could make a vector of
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
                datenum('063019','mmddyy'),...
                datenum('071019','mmddyy'),datenum('072019','mmddyy'),...
                datenum('073019','mmddyy')];

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
                datenum('052119','mmddyy'),datenum('053119','mmddyy')];

            if titl == 'Float 12700'
                xtickangle(45);
                xticks(xdates12700);             %% this uses the list of dates above
                datetick('x','dd/mm/yy','keepticks') ;
            else
                xtickangle(45);
                xticks(xdates12881);             %% this uses the list of dates above
                datetick('x','dd/mm/yy','keepticks') ;
            end
            

            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'fontsize',13)
            set(gcf,'color','w')

            grid on

            title([titl,' : DIC Integrated down to 1500m'],'fontsize',18)
            xlabel('time','fontsize',18)
            ylabel('depth-integrated heat content [J/m^2]','fontsize',18)
            %legend([TOT],'Heat Integration');
            set(gca,'color','none')


            hold off
    end 



end