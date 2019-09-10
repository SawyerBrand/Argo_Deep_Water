function heat_content_plots(filn,titl)

heat_cont(filn,titl)
diffheat_tot_col(filn,titl)
heat_tot_col(filn,titl)

% % 

    function heat_cont(filn,titl)
        d = load(filn);

        figure(1)
        %% plot most of column
        zmin = 0;
        zmax = 1800;
        dzlayer = 200;
        depths = fliplr(-[zmin:dzlayer:zmax-dzlayer; zmin+dzlayer:dzlayer:zmax]);
        nlayers = size(depths,2);

        heatLayers = nan(nlayers,1);
        cols = jet(nlayers);
        fig1 = figure(1);
        set(fig1,'units','normalized','outerposition',[0 0 1 1])
        hold on
        %lgdtxt = cell(1,nlayers);
        minval = inf;
        maxval = -inf;

        for i = 1:1:size(d.idat.layerHeat)
            heat(i,:) = sum(d.idat.layerHeat(i,1:100));
            x(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
        end
        ONE = plot(x(:,1),heat,'color',cols(1,:),'LineWidth',1.5);

        for i = 1:1:size(d.idat.layerHeat)
            heat(i,:) = sum(d.idat.layerHeat(i,100:200));
            x(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
        end
        TWO = plot(x(:,1),heat,'color',cols(2,:),'LineWidth',1.5);

        for i = 1:1:size(d.idat.layerHeat)
            heat(i,:) = sum(d.idat.layerHeat(i,200:300));
            x(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
        end
        THREE = plot(x(:,1),heat,'color',cols(3,:),'LineWidth',1.5);

        for i = 1:1:size(d.idat.layerHeat)
            heat(i,:) = sum(d.idat.layerHeat(i,300:400));
            x(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
        end
        FOUR = plot(x(:,1),heat,'color',cols(4,:),'LineWidth',1.5);

        for i = 1:1:size(d.idat.layerHeat)
            heat(i,:) = sum(d.idat.layerHeat(i,400:500));
            x(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
        end
        FIVE = plot(x(:,1),heat,'color',cols(5,:),'LineWidth',1.5);

        for i = 1:1:size(d.idat.layerHeat)
            heat(i,:) = sum(d.idat.layerHeat(i,500:600));
            x(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
        end
        SIX = plot(x(:,1),heat,'color',cols(6,:),'LineWidth',1.5);

        for i = 1:1:size(d.idat.layerHeat)
            heat(i,:) = sum(d.idat.layerHeat(i,600:700));
            x(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
        end
        SEVEN = plot(x(:,1),heat,'color',cols(7,:),'LineWidth',1.5);

%         for i = 1:1:size(d.idat.layerHeat)
%             heat(i,:) = sum(d.idat.layerHeat(i,700:800));
%             x(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
%         end
%         EIGHT = plot(x(:,1),heat,'color',cols(8,:),'LineWidth',1.5);
% 
%         for i = 1:1:size(d.idat.layerHeat)
%             heat(i,:) = sum(d.idat.layerHeat(i,800:899));
%             x(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
%         end
%         NINE = plot(x(:,1),heat,'color',cols(9,:),'LineWidth',1.5);


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
                datenum('061019','mmddyy')];
            
        if titl == 'Float 12881'
            xtickangle(45);
            xticks(xdates12881);             %% this uses the list of dates above
        end

        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',13)
        set(gcf,'color','w')

        grid on

        title([titl,': Depth Integrated over 200m Layers'],'fontsize',18)
        xlabel('time','fontsize',18)
        ylabel('depth-integrated heat content [J/m^2]','fontsize',18)
        legend([ONE TWO THREE FOUR FIVE SIX SEVEN],'0-200m','200-400m','400-600m','600-800m','800-1000m','1000-1200m','1200-1400m');
        set(gca,'color','none')


        hold off



        %%
        %%%%
        %%%%%%%

        figure(2)
        %% plot most of column
        zmin = 0;
        zmax = 2400;
        dzlayer = 100;
        depths = fliplr(-[zmin:dzlayer:zmax-dzlayer; zmin+dzlayer:dzlayer:zmax]);
        nlayers = size(depths,2);
        heatLayers = nan(nlayers,1);
        cols = jet(nlayers);
        fig2 = figure(2);
        set(fig2,'units','normalized','outerposition',[0 0 1 1])
        hold on
        %lgdtxt = cell(1,nlayers);
        minval = inf;
        maxval = -inf;

        for i = 1:1:size(d.idat.layerHeat)
            heat(i,:) = sum(d.idat.layerHeat(i,1:50));
            x(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
        end
        One = plot(x(:,1),heat,'r-','LineWidth',1.5);

        for i = 1:1:size(d.idat.layerHeat)
            heat(i,:) = sum(d.idat.layerHeat(i,50:500));
            y(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
        end
        Two = plot(y(:,1),heat,'b-','LineWidth',1.5);

        for i = 1:1:size(d.idat.layerHeat)
            heat(i,:) = sum(d.idat.layerHeat(i,500:899));
            y(i,:) = d.idat.time(i,1);%'color',cols(i,1),'linewidth',1)
        end
        Three = plot(y(:,1),heat,'g-','LineWidth',1.5);


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
                datenum('061019','mmddyy')];

        if titl == 'Float 12881'
            xtickangle(45);
            xticks(xdates12881);             %% this uses the list of dates above
        end

        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',13)

        grid on

        title('Float 12881: Heat Integrated over 3 Layers','fontsize',18)
        xlabel('time','fontsize',18)
        ylabel('depth-integrated heat content [J/m^2]','fontsize',18)
        legend([One Two Three],'0-200m','200-1000m','1000-1800m')
        set(gca,'color','none')


        hold off
        
    end

    function heat_tot_col(fn,titl)
        
        d = load(fn);
        %% plot most of column
        zmin = 0;
        zmax = 1800;
        dzlayer = 200;
        depths = fliplr(-[zmin:dzlayer:zmax-dzlayer; zmin+dzlayer:dzlayer:zmax]);
        nlayers = size(depths,2);

        heatLayers = nan(nlayers,1);
        cols = jet(nlayers);
        fig1 = figure(3);
        set(fig1,'units','normalized','outerposition',[0 0 1 1])
        hold on
        %lgdtxt = cell(1,nlayers);
        minval = inf;
        maxval = -inf;

        for i = 1:1:size(d.idat.layerHeat)
            heat(i,:) = 2*sum(d.idat.layerHeat(i,1:751));
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
        elseif titl == 'Float 12881'
            xtickangle(45);
            xticks(xdates12881);             %% this uses the list of dates above
            datetick('x','dd/mm/yy','keepticks') ;
        end  

        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',13)
        set(gcf,'color','w')

        grid on

        title([titl,' : Depth Integrated down to 1500m'],'fontsize',18)
        xlabel('time','fontsize',18)
        ylabel('depth-integrated heat content [J/m^2]','fontsize',18)
        legend([TOT],'Heat Integration');
        set(gca,'color','none')


        hold off
        
    end

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
        elseif titl == 'Float 12881'
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

end