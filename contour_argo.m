function contour_argo(filenam,titl)
        d = load(filenam);
        
        figure(7)
        
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
        hold on
        
        set(gca,'ydir','reverse')
        set(gca, 'FontSize',18)
        
        ax = gca;
        ax.TickDir = 'out';
        
        title([titl ': Conservative Temp Contour'])
        ylabel('Pressure [dbars]');
        xlabel('Profiles');
        contourf((d.idat.time),d.idat.P,d.idat.ct);
        
        pressure_levels = [0,4,6,8,9,10];
        
        h = colorbar;
        
        [M,h2] = contour((d.idat.time),d.idat.P,d.idat.ct,'k','ShowText','on');
        h2.LineWidth = 0.5;
        clabel(M,h2,'FontSize',9,'Color','k','FontName','Ubuntu')
        
        xtickangle(45);
        xticks(xdates12881);             %% this uses the list of dates above
        datetick('x','dd/mm/yy','keepticks') ;  
        
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',13)
        
        %         dateFormat = 'dd mmm yyyy';
        %         datetick('x',dateFormat)
        
        hold off
        
        figure(8)
        
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
            
        hold on
        
        set(gca,'ydir','reverse')
        set(gca,'FontSize',18)
        
        ax = gca;
        ax.TickDir = 'out';
        
        title([titl ': Heat Contour'])
        ylabel('Pressure [dbars]');
        xlabel('Profiles');
        
        colormap('spring')
        h = colorbar;
        contourf((d.idat.time),d.idat.P,d.idat.heat);
        
        pressure_levels = [0,4,6,8,9,10];
        
        [M,h2] = contour((d.idat.time),d.idat.P,d.idat.heat,'k','ShowText','on');
        h2.LineWidth = 0.5;
        clabel(M,h2,'FontSize',9,'Color','k','FontName','Ubuntu')
            
        xtickangle(45);
        xticks(xdates12881);             %% this uses the list of dates above
        datetick('x','dd/mm/yy','keepticks') ;       %% this tells it to really use the list you gave it

        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',13)

        hold off 
end 
       