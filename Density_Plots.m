function Density_Plots(fn,titl)

        dd = load(fn);
        
        zmin = 0;
        zmax = (size(dd.ipdat.T)).^2;
        dz = size(dd.ipdat.T);
        depths = fliplr(-[zmin:dz:zmax-dz; zmin+dz:dz:zmax]);
        nlayers = size(depths,2);
        cols = hsv(nlayers);
        
        figure(1)
        
        pos1 = [0.08 0.1 0.4 0.85];
        subplot('Position',pos1)
        
        hold on
        grid on
        set(gca,'ydir','reverse')
        
        title([titl,' : Temperature Profiles'],'FontSize',20)
        xlabel('Temp [deg C]','FontSize',18)
        ylabel('Pot Density [kg/m^3]','FontSize',18)
        
        xmin = 0;
        xmax = (size(dd.ipdat.T)).^2;
        dx = size(dd.ipdat.T);
        depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
        nlayers = size(depths,2);
        cols = hsv(nlayers);
        
        for i = 1:size(dd.ipdat.T)
            plot(dd.ipdat.T(i,:),dd.ipdat.pden(i,:),'color',cols(i,:))
        end
        
        if titl == '12700'
            plot(dd.ipdat.T(17,:),dd.ipdat.pden(17,:),'r.-','MarkerSize',10)
        end
        
        if titl == '12778'
            plot(dd.ipdat.T(6,:),dd.ipdat.pden(6,:),'r.-','MarkerSize',10)
        end
        
        hold off
        
        %saveas(gcf,[titl,'_T_pden_Profiles.jpg'])
        
        pos2 = [0.54 0.1 0.45 0.85];
        %left bottom width height
        subplot('Position',pos2)
        
        hold on
        grid on
        set(gca,'ydir','reverse')
        
        title([titl,' : Absolute Salinity Profiles'],'FontSize',20)
        xlabel('Abs Sal [g/kg]','FontSize',18)
        ylabel('Pot Density [kg/m^3]','FontSize',18)
        
        xmin = 0;
        xmax = (size(dd.ipdat.sa)).^2;
        dx = size(dd.ipdat.sa);
        depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
        nlayers = size(depths,2);
        cols = hsv(nlayers);
        
        for i = 1:size(dd.idat.sa)
            plot(dd.ipdat.sa(i,:),dd.ipdat.pden(i,:),'color',cols(i,:))
        end
        
        if titl == '12700'
            plot(dd.ipdat.sa(17,:),dd.ipdat.pden(17,:),'r.-','MarkerSize',10)
        end
        
        if titl == '12778'
            plot(dd.ipdat.sa(6,:),dd.ipdat.pden(6,:),'r.-','MarkerSize',10)
        end
       
        saveas(gcf,[titl,'_SA_CT_pden_Profiles.jpg']) 
        hold off
        
        figure(2)
        hold on
        grid on
        set(gca,'ydir','reverse')
        
        title([titl,': Dissolved Oxygen Profiles in Argentine Basin'],'FontSize',20)
        xlabel('Diss Oxy [micro-mol/kg]','FontSize',18)
        ylabel('Pot Density [kg/m^3]','FontSize',18)
        
        xmin = 0;
        xmax = (size(dd.ipdat.DO)).^2;
        dx = size(dd.ipdat.DO);
        depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
        nlayers = size(depths,2);
        cols = hsv(nlayers);
        
        for i = 1:size(dd.idat.DO)
            plot(dd.ipdat.DO(i,:),dd.ipdat.pden(i,:),'color',cols(i,:))
        end
        
        if titl == '12700'
            plot(dd.ipdat.DO(17,:),dd.ipdat.pden(17,:),'r.-','MarkerSize',10)
        end
        
        if titl == '12778'
            plot(dd.ipdat.DO(6,:),dd.ipdat.pden(6,:),'r.-','MarkerSize',10)
        end
      
        saveas(gcf,[titl,'_DO_pden_Profiles.jpg'])

        hold off


end