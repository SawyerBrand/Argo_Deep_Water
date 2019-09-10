function pCO2_plots(fn,titl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was made by Sawyer Brand, SIO/SOCCOM            %
% during research into the Argentine Basin mesoscale activities %
% and CO2/Heat Flux in that region                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The functions that are used in this mfile are called below:
%pCO2_profile(fn,titl)
%TALK_profile(fn,titl)
%DIC_profile(fn,titl)
%pCO2_Contour_Plots(fn,titl)

%function to display basic pCO2 profiles for each float:

    function pCO2_profile(fn,titl)
        
        f = load(fn); 
        
        ColOrd = get(gca,'ColorOrder');
        % Determine the number of colors in
        % the matrix
        [m,n] = size(ColOrd);

        
        hold on
        %subplot(1,2,1)
        title([titl,' : pCO2 profiles'],'FontSize',20)
        set(gca,'ydir','reverse')
        ylabel('Pressure [dbar]','FontSize',17)
        xlabel('pCO2 [mmHg]','FontSize',17)
        
        
        for i = 1:1:size(f.idat.pCO2,1)
            ColRow = rem(i,m);
            if ColRow == 0
              ColRow = m;
            end
            % Get the color
            Col = ColOrd(ColRow,:);
            plot(f.idat.pCO2(i,:),f.idat.P(i,:),'linewidth', 1, 'Color',Col)
        end

       

        hold off
    end


%%
function TALK_profile(fn,titl)
        
        f = load(fn); 
        
        ColOrd = get(gca,'ColorOrder');
        % Determine the number of colors in
        % the matrix
        [m,n] = size(ColOrd);

        figure(2)
        hold on
        %subplot(1,2,1)
        title([titl,' : Alkalinity profiles'],'FontSize',20)
        set(gca,'ydir','reverse')
        ylabel('Pressure [dbar]','FontSize',17)
        xlabel('Alkalinity [micro-mol kg/ soln]','FontSize',17)
        
        
        for i = 1:1:size(f.idat.TALK,1)
            ColRow = rem(i,m);
            if ColRow == 0
              ColRow = m;
            end
            % Get the color
            Col = ColOrd(ColRow,:);
            plot(f.idat.TALK(i,:),f.idat.P(i,:),'linewidth', 1, 'Color',Col)
        end

       

        hold off
end

%%
function DIC_profile(fn,titl)
        
        f = load(fn); 
        figure(3)
        ColOrd = get(gca,'ColorOrder');
        % Determine the number of colors in
        % the matrix
        [m,n] = size(ColOrd);

        
        hold on
        %subplot(1,2,1)
        title([titl,' : Dissolved Inorganic Carbon profiles'],'FontSize',20)
        set(gca,'ydir','reverse')
        ylabel('Pressure [dbar]','FontSize',17)
        xlabel('DIC [micro-meters]','FontSize',17)
        
        
        for i = 1:1:size(f.idat.DIC,1)
            ColRow = rem(i,m);
            if ColRow == 0
              ColRow = m;
            end
            % Get the color
            Col = ColOrd(ColRow,:);
            plot(f.idat.DIC(i,:),f.idat.P(i,:),'linewidth', 1, 'Color',Col)
        end

       

        hold off
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%CONTOUR BELOW%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 
    function pCO2_Contour_Plots(fn,titl)
        
        d = load(fn);
        
        figure(2)
        
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
            
        hold on
        
        set(gca,'ydir','reverse')
        ax = gca;
        ax.TickDir = 'out';
        
        contourf((d.idat.time),d.idat.P,d.idat.pCO2);
        
        pressure_levels = [0,4,6,8,9,10];
        
        h = colorbar;
        
        [M,h2] = contour((d.idat.time),d.idat.P,d.idat.pCO2,'k','ShowText','on');
        h2.LineWidth = 0.75;
        clabel(M,h2,'FontSize',11,'Color','k','FontWeight','bold')
        
        if titl == 'Float 12700'
            plot((d.idat.time(17,:)),d.idat.P(17,:),'k-','LineWidth',2)
            plot((d.idat.time(6,:)),d.idat.P(6,:),'k-','LineWidth',2)
        else 
            plot((d.idat.time(21,:)),d.idat.P(21,:),'k-','LineWidth',2)
        end
                
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
        
        title([titl ': pCO2 Contour'],'FontSize',22)
        ylabel('Pressure [dbars]','FontSize',17);
        xlabel('Profiles','FontSize',17);
        
        %         dateFormat = 'dd mmm yyyy';
        %         datetick('x',dateFormat)
        
        hold off

    end







end