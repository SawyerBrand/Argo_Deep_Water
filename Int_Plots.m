function Int_Plots(fn,titl)

big = process_argo(fn,titl);
%basic_int_plots([titl,'.mat'],titl)
%contour_argo([titl,'.mat'],titl)
%contour_pden_argo[titl,'.mat'],titl)
%heat_diff_profiles[titl,'.mat'],titl)


%% get the interpolated data:

    function data = process_argo(fn,titl)
  
        outfile = [titl,'.mat'];
        addpath(genpath('/Volumes/SOCCOM/Old_Research/GSW/Toolbox...'));

        % Read the variables in the netcdf (.nc) file and save in struct
        ncid = netcdf.open(fn,'NC_NOWRITE');
        [~, numvars, ~, ~] = netcdf.inq(ncid);
        for i = 1:numvars
            varname = netcdf.inqVar(ncid,i-1);
            raw_data.(varname) = netcdf.getVar(ncid,i-1);
        end

        data = struct();
        data.raw = raw_data; %retain raw data for future referece
        t = datenum('1950-01-01 00:00:00');
%         t = raw_data.REFERENCE_DATE_TIME';
%         ref_t = [t(1:4) '-' t(5:6) '-' t(7:8) ' ' t(9:10) ':' t(11:12) ':' t(13:14)];
        data.time = (raw_data.JULD+t);

        data.lat  = raw_data.Lat;
        data.lon = raw_data.Lon;
%         for i =1:size(raw_data.Lon)
%             x = raw_data.Lon(i,:);
%             data.lonA = x(find(x > 180))-360;
%             data.lonB = x(find(x <= 180));
%             data.lon(:,i) = [data.lonA data.lonB];
%         end
%         data.lon = data.lon';
        
%         data.lon(1:(x/2))  = raw_data.Lon(1:(x/2))-360;
        data.P  = fliplr(abs(raw_data.Pressure)'); %chose to use the nonadjusted vairables because they had all values of 99999
        data.T  = fliplr(raw_data.Temperature');
        data.SP  = fliplr(raw_data.Salinity');
        data.DO  = fliplr(raw_data.Oxygen');
        data.pCO2 = raw_data.pCO2_LIAR';
        data.TALK = raw_data.TALK_LIAR';
        data.DIC = raw_data.DIC_LIAR';
        data.N2 = raw_data.Nitrate;
        data.b_bp700 = raw_data.b_bp700';
        data.Chl_a = raw_data.Chl_a';


        % missing data
        data.P(data.P >= 9000 | data.P <=-1000) = NaN;
        data.T(data.T >= 9000 | data.T <= -1000) = NaN;
        data.SP(data.SP >= 9000 | data.SP <= -1000) = NaN;
        data.DO(data.DO >= 9000 | data.DO <= -1000) = NaN;
        data.pCO2(data.pCO2 >= 9000 | data.pCO2 <= -1000) = NaN;
        data.TALK(data.TALK >= 9000 | data.TALK <= -1000) = NaN;
        data.DIC(data.DIC >= 9000 | data.DIC <= -1000) = NaN;
        data.DIC(isnan(data.DIC)) = 0;
        %data.N2(data.N2 >= 9000 | data.N2 <= -1000) = NaN;
        data.b_bp700(data.b_bp700 >= 9000 | data.b_bp700 <= -1000) = NaN;
        data.Chl_a(data.Chl_a >= 9000 | data.Chl_a <= -1000) = NaN;
        
        if titl == '12700'
            zint = 1:1:1900;
        else
            zint = 1:2:2000;
        end

        
          % depth
          data.z = gsw_z_from_p(data.P, data.lat);

          % absolute salinity
          [data.sa, data.sstar, data.in_ocean] = gsw_SA_Sstar_from_SP(data.SP, data.P, data.lon, data.lat);

          % potential temperature
          data.pt = gsw_pt_from_t(data.sa, data.T, data.P, 0);

          % conservative temperature
          data.ct = gsw_CT_from_t(data.sa, data.T, data.P);

          % potential density
          data.pden = gsw_rho(data.sa, data.ct, 0);

          % in situ density
          data.rho = gsw_rho(data.sa, data.ct, data.P);

          % brunt-väisälä frequency
          [data.N2, data.p_mid] = gsw_Nsquared(data.sa, data.ct, data.P, data.lat);
          
          data.cp = gsw_cp_t_exact(data.sa,data.ct, data.P);
          
          data.heat = data.rho .* data.cp .* (data.ct+273.15);
          
          [C,idx]=unique(data.T,'stable');%don't sort
          idx=setxor(idx,1:numel(data.T));
          data.T(idx) = NaN;
          
          [C,idx]=unique(data.P,'stable');%don't sort
          idx=setxor(idx,1:numel(data.P));
          data.P(idx) = NaN;
          


       for i = 1:1:size(data.T,1)
          tmp = data.P(i,:);
          tmpT = data.T(i,:);
          tmpSP = data.SP(i,:);
          tmpDO = data.DO(i,:);
          tmpSA = data.sa(i,:);
          tmpPT = data.pt(i,:);
          tmpCT = data.ct(i,:);
          tmpPD = data.pden(i,:);
          tmpR = data.rho(i,:);
          tmpH = data.heat(i,:);
          tmppCO2 = data.pCO2(i,:);
          tmpTALK = data.TALK(i,:);
          tmpDIC = data.DIC(i,:);
          tmpb_bp700 = data.b_bp700(i,:);
          tmpChl = data.Chl_a(i,:);
          
          I = find(isfinite(tmp)==1 & isfinite(tmpT)==1);
          idat.T(i,:) = interp1(tmp(I),tmpT(I),zint);

          I2 = find(isfinite(tmp)==1 & isfinite(tmpSP)==1);
          idat.SP(i,:) = interp1(tmp(I2),tmpSP(I2),zint);
          
          I4 = find(isfinite(tmp)==1 & isfinite(tmp)==1);
          idat.P(i,:) = interp1(tmp(I4),tmp(I4),zint);  
          
          I5 = find(isfinite(tmp)==1 & isfinite(tmpSA)==1);
          idat.sa(i,:) = interp1(tmp(I5),tmpSA(I5),zint);
          
          I6 = find(isfinite(tmp)==1 & isfinite(tmpPT)==1);
          idat.pt(i,:) = interp1(tmp(I6),tmpPT(I6),zint);
          
          I7 = find(isfinite(tmp)==1 & isfinite(tmpCT)==1);
          idat.ct(i,:) = interp1(tmp(I7),tmpCT(I7),zint);
          
          I8 = find(isfinite(tmp)==1 & isfinite(tmpPD)==1);
          idat.pden(i,:) = interp1(tmp(I8),tmpPD(I8),zint);
          
          I9 = find(isfinite(tmp)==1 & isfinite(tmpR)==1);
          idat.rho(i,:) = interp1(tmp(I9),tmpR(I9),zint);
          
          I10 = find(isfinite(tmp)==1 & isfinite(tmpH)==1);
          idat.heat(i,:) = interp1(tmp(I10),tmpH(I10),zint); 
          
%           I11 = find(isfinite(tmp)==1 & isfinite(tmppCO2)==1);
%           idat.pCO2(i,:) = interp1(tmp(I11),tmppCO2(I11),zint);
%           
%           I12 = find(isfinite(tmp)==1 & isfinite(tmpTALK)==1);
%           idat.TALK(i,:) = interp1(tmp(I12),tmpTALK(I12),zint);
%           
%           I13 = find(isfinite(tmp)==1 & isfinite(tmpDIC)==1);
%           idat.DIC(i,:) = interp1(tmp(I13),tmpDIC(I13),zint);

         I15 = find(isfinite(tmp)==1 & isfinite(tmpb_bp700)==1);
         idat.b_bp700(i,:) = interp1(tmp(I15),tmpb_bp700(I15),zint);
         
         I16 = find(isfinite(tmp)==1 & isfinite(tmpChl)==1);
         idat.Chl(i,:) = interp1(tmp(I16),tmpChl(I16),zint);
%           
          
       end
       
       for i = 1:size(data.DO,1)
          tmp= data.P(i,:);
          tmpDO = data.DO(i,:);
          I3 = find(isfinite(tmp)==1 & isfinite(tmpDO)==1);
          idat.DO(i,:) = interp1(tmp(I3),tmpDO(I3),zint);
       end
       
       for i = 1:size(data.N2,1)
          tmp= data.P(i,:);
          tmpN2 = data.N2(i,:);
          I3 = find(isfinite(tmp)==1 & isfinite(tmpN2)==1);
          idat.N2(i,:) = interp1(tmp(I3),tmpN2(I3),zint);
       end
       
       for i = 1:1:size(idat.ct,2)
           idat.time(:,i) = (data.time);
       end
       
       %%% PDEN 
       
       if titl == '12700'
            pint = 1027.45:0.001:1027.7;
       elseif titl == '12778'
            pint = 1027.62:0.001:1027.75;
       else 
            pint = 1027:0.001:1028;
       end
        
       pint = 1027.6:0.001:1027.8;
       
       for i = 1:size(data.T,1)
          tmp= data.pden(i,:);
          tmpT = data.T(i,:);
          tmpSP = data.SP(i,:);
          tmpDO = data.DO(i,:);
          tmpSA = data.sa(i,:);
          tmpPT = data.pt(i,:);
          tmpCT = data.ct(i,:);
          tmpPD = data.P(i,:);
          tmpR = data.rho(i,:);
          tmpH = data.heat(i,:);
          
          I = find(isfinite(tmp)==1 & isfinite(tmpT)==1);
          ipdat.T(i,:) = interp1(tmp(I),tmpT(I),pint);
          
          I2 = find(isfinite(tmp)==1 & isfinite(tmpSP)==1);
          ipdat.SP(i,:) = interp1(tmp(I2),tmpSP(I2),pint);
          
          I4 = find(isfinite(tmp)==1 & isfinite(tmp)==1);
          ipdat.pden(i,:) = interp1(tmp(I4),tmp(I4),pint);  
          
          I5 = find(isfinite(tmp)==1 & isfinite(tmpSA)==1);
          ipdat.sa(i,:) = interp1(tmp(I5),tmpSA(I5),pint);
          
          I6 = find(isfinite(tmp)==1 & isfinite(tmpPT)==1);
          ipdat.pt(i,:) = interp1(tmp(I6),tmpPT(I6),pint);
          
          I7 = find(isfinite(tmp)==1 & isfinite(tmpCT)==1);
          ipdat.ct(i,:) = interp1(tmp(I7),tmpCT(I7),pint);
          
          I8 = find(isfinite(tmp)==1 & isfinite(tmpPD)==1);
          ipdat.P(i,:) = interp1(tmp(I8),tmpPD(I8),pint);
          
          I9 = find(isfinite(tmp)==1 & isfinite(tmpR)==1);
          ipdat.rho(i,:) = interp1(tmp(I9),tmpR(I9),pint);
          
          I10 = find(isfinite(tmp)==1 & isfinite(tmpH)==1);
          ipdat.heat(i,:) = interp1(tmp(I10),tmpH(I10),pint);       
          
       end
       
       for i = 1:size(data.DO,1)
          tmp= data.pden(i,:);
          tmpDO = data.DO(i,:);
          I3 = find(isfinite(tmp)==1 & isfinite(tmpSP)==1);
          ipdat.DO(i,:) = interp1(tmp(I3),tmpDO(I3),pint);
       end
       
       for i = 1:1:size(ipdat.ct,2)
           ipdat.time(:,i) = (data.time);
       end
       
       dz = 100;
       
       layermeanSA = (idat.sa(:,1:end-1) + idat.sa(:,2:end)) ./ 2;
       layermeanP = (idat.P(:,1:end-1) + idat.P(:,2:end)) ./ 2;
       layermeanT = (idat.T(:,1:end-1) + idat.T(:,2:end)) ./ 2;
       layermeanRho = (idat.rho(:,1:end-1) + idat.rho(:,2:end)) ./ 2;
       layermeanCT = (idat.ct(:,1:end-1) + idat.ct(:,2:end)) ./ 2;
       idat.layermeanSP = (idat.SP(:,1:end-1) + idat.SP(:,2:end)) ./ 2;
       idat.layerSP = layermeanRho .* idat.layermeanSP .* dz ./ 1000; % [kg]
       idat.cp = gsw_cp_t_exact(layermeanSA,layermeanT,layermeanP);
       idat.layerHeat = layermeanRho .* idat.cp .* (layermeanCT+273.15) .* dz; % [J / m^2]

       idat.layerHeat(isnan(idat.layerHeat)) = 0;
    
%           keyboard
%               if length(I)>2
%                   idat.T(i,:) = interp1(data.P(I(i,:)),data.T(I(i,:)),zint)
%               end
%               if length(I2)>2
%                   idat.SP(i,:) = interp1(data.P(I(i,:)),data.SP(I(i,:)),zint)
%               end


        save(outfile,'idat','data','raw_data','zint','ipdat');

    end
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% 
    function basic_int_plots(filenam,titl)
        dd = load(filenam);
        
        zmin = 0;
        zmax = (size(dd.idat.T)).^2;
        dz = size(dd.idat.T);
        depths = fliplr(-[zmin:dz:zmax-dz; zmin+dz:dz:zmax]);
        nlayers = size(depths,2);
        cols = hsv(nlayers);
        
        figure
        hold on
        grid on
        set(gca,'ydir','reverse')
        
        title([titl,' : Temperature Profiles'],'FontSize',20)
        xlabel('Temp [deg C]','FontSize',18)
        ylabel('Pressure [dbars]','FontSize',18)
        
        xmin = 0;
        xmax = (size(dd.idat.T)).^2;
        dx = size(dd.idat.T);
        depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
        nlayers = size(depths,2);
        cols = hsv(nlayers);
        
        for i = 1:size(dd.idat.T)
            plot(dd.idat.T(i,:),dd.idat.P(i,:),'color',cols(i,:))
        end
        
        if titl == '12700'
            plot(dd.idat.T(17,:),dd.idat.P(17,:),'k.')
        end
        
        if titl == '12778'
            plot(dd.idat.T(6,:),dd.idat.P(6,:),'k.')
        end
        
        saveas(gcf,[titl,'_T_Profiles.jpg'])
        hold off
        
        figure
        hold on
        grid on
        set(gca,'ydir','reverse')
        
        title([titl,' : Absolute Salinity Profiles'],'FontSize',20)
        xlabel('Abs Sal [g/kg]','FontSize',18)
        ylabel('Pressure [dbars]','FontSize',18)
        
        xmin = 0;
        xmax = (size(dd.idat.sa)).^2;
        dx = size(dd.idat.sa);
        depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
        nlayers = size(depths,2);
        cols = hsv(nlayers);
        
        for i = 1:size(dd.idat.sa)
            plot(dd.idat.sa(i,:),dd.idat.P(i,:),'color',cols(i,:))
        end
        
        if titl == '12700'
            plot(dd.idat.sa(17,:),dd.idat.P(17,:),'k.')
        end
        
        if titl == '12778'
            plot(dd.idat.sa(6,:),dd.idat.P(6,:),'k.')
        end
       
        saveas(gcf,[titl,'_SA_Profiles.jpg']) 
        hold off
        
        figure
        hold on
        grid on
        set(gca,'ydir','reverse')
        
        title([titl,': Dissolved Oxygen Profiles in Argentine Basin'],'FontSize',20)
        xlabel('Diss Oxy [micro-mol/kg]','FontSize',18)
        ylabel('Pressure [dbars]','FontSize',18)
        
        xmin = 0;
        xmax = (size(dd.idat.DO)).^2;
        dx = size(dd.idat.DO);
        depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
        nlayers = size(depths,2);
        cols = hsv(nlayers);
        
        for i = 1:size(dd.idat.DO)
            plot(dd.idat.DO(i,:),dd.idat.P(i,:),'color',cols(i,:))
        end
        
        if titl == '12700'
            plot(dd.idat.DO(17,:),dd.idat.P(17,:),'k.')
        end
        
        if titl == '12778'
            plot(dd.idat.DO(6,:),dd.idat.P(6,:),'k.')
        end
      
        saveas(gcf,[titl,'_DO_Profiles.jpg'])

        hold off
        
%         figure
%         hold on
%         grid on
%         set(gca,'ydir','reverse')
%         
%         title([titl,': Nitrate Profiles in Argentine Basin'],'FontSize',20)
%         xlabel('Diss Oxy [micro-mol/kg]','FontSize',18)
%         ylabel('Pressure [dbars]','FontSize',18)
%         
%         xmin = 0;
%         xmax = (size(dd.data.N2)).^2;
%         dx = size(dd.data.N2);
%         depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
%         nlayers = size(depths,2);
%         cols = hsv(nlayers);
%         
%         for i = 1:size(dd.raw_data.Nitrate)
%             plot(dd.raw_data.Nitrate(i,:),dd.raw_data.Pressure(i,:),'color',cols(i,:))
%         end
%         
%         if titl == '12700'
%             plot(dd.raw_data.Nitrate(17,:),dd.raw_data.Pressure(17,:),'k.')
%         end
%         
%         if titl == '12778'
%             plot(dd.raw_data.Nitrate(6,:),dd.raw_data.Pressure(6,:),'k.')
%         end
% 
% %         dec20th = plot(dd.idat.DO(106,:),dd.idat.P(106,:),'r-','LineWidth',2)
% %         
% %         mar2nd = plot(dd.idat.DO(113,:),dd.idat.P(113,:),'k-','LineWidth',1.5)
% %         legend([dec20th mar2nd],'Profile taken on Dec 20th','Profile March 2nd');
%         
%         saveas(gcf,[titl,'_N2_Profiles.jpg'])
% 
%         hold off
%         
        
        figure
        hold on
        grid on
        set(gca,'ydir','reverse')
        
        title([titl,': Flouresence Profiles in Argentine Basin'],'FontSize',20)
        xlabel('b_bp700 [micro-mol/kg]','FontSize',18)
        ylabel('Pressure [dbars]','FontSize',18)
        
        xmin = 0;
        xmax = (size(dd.idat.b_bp700)).^2;
        dx = size(dd.idat.b_bp700);
        depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
        nlayers = size(depths,2);
        cols = hsv(nlayers);
        
        for i = 1:size(dd.idat.b_bp700)
            plot(dd.idat.b_bp700(i,:),dd.idat.P(i,:),'color',cols(i,:))
        end
        
        if titl == '12700'
            plot(dd.idat.b_bp700(17,:),dd.idat.P(17,:),'k.')
        end
        
        if titl == '12778'
            plot(dd.idat.b_bp700(6,:),dd.idat.P(6,:),'k.')
        end
      
        saveas(gcf,[titl,'_b_bp700_Profiles.jpg'])

        hold off
        
        figure
        hold on
        grid on
        set(gca,'ydir','reverse')
        
        title([titl,': Chlorophyll Profiles in Argentine Basin'],'FontSize',20)
        xlabel('Chl [mg chloro/kg]','FontSize',18)
        ylabel('Pressure [dbars]','FontSize',18)
        
        xmin = 0;
        xmax = (size(dd.idat.Chl)).^2;
        dx = size(dd.idat.Chl);
        depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
        nlayers = size(depths,2);
        cols = hsv(nlayers);
        
        for i = 1:size(dd.idat.Chl)
            plot(dd.idat.Chl(i,:),dd.idat.P(i,:),'color',cols(i,:))
        end
        
        if titl == '12700'
            plot(dd.idat.Chl(17,:),dd.idat.P(17,:),'k.')
        end
        
        if titl == '12778'
            plot(dd.idat.Chl(6,:),dd.idat.P(6,:),'k.')
        end
      
        saveas(gcf,[titl,'_Chl_Profiles.jpg'])

        hold off
        
    end

%% 
    function heat_diff_profiles(filen,titl)
        d = load(filen);
        
        figure(5)
        
        hold on
        set(gca,'ydir','reverse')
        
        for i = 3:2:size(d.idat.heat)
            plot((d.idat.heat(i,:)-d.idat.heat(i-2,:)),d.idat.P)
        end 
        hold off
        
    end

%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
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
        ax.FontSize = 25;
        
        title([titl ': Conservative Temp Contour'])
        ylabel('Pressure [dbars]');
        xlabel('Profiles');
        contourf((d.idat.time),d.idat.P,d.idat.ct);
        
        pressure_levels = [0,4,6,8,9,10];
        
        h = colorbar;
        
        %[M,h2] = contour((d.idat.time),d.idat.P,d.idat.ct,'k','ShowText','on');
        %h2.LineWidth = 0.75;
        %clabel(M,h2,'FontSize',11,'Color','k','FontWeight','bold')
        
        [tt, zz] = meshgrid(d.idat.T,d.idat.P);
        pden = d.idat.pden - 1000;
        [C,h] = contour(d.idat.time,d.idat.P,pden,10,'k-');
        clabel(C,h,'FontSize',14,'Color','k')
        contour(d.idat.time,d.idat.P,pden,15,'k:');
        
        if titl == '12700'
            plot((d.idat.time(17,:)),d.idat.P(17,:),'k-','LineWidth',2);
        end
                
        if titl == '12700'
            xtickangle(45);
            xticks(xdates12700);             %% this uses the list of dates above
            datetick('x','dd/mm/yy','keepticks') ;
        elseif titl == '12881'
            xtickangle(45);
            xticks(xdates12881);             %% this uses the list of dates above
            datetick('x','dd/mm/yy','keepticks') ;
        end
        
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
        
        title([titl ': Heat Contour'],'FontSize',25)
        ylabel('Pressure [dbars]');
        xlabel('Profiles');
        
        colormap('spring')
        h = colorbar;
        set(get(h,'label'),'string','x 10^9','FontSize',19) %creates the label on the colorbar
        contourf((d.idat.time),d.idat.P,(d.idat.heat/1000000000));
        
        if titl == '12700'
            plot((d.idat.time(17,:)),d.idat.P(17,:),'k-','LineWidth',2);
        end
        
        pressure_levels = [0,4,6,8,9,10];
%         
%         [M,h2] = contour((d.idat.time),d.idat.P,(d.idat.heat/1000000000),'k','ShowText','on');
%         h2.LineWidth = 0.5;
%         clabel(M,h2,'FontSize',11,'Color','k','FontWeight','bold')
        
        pden = d.idat.pden - 1000;
        [C,h] = contour(d.idat.time,d.idat.P,pden,10,'k-');
        clabel(C,h,'FontSize',14,'Color','k')
        contour(d.idat.time,d.idat.P,pden,10,'k:');
                    
        if titl == '12700'
            xtickangle(45);
            xticks(xdates12700);             %% this uses the list of dates above
            datetick('x','dd/mm/yy','keepticks') ;
        elseif titl == '12881'
            xtickangle(45);
            xticks(xdates12881);             %% this uses the list of dates above
            datetick('x','dd/mm/yy','keepticks') ;
        end
        
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',13)

        hold off 
        
 %% POTENTIAL DENSITY CONTOUR       
       figure(9)
        
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
        
        title([titl ': Potential Density Contour'],'FontSize',25)
        ylabel('Pressure [dbars]');
        xlabel('Profiles');
        
        colormap('cool')
        h = colorbar;
        contourf((d.idat.time),d.idat.P,d.idat.pden);
        
        if titl == '12700'
            plot((d.idat.time(17,:)),d.idat.P(17,:),'k-','LineWidth',2)
        end
        
        pressure_levels = [0,4,6,8,9,10];
        
%         [M,h2] = contour((d.idat.time),d.idat.P,d.idat.pden,'k','ShowText','on');
%         h2.LineWidth = 0.5;
%         clabel(M,h2,'FontSize',11,'Color','k','FontWeight','bold')

        pden = d.idat.pden - 1000;
        [C,h] = contour(d.idat.time,d.idat.P,pden,10,'k-');
        clabel(C,h,'FontSize',14,'Color','k')
        contour(d.idat.time,d.idat.P,pden,10,'k:');
            
        if titl == '12700'
            xtickangle(45);
            xticks(xdates12700);             %% this uses the list of dates above
            datetick('x','dd/mm/yy','keepticks') ;
        elseif titl == '12881'
            xtickangle(45);
            xticks(xdates12881);             %% this uses the list of dates above
            datetick('x','dd/mm/yy','keepticks') ;
        end

        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',13)

        hold off 
        
%% DIFF HEAT CONTOUR       
       figure(10)
        
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
        
        title([titl ': Absolute Salinity Contour'],'FontSize',25)
        ylabel('Pressure [dbars]');
        xlabel('Profiles');
        
        colormap('parula')
        h = colorbar;
        set(get(h,'label'),'string','  ','FontSize',17) %creates the label on the colorbar
        contourf((d.idat.time),d.idat.P,d.idat.sa);
        
        if titl == '12700'
            plot((d.idat.time(17,:)),d.idat.P(17,:),'k-','LineWidth',2)
            %plot((d.idat.time(6,:)),d.idat.P(6,:),'k-','LineWidth',2)
        else 
            plot((d.idat.time(21,:)),d.idat.P(21,:),'k-','LineWidth',2)
        end
        
        pressure_levels = [0,4,6,8,9,10];
        
%         [M,h2] = contour((d.idat.time),d.idat.P,d.idat.sa,'k','ShowText','on');
%         h2.LineWidth = 0.5;
%         clabel(M,h2,'FontSize',11,'Color','k','FontWeight','bold')

        pden = d.idat.pden - 1000;
        [C,h] = contour(d.idat.time,d.idat.P,pden,10,'k-');
        clabel(C,h,'FontSize',14,'Color','k')
        contour(d.idat.time,d.idat.P,pden,10,'k:');
            
        if titl == '12700'
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

        hold off 

 
        figure(11)
        
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
        ax.FontSize = 25;
        pden = d.idat.pden - 1000;
        [C,h] = contour(d.idat.time,d.idat.P,pden,10,'k-');
        clabel(C,h,'FontSize',14,'Color','k')
        contour(d.idat.time,d.idat.P,pden,10,'k:');
        title([titl ': DIFF Heat Contour'])
        ylabel('Pressure [dbars]');
        xlabel('Profiles');
        
        for i = 2:1:size(d.idat.heat)
            heatDiff(1,:) = d.idat.heat(1,:);
            heatDiff(i,:) = d.idat.heat(i,:)-d.idat.heat(i-1,:);
        end
        
        zmin = 0;
        zmax = 1800;
        dzlayer = 200;
        depths = fliplr(-[zmin:dzlayer:zmax-dzlayer; zmin+dzlayer:dzlayer:zmax]);
        nlayers = size(depths,2);
        
        heatLayers = nan(nlayers,1);
        cols = jet(nlayers);
        fig1 = figure(11);
        set(fig1,'units','normalized','outerposition',[0 0 1 1])
        hold on
        %lgdtxt = cell(1,nlayers);
        minval = inf;
        maxval = -inf;
        
        contourf((d.idat.time),d.idat.P,heatDiff);
        
        pressure_levels = [0,4,6,8,9,10];
        
        h = colorbar;
        
%         [M,h2] = contour((d.idat.time),d.idat.P,d.idat.ct,'k','ShowText','on');
%         h2.LineWidth = 0.75;
%         clabel(M,h2,'FontSize',11,'Color','k','FontWeight','bold')
        
        pden = d.idat.pden - 1000;
        [C,h] = contour(d.idat.time,d.idat.P,pden,10,'k-');
        clabel(C,h,'FontSize',14,'Color','k')
        contour(d.idat.time,d.idat.P,pden,10,'k:');
        
        if titl == '12700'
            plot((d.idat.time(17,:)),d.idat.P(17,:),'k-','LineWidth',2);
        end
                
        if titl == '12700'
            xtickangle(45);
            xticks(xdates12700);             %% this uses the list of dates above
            datetick('x','dd/mm/yy','keepticks') ;
        elseif titl == '12881'
            xtickangle(45);
            xticks(xdates12881);             %% this uses the list of dates above
            datetick('x','dd/mm/yy','keepticks') ;
        end
        
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',13)
        
        %         dateFormat = 'dd mmm yyyy';
        %         datetick('x',dateFormat)
        
        hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%
    function contour_pden_argo(filenam,titl)
        d = load(filenam);
        
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
        
        figure(5)
        hold on
        
        set(gca,'ydir','reverse')
        set(gca, 'FontSize',18)

        ax = gca;
        ax.TickDir = 'out';
        
        title([titl ': Conservative Temp Contour'],'FontSize',20)
        ylabel('Potential Density [kg/m^3]');
        xlabel('Profiles');
        contourf((d.ipdat.time),d.ipdat.pden,d.ipdat.ct);
        
        if titl == '12700'
            plot((d.ipdat.time(17,:)),d.ipdat.pden(17,:),'k-')
        end
        
        h = colorbar;
        
        [M,h2] = contour((d.ipdat.time),d.ipdat.pden,d.ipdat.ct,'k','ShowText','on');
        h2.LineWidth = 0.5;
        clabel(M,h2,'FontSize',9,'Color','k','FontName','Ubuntu')
        
        if titl == '12700'
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
        
%         dateFormat = 'dd mmm yyyy';
%         datetick('x',dateFormat)

        hold off
        
        figure(6)
        hold on
        
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
            
        set(gca,'ydir','reverse')
        set(gca,'FontSize',18)
        
        ax = gca;
        ax.TickDir = 'out';
        
        title([titl ': Heat Contour','FontSize',20])
        ylabel('Potential Density [kg/m^3]');
        xlabel('Profiles');
        
        colormap('spring')
        h = colorbar;
        set(get(h,'label'),'string','x 10^9','FontSize',19) %creates the label on the colorbar

        contourf((d.ipdat.time),d.ipdat.pden,(d.ipdat.heat/1000000000));
        
        pressure_levels = [0,4,6,8,9,10];
        
        [M,h2] = contour((d.ipdat.time),d.ipdat.pden,(d.ipdat.heat/1000000000),'k','ShowText','on');
        h2.LineWidth = 0.5;
        clabel(M,h2,'FontSize',9,'Color','k','FontName','Ubuntu')

        if titl == '12700'
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
        
        hold off
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function  add_pden(intp,major,minor)

      % hold on
      [tt, zz] = meshgrid(d.idat.T,d.idat.P);
      pden = idat.pden - 1000;
      [C,h] = contour(datenum(tt),zz,vpden',major,'k-');
      clabel(C,h,'FontSize',14,'Color','k')
      contour(datenum(tt),zz,pden',minor,'k:');
      %lgd = legend(leg);
      %set(lgd,'fontsize',16,'location','southwest')

    end

end