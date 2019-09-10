function Heat_Eval(fn)
data = process_argo(fn)
basic_plots('12700.mat')
density_plots('12700.mat')
highlight_profile('12700.mat',33)
highlight_profile_density('12700.mat',33)

%% data processing from the .nc file to a mat file, which is interpolated
    function data = process_argo(fn)
  
        outfile = '12700.mat';
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
        t = raw_data.REFERENCE_DATE_TIME';
        ref_t = [t(1:4) '-' t(5:6) '-' t(7:8) ' ' t(9:10) ':' t(11:12) ':' t(13:14)];
        data.time = datetime(ref_t) + days(raw_data.JULD);

        data.lat  = raw_data.LATITUDE;
        data.lon  = raw_data.LONGITUDE;
        data.P  = raw_data.PRES'; %chose to use the nonadjusted vairables because they had all values of 99999
        data.T  = raw_data.TEMP';
        data.SP  = raw_data.PSAL';
        data.DO  = raw_data.DOXY';

        % missing data
        data.P(data.P >= 90000) = NaN;
        data.T(data.T >= 9000) = NaN;
        data.SP(data.SP >= 9000) = NaN;
        data.DO(data.DO >= 9000) = NaN;
        zint = 1:1:500;

        
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

%           tmp(1,:)= data.P(1,:);
%           tmpT(1,:) = data.T(1,:);
%           tmpSP(1,:) = data.SP(1,:);
%           I = find(isfinite(tmp(1,:))==1 & isfinite(tmpT(1,:))==1);
%           I
%           idat.T(1,:) = interp1(tmp(I),tmpT(I),zint)
%           
%           tmp= data.P(3,:);
%           tmpT = data.T(3,:);
%           tmpSP = data.SP(3,:);
%           I = find(isfinite(tmp)==1 & isfinite(tmpT)==1);
%           I
%           idat.T(2,:) = interp1(tmp(I),tmpT(I),zint)

       for i = 1:2:50
          tmp= data.P(i,:);
          tmpT = data.T(i,:);
          tmpSP = data.SP(i,:);
          I = find(isfinite(tmp)==1 & isfinite(tmpT)==1);
          idat.T(i,:) = interp1(tmp(I),tmpT(I),zint)
       end

%           keyboard
%               if length(I)>2
%                   idat.T(i,:) = interp1(data.P(I(i,:)),data.T(I(i,:)),zint)
%               end
%               if length(I2)>2
%                   idat.SP(i,:) = interp1(data.P(I(i,:)),data.SP(I(i,:)),zint)
%               end


        save(outfile);

    end

%% 
    function basic_plots(filen)
        
        f = load(filen);
        
        figure(1) 
        hold on 
       
        title('Temperature vs Pressure 12700')
        xlabel('In-Situ Temp (C)')
        ylabel('In-situ Pressure')

        set(gca,'fontsize',18)
        set(gca,'Ydir','reverse')

        for i = 1:size(f.data.T,1)
            plot(f.data.T(i,:),f.data.P(i,:))
        end 
        %acc_plots

        hold off
        
%         figure(2) 
%         hold on 
%        
%         title('Temperature vs Pressure 12700')
%         xlabel('In-Situ Temp (C)')
%         ylabel('In-situ Pressure')
% 
%         set(gca,'Ydir','reverse')
% 
%         for i = 1:size(f.data.T,1)
%             plot(f.data.T(i,:),f.data.P(i,:))
%         end 
% 
%         hold off
        
        
        figure(3) 
        hold on 
       
        title('Absolute Salinity vs Pressure 12700')
        xlabel('Absolute Salinity (g/kg)')
        ylabel('In-situ Pressure')

        set(gca,'fontsize',18)
        set(gca,'Ydir','reverse')

        for i = 1:size(f.data.sa,1)
            plot((f.data.sa(i,:)),f.data.P(i,:))
        end 
        %acc_plots

        hold off
        
        figure(4) 
        hold on 
       
        title('Dissolved Oxygen vs Pressure 12700')
        xlabel('Dissolved Oxygen (???)')
        ylabel('In-situ Pressure')

        set(gca,'fontsize',18)
        set(gca,'Ydir','reverse')

        for i = 1:size(f.data.DO,1)
            plot((f.data.DO(i,:)),f.data.P(i,:))
        end 
        %acc_plots

        hold off
        
    end    
  
%%
    function density_plots(fn)
        f = load(fn);
        
        figure(5) 
        hold on 
       
        title('Temperature vs Potential Density 12700')
        xlabel('In-Situ Temp (C)')
        ylabel('Potential Density')

        set(gca,'fontsize',18)
        set(gca,'Ydir','reverse')

        for i = 1:size(f.data.T,1)
            plot(f.data.T(i,:),f.data.pden(i,:))
        end 
        %acc_plots

        hold off
        
        
        figure(6) 
        hold on 
       
        title('Absolute Salinity vs Potential Density 12700')
        xlabel('Absolute Salinity (g/kg)')
        ylabel('Potential Desntiy')

        set(gca,'fontsize',18)
        set(gca,'Ydir','reverse')

        for i = 1:size(f.data.sa,1)
            plot((f.data.sa(i,:)),f.data.pden(i,:))
        end 
        %acc_plots

        hold off
        
        figure(7) 
        hold on 
       
        title('Dissolved Oxygen vs Potential Density 12700')
        xlabel('Dissolved Oxygen (???)')
        ylabel('Potential Density')

        set(gca,'fontsize',18)
        set(gca,'Ydir','reverse')

        for i = 1:size(f.data.pden,1)
            plot((f.data.DO(i,:)),f.data.pden(i,:))
        end 
        %acc_plots

        hold off
    end
     
%%
    function highlight_profile(fn,profile) 

        % This script highlights a specific profile within a whole argo .mat file
        % you input, highlighting the profile number you specify 
        %for 12700, the weird profile correlates to 33

        h = profile
        d = load(fn)

        figure(8)
        hold on
        title('Temperature Profiles')

        set(gca,'fontsize',18)
        set(gca,'ydir','reverse')
        for i = (profile-8):2:(profile+8)
            x = plot(d.data.T(i,:),d.data.P(i,:),'k.','markersize',4)
            plot(d.data.T(i,:),d.data.P(i,:),'k-','linewidth',0.5)
        end

        y = plot(d.data.T(profile,:),d.data.P(profile,:),'r.','markersize',6)
        plot(d.data.T(profile,:),d.data.P(profile,:),'r-','linewidth',0.5)

        z = plot(d.data.T(profile+2,:),d.data.P(profile+2,:),'c.','markersize',6)
        plot(d.data.T(profile+2,:),d.data.P(profile+2,:),'c-','linewidth',0.5)

        plot(d.data.T(profile-2,:),d.data.P(profile-2,:),'c.','markersize',6)
        plot(d.data.T(profile-2,:),d.data.P(profile-2,:),'c-','linewidth',0.5)

        xlabel('Temperature [C]')
        ylabel('Pressure [dbar]')
        legend([x z y],'Normal Profiles','3 Before and 3 After Profiles','Weird Profile','location','northwest')

        hold off

        figure(9)
        hold on

        title('Absolute Salinity Profiles')

        set(gca,'fontsize',18)
        set(gca,'ydir','reverse')
        for i = (profile-8):2:(profile+8)
            x= plot(d.data.sa(i,:),d.data.P(i,:),'k.','markersize',4)
            plot(d.data.sa(i,:),d.data.P(i,:),'k-','linewidth',0.5)
        end

        y = plot(d.data.sa(profile,:),d.data.P(profile,:),'r.','markersize',6)
        plot(d.data.sa(profile,:),d.data.P(profile,:),'r-','linewidth',0.5)

        z = plot(d.data.sa(profile+2,:),d.data.P(profile+2,:),'c.','markersize',6)
        plot(d.data.sa(profile+2,:),d.data.P(profile+2,:),'c-','linewidth',0.5)

        plot(d.data.sa(profile-2,:),d.data.P(profile-2,:),'c.','markersize',6)
        plot(d.data.sa(profile-2,:),d.data.P(profile-2,:),'c-','linewidth',0.5)

        ylabel('Pressure [dbar]')
        xlabel('Absolute Salinity [g/kg]')
        legend([x z y],'Normal Profiles','3 Before and 3 After Profiles','Weird Profile','location','northwest')

        %acc_plots

        hold off

        figure(10)
        hold on

        title('Dissolved Oxygen Profiles')

        set(gca,'fontsize',18)
        set(gca,'ydir','reverse')

        prof = profile+1; 

        for i = (prof-8):2:(prof+8)
            x= plot(d.data.DO(i,:),d.data.P(i,:),'k.','markersize',4)
            plot(d.data.DO(i,:),d.data.P(i,:),'k-','linewidth',0.5)
        end

        y = plot(d.data.DO(prof,:),d.data.P(prof,:),'r.','markersize',6)
        plot(d.data.DO(prof,:),d.data.P(prof,:),'r-','linewidth',0.5)

        z = plot(d.data.DO(prof+2,:),d.data.P(prof+2,:),'c.','markersize',6)
        plot(d.data.DO(prof+2,:),d.data.P(prof+2,:),'c-','linewidth',0.5)

        plot(d.data.DO(prof-2,:),d.data.P(prof-2,:),'c.','markersize',6)
        plot(d.data.DO(prof-2,:),d.data.P(prof-2,:),'c-','linewidth',0.5)

        ylabel('Pressure [dbar]')
        xlabel('Dissolved Oxygen [mol/m^2]')
        legend([x z y],'Normal Profiles','3 Before and 3 After Profiles','Weird Profile','location','northwest')

        %acc_plots

        hold off

    end


%% 
    function highlight_profile_density(fn,profile) 

        % This script highlights a specific profile within a whole argo .mat file
        % you input, highlighting the profile number you specify 
        %for 12700, the weird profile correlates to 33

        h = profile
        d = load(fn)

        figure(11)
        hold on
        title('Temperature Profiles')

        set(gca,'fontsize',18)
        set(gca,'ydir','reverse')
        for i = (profile-8):2:(profile+8)
            x = plot(d.data.T(i,:),d.data.pden(i,:),'k.','markersize',4)
            plot(d.data.T(i,:),d.data.pden(i,:),'k-','linewidth',0.5)
        end

        y = plot(d.data.T(profile,:),d.data.pden(profile,:),'r.','markersize',6)
        plot(d.data.T(profile,:),d.data.pden(profile,:),'r-','linewidth',0.5)

        z = plot(d.data.T(profile+2,:),d.data.pden(profile+2,:),'c.','markersize',6)
        plot(d.data.T(profile+2,:),d.data.pden(profile+2,:),'c-','linewidth',0.5)

        plot(d.data.T(profile-2,:),d.data.pden(profile-2,:),'c.','markersize',6)
        plot(d.data.T(profile-2,:),d.data.pden(profile-2,:),'c-','linewidth',0.5)

        xlabel('Temperature [C]')
        ylabel('Pressure [dbar]')
        legend([x z y],'Normal Profiles','3 Before and 3 After Profiles','Weird Profile','location','northwest')

        hold off

        figure(12)
        hold on

        title('Absolute Salinity Profiles')

        set(gca,'fontsize',18)
        set(gca,'ydir','reverse')
        for i = (profile-8):2:(profile+8)
            x= plot(d.data.sa(i,:),d.data.pden(i,:),'k.','markersize',4)
            plot(d.data.sa(i,:),d.data.pden(i,:),'k-','linewidth',0.5)
        end

        y = plot(d.data.sa(profile,:),d.data.pden(profile,:),'r.','markersize',6)
        plot(d.data.sa(profile,:),d.data.pden(profile,:),'r-','linewidth',0.5)

        z = plot(d.data.sa(profile+2,:),d.data.pden(profile+2,:),'c.','markersize',6)
        plot(d.data.sa(profile+2,:),d.data.pden(profile+2,:),'c-','linewidth',0.5)

        plot(d.data.sa(profile-2,:),d.data.pden(profile-2,:),'c.','markersize',6)
        plot(d.data.sa(profile-2,:),d.data.pden(profile-2,:),'c-','linewidth',0.5)

        ylabel('Pressure [dbar]')
        xlabel('Absolute Salinity [g/kg]')
        legend([x z y],'Normal Profiles','3 Before and 3 After Profiles','Weird Profile','location','northwest')

        %acc_plots

        hold off

        figure(13)
        hold on

        title('Dissolved Oxygen Profiles')

        set(gca,'fontsize',18)
        set(gca,'ydir','reverse')

        prof = profile +1; 

        for i = (prof-8):2:(prof+8)
            x= plot(d.data.DO(i,:),d.data.pden(i+1,:),'k.','markersize',4)
            plot(d.data.DO(i,:),d.data.pden(i+1,:),'k-','linewidth',0.5)
        end

        y = plot(d.data.DO(prof,:),d.data.pden(profile,:),'r.','markersize',6)
        plot(d.data.DO(prof,:),d.data.pden(profile,:),'r-','linewidth',0.5)

        z = plot(d.data.DO(prof+2,:),d.data.pden(profile+2,:),'c.','markersize',6)
        plot(d.data.DO(prof+2,:),d.data.pden(profile+2,:),'c-','linewidth',0.5)

        plot(d.data.DO(prof-2,:),d.data.pden(profile-2,:),'c.','markersize',6)
        plot(d.data.DO(prof-2,:),d.data.pden(profile-2,:),'c-','linewidth',0.5)

        ylabel('Pressure [dbar]')
        xlabel('Dissolved Oxygen [mol/m^2]')
        legend([x z y],'Normal Profiles','3 Before and 3 After Profiles','Weird Profile','location','northwest')

        %acc_plots

        hold off

    end
end