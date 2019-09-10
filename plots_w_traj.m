function plots_w_traj(cruise,titl)

DO_w_traj(cruise,titl)
%T_w_traj(cruise,titl)
%SA_w_traj(cruise,titl)

    function DO_w_traj(cruise,titl)    
        dd = load([titl,'.mat']);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% DO plot %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        close all
        
        figure(1)
        
        pos2 = [0.48 0.1 0.49 0.85];
        %left bottom width height
        subplot('Position',pos2)
        
        hold on
        grid on
        set(gca,'ydir','reverse')
        set(gca,'FontSize',25)
        
        title([titl,': Dissolved Oxygen Profiles in Argentine Basin'],'FontSize',25)
        xlabel('Diss Oxy [micro-mol/kg]','FontSize',25)
        ylabel('Pressure [dbars]','FontSize',25)
        
        xmin = 0;
        xmax = (size(dd.data.lat)).^2;
        dx = size(dd.data.lat);
        depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
        nlayers = size(depths,2);
        cols = hsv(nlayers);
        
        for i = 1:size(dd.idat.DO)
            plot(dd.idat.DO(i,:),dd.idat.P(i,:),'color',cols(i,:))
        end
        
        if titl == '12700'
            mar4th = plot(dd.idat.DO(17,:),dd.idat.P(17,:),'r-','LineWidth',2);
        elseif titl == '12778'
            dec20th = plot(dd.idat.DO(6,:),dd.idat.P(6,:),'r-','LineWidth',2);
        end
%         dec20th = plot(dd.idat.DO(106,:),dd.idat.P(106,:),'r-','LineWidth',2)
%         
%         mar2nd = plot(dd.idat.DO(113,:),dd.idat.P(113,:),'k-','LineWidth',1.5)
%         legend([dec20th mar2nd],'Profile taken on Dec 20th','Profile March 2nd');
        
        hold off
        
        pos1 = [0.05 0.1 0.35 0.85];
        subplot('Position',pos1)
        
         %   %closes all open figure windows 
        %=========================================
        Meta = MetaFile(cruise);

        hold on   %keep figure up throughout the script so that we can keep plotting awesome data over it!
        title([titl,': Trajectory from ',datestr(dd.data.time(1,:)),' to ',datestr(dd.data.time(end,:))],'FontSize',21)

        Bath = load(Meta.BathyFile);         %loads Bathymetry file for cruise specified above
        Bathf = fields(Bath);           %gets the fields from the loaded Bathymetry file
        BathFile = Bath.([Bathf{:}]);   %creates final bathfile

        [numrows, numcols] = size(BathFile.z);  %finds the size of the depth part of the Bathymetry file and saves the number of rows/collumns
        for r = 1:numrows        %double for loop goes through all of the data in the bath file
            for c = 1:numcols
                if (BathFile.z(r,c)) >= 0       %if the elevation is greater or equal to zero
                    BathFile.z(r,c) = 51;      %set the elevation to 51 (the colorscale only colors it in as land if it is greater than 50 so this makes all the land colored in)
                end
            end
        end


        [ahomap,cblev] = AHO_gebco_matlab;      %makes the colormap the gebco map and assigns rgb color code values to ahomap and assigns depths to cblev
        caxis([cblev(1) cblev(end)])            %makes the axis go from 100 (tan land color) to -6000 (dark blue)
        cba = colorbar;                         %cba is our colorbar variable
        % se added comments cba.Label.String = 'Depth (meters)';
        set(cba,'YTick', cblev(1:end-1))        %set tick marks on cba to show depths from cblev (all except the last the land one (100): 0, -200, -500, -1000, -2000, -3000, -4000, -5000, -6000), these ticks correspond to changing color levels
        colormap(ahomap)                        %makes the colormap use the ahomap rgb codes


        [X,Y] = meshgrid(BathFile.x,BathFile.y);      %makes grid using longitudes and latitudes from the Bathymetry file
        imagesc(BathFile.x,BathFile.y,BathFile.z');   %plots image using corresponding lon/lat and the depths from the Bathymetry file (rotated so that they match up correctly)
        set(gca,'YDir','normal');        %gets handle of current axis (gca) and sets it to a normal y-direction
        set(gca,'YLim',[Meta.LatMin Meta.LatMax]); % sets map y-limits to the Latitudes specified for your area
        set(gca,'XLim',[Meta.LonMin Meta.LonMax]); % sets map x-limits to the Longitudes specified for you area
        set(gca,'FontSize',21)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %==================================================================
        d = load([titl,'.mat']);

        Latitude = dd.data.lat;   %get Latitude variable
        Longitude = dd.data.lon;
        Cruise = line(Longitude(2:end),Latitude(2:end),'LineWidth',2,'Color',[0 0 0]);
        
        xmin = 0;
        xmax = (size(dd.data.lat)).^2;
        dx = size(dd.data.lat);
        depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
        nlayers = size(depths,2);
        cols = jet(nlayers);
        
        for row = 1:size(dd.data.lat)    %loop through the positions in the file
             %get Longitude variable
            Flo = plot(Longitude(row),Latitude(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',8,'MarkerEdgeColor',cols(row,:)); 
            if titl == '12700'
                mar4th = plot(Longitude(17,:),Latitude(17,:),'.','MarkerSize',Meta.MarkerSize,'LineWidth',6,'MarkerEdgeColor',[0 0 0]);
            elseif titl == '12778'
                dec20th = plot(Longitude(6,:),Latitude(6,:),'.','MarkerSize',Meta.MarkerSize,'LineWidth',6,'MarkerEdgeColor',[0 0 0]);;
            end 
        end
        
        %For 12700:
%         for row = 17 %length(d.data.lat)    %loop through the positions in the file
%              %get Longitude variable
%             mar2nd = plot(Longitude(row),Latitude(row),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
%             %WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);
%         end
        
        
%         for row = 113 %length(d.data.lat)    %loop through the positions in the file
%              %get Longitude variable
%             mar2nd = plot(Longitude(row),Latitude(row),'.','MarkerSize',2.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]); 
%             %WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]); 
%         end
%         
%         for row = 106 %length(d.data.lat)    %loop through the positions in the file
%              %get Longitude variable
%             dec20th = plot(Longitude(row),Latitude(row),'.','MarkerSize',2.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]); 
%             %WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]); 
%         end


        %legend([WeirdFlo WeirdFlo2],'Profile 17 of 12700','Profile 6 of 12788')

        %==================================================================
        
        ncfname   = fullfile(Meta.common_loc,Meta.SeaIceFile1); %gets ice file by calling fullfile(FolderName,'FileName') and referring to variables specified at the top
        ncid = netcdf.open(ncfname,'NC_NOWRITE');     %opens ice file so that you can't edit the file and assigns file to "ncid"
        [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);   %keeps track of netcdf variables
        for i = 1:numvars                             %loops through variables
            [varnames{i}, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,i-1);   %saves variables
            Data.([varnames{i}]) = netcdf.getVar(ncid,i-1);   %finds what variables are and saves all of them into Data so that you can access them
        end
        netcdf.close(ncid);  %closes netcdf file
        [numrows,numcols] = size(Data.longitude);  %finds size of longitude,latitude,and Ice Conc arrays (they will each be same size) and assigns those values to numrows and numcols
        for r = 1:numrows  %double for loop goes through each row and collumn or Ice Concentration array
            for c = 1:numcols  %goes through all the collumns
                if (isnan(Data.IceConc(r,c)) == 1) || (Data.IceConc(r,c) == 0) %checks for NaN data and puts 0 in instead
                    Data.IceConc(r,c)=0;
                elseif ((Data.longitude(r,c) > 180) && (Data.longitude(r,c) < 360)) == 1   %if longitude values are between 180 and 360 in the ice file, subtract 360 from them to make the value match up to our -180 to 180 scale
                    Ice = plot(Data.longitude(r,c)-360, Data.latitude(r,c),'o','MarkerSize',4.5,'LineWidth',1,'MarkerEdgeColor',[0 1 1]); %plot where the ice is at the correct lat lon with cyan circles (keep middle of cirles clear so that you can see bathymetry levels) and assign to "ice" so that we can use in the legend
                else   %if the ice location is between 0 and 180 degrees east
                    Ice = plot(Data.longitude(r,c), Data.latitude(r,c),'o','MarkerSize',4.5,'LineWidth',1,'MarkerEdgeColor',[0 1 1]); %plot where the ice is at the correct lat lon with cyan circles (keep middle of circles clear so that you can see bathymetry) and assign to "ice" so that we can use in legend
                end
            end
        end

        ncfname   = fullfile(Meta.common_loc,Meta.SeaIceFile2); %gets second ice file by calling fullfile(FolderName,'FileName') and referring to variables specified above
        ncid = netcdf.open(ncfname,'NC_NOWRITE');     %opens ice file so that you can't edit the file and assigns file to "ncid"
        [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);  %keeps track of netcdf variables
        for i = 1:numvars                                                  %loops through variables
            [varnames{i}, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,i-1);    %saves variables
            Data.([varnames{i}]) = netcdf.getVar(ncid,i-1);   %finds what variables are and saves all of them into Data so that you can access them
        end
        netcdf.close(ncid);
        [numrows,numcols] = size(Data.longitude);  %finds size of longitude,latitude,and Ice Conc arrays (they will each be same size) and assigns those values to numrows and numcols
        for r = 1:numrows  %double for loop goes through each row and each collumn
            for c = 1:numcols
                if (isnan(Data.IceConc(r,c)) == 1) || (Data.IceConc(r,c) == 0) %checks for NaN data and puts 0 in instead
                    Data.IceConc(r,c)=0;
                elseif ((Data.longitude(r,c) > 180) && (Data.longitude(r,c) < 360)) == 1    %if longitude values are between 180 and 360 in the ice file, subtract 360 from them to make the value match up to our -180 to 180 scale
                    Feb = plot(Data.longitude(r,c)-360, Data.latitude(r,c),'^','MarkerSize',5,'LineWidth',.45,'MarkerEdgeColor',[0 0 1]); %plot where the ice is at the correct lat lon with a blue triangle (middle open so that you can see bathymetry) and assign to "Feb" so that we can use in legend
                else              %if longitude values for ice are between 0 and 180 east
                    Feb = plot(Data.longitude(r,c), Data.latitude(r,c),'^','MarkerSize',5,'LineWidth',.45,'MarkerEdgeColor',[0 0 1]);    %plot where the ice is at the correct lat lon with blue triangle (middle open so that you can see bathymetry) and assign to "Feb so that we can use in legend
                end
            end
        end



        SAF = load(fullfile(Meta.common_loc,'saf.asc.txt'));                            %load Subantarctic Front file
        plot(SAF(:,1),SAF(:,2),'-','Color',[.9 .87 .88],'LineWidth', 1.5)      %plot the longitude/latitude of this front in a tan color
        % This code does not work:  see working code below commented lines
        %for row = 1:80:length(SAF)                            %loop through length of SAF file so that you label every nth point (adjust n to make make more/less cluttered)
        %    x = SAF(row,1);                                   %assign x to longitude of SAF
        %    y = SAF(row,2);                                   %ssign y to latitude of SAF
        %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))  %only add text if the text will stay inside the range of the map
        %        text(x,y,'SAF','FontSize',8)                  %labels the fronts
        %    end
        %end
        % use below instead
        %    x = SAF(row,1);                                   %assign x to longitude of SAF
        %    y = SAF(row,2);                                   %ssign y to latitude of SAF
        ipos = find( (SAF(:,1) > (Meta.LonMin+1)) & (SAF(:,1) < (Meta.LonMax-1)) & (SAF(:,2) > (Meta.LatMin+1)) & (SAF(:,2) < (Meta.LatMax-1)));  %only add text if the text will stay inside the range of the map
        text(SAF(ipos(1:40:end),1),SAF(ipos(1:40:end),2),'SAF','FontSize',8)  %label the front only once
        clear ipos

        SBDY = load(fullfile(Meta.common_loc,'sbdy.asc.txt'));                           %load the southern boundary file
        plot(SBDY(:,1),SBDY(:,2),'-','Color',[.9 .87 .88],'LineWidth', 1.5)     %plot the longitude/latitude of this in a tan color
        %for row = 1:18:length(SBDY)                            %loop through length of SBDY file so that you label every nth point (adjust n to make make more/less cluttered)
        %    x = SBDY(row,1);                                   %assign x to longitude of SBDY
        %    y = SBDY(row,2);                                   %assign y to latitude of SBDY
        %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))  %only add text if the text will stay inside the range of the map
        %        text(x,y,'SBDY','FontSize',8)                  %add text
        %    end
        %end
        ipos = find( (SBDY(:,1) > (Meta.LonMin+1)) & (SBDY(:,1) < (Meta.LonMax-1)) & (SBDY(:,2) > (Meta.LatMin+1)) & (SBDY(:,2) < (Meta.LatMax-1))) ; %only add text if the text will stay inside the range of the map
        text(SBDY(ipos(1:40:end),1),SBDY(ipos(1:40:end),2),'SBDY','FontSize',8)  %label the front only once
        clear ipos

        STF = load(fullfile(Meta.common_loc,'stf.asc.txt'));                              %loads Subtropical Front file (only appears on a few maps because it is farther north)
        plot(STF(:,1),STF(:,2),'-','Color', [.9 .87 .88],'LineWidth', 1.5)       %plot the longitude/latitude of this in a tan color
        %for row = 1:30:length(STF)                              %loop through length of STF file so that you label every nth point (adjust n to make make more/less cluttered)
        %    x = STF(row,1);                                     %assign x to longitude of STF
        %    y = STF(row,2);                                     %assign y to latitude of STF
        %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))   %only add text if the text will stay inside the range of the map
        %        text(x,y,'STF','FontSize',8)                    %add text
        %    end
        %end
        ipos = find( (STF(:,1) > (Meta.LonMin+1)) & (STF(:,1) < (Meta.LonMax-1)) & (STF(:,2) > (Meta.LatMin+1)) & (STF(:,2) < (Meta.LatMax-1))) ; %only add text if the text will stay inside the range of the map
        if ipos    
           text(STF(ipos(1:40:end),1),STF(ipos(1:40:end),2),'STF','FontSize',8)  %label the front only once
        end
        clear ipos

        PF = load(fullfile(Meta.common_loc,'pf.asc.txt'));                               %loads Polar Front file
        plot(PF(:,1),PF(:,2),'-','Color',[.9 .87 .88],'LineWidth', 1.5)         %plot the longitude/latitude of this in a tan color
        %for row = 1:90:length(PF)                              %loop through length of PF file so that you label every nth point (adjust n to make make more/less cluttered)
        %    x = PF(row,1);                                     %assign x to longitude of STF file
        %    y = PF(row,2);                                     %assign y to latitude of STF file
        %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))  %only add text if the text will stay inside the range of the map
        %        text(x,y,'PF','FontSize',8)                    %add text
        %    end
        %end
        ipos = find( (PF(:,1) > (Meta.LonMin+1)) & (PF(:,1) < (Meta.LonMax-1)) & (PF(:,2) > (Meta.LatMin+1)) & (PF(:,2) < (Meta.LatMax-1)));  %only add text if the text will stay inside the range of the map
        text(PF(ipos(1:40:end),1),PF(ipos(1:40:end),2),'PF','FontSize',8)  %label the front 
        clear ipos

        SACCF = load(fullfile(Meta.common_loc,'saccf.asc.txt'));                        %load the Southern Atlantic Circumpolar Current Front file
        plot(SACCF(:,1),SACCF(:,2),'-','Color',[.9 .87 .88],'LineWidth', 1.5)  %plot the longitude/latitude of this front in a tan color
        %for row = 1:25:length(SACCF)                          %loop through length of SACCF file so that you label every nth point (adjust n to make make more/less cluttered)
        %    x = SACCF(row,1);                                 %assign x to longitude of SACCF
        %    y = SACCF(row,2);                                 %assign y to latitude of SACCF
        %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))     %only add text if the text will stay inside the range of the map
        %        text(x,y,'SACCF','FontSize', 8)               %labels the fronts
        %    end
        %end
        ipos = find( (SACCF(:,1) > (Meta.LonMin+1)) & (SACCF(:,1) < (Meta.LonMax-1)) & (SACCF(:,2) > (Meta.LatMin+1)) & (SACCF(:,2) < (Meta.LatMax-1))) ; %only add text if the text will stay inside the range of the map
        text(SACCF(ipos(1:40:end),1),SACCF(ipos(1:40:end),2),'SACCF','FontSize',8)  %label the front only once
        clear ipos

        contour(BathFile.x,BathFile.y,BathFile.z',[0 0],'Color','k')                  % outline the coasts with black contour by accessing Bathymetry file lon/lat/dep
        [L,m] = contour(BathFile.x,BathFile.y,BathFile.z',[-4000 -3000 -2000 -1000],'Color','k'); % -2000m isobath (farthest floats go to) black countour and assign lon/lat positions to L/M
        clabel(L,m,'LabelSpacing',400);                                               %label contour with "2000" every nth datapoint

        ylabel('Latitude(degrees)','FontSize',16)    %add y label for Latitude
        xlabel('Longitude(degrees)','FontSize',16)   %add x label for Longitude


        % if  strcmp(area,'OOI')      %if your area is OOI
        %     legend([Ice Feb More Flo], 'Ice on September 20, 2006', 'Ice on February 20, 2006','Station Locations','Mooring','Proposed Floats','Location','Northwest')  %add legend
        % elseif strcmp(area,'CSIRO') %if your area is CSIRO
        %     legend([Ice Feb Flo], 'Ice on September 20, 2006', 'Ice on February 20, 2006','Proposed Floats','Location','Northwest')  %add legend

        saveas(gcf,[titl,'_DO_w_Traj.jpg'])
        
        hold off    %we are done plotting
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% T plot %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%     function T_w_traj(cruise,titl)    
%         figure(2)
%         
%         dd = load([titl,'.mat']);
%         
%         pos2 = [0.48 0.1 0.49 0.85];
%         %left bottom width height
%         subplot('Position',pos2)
%         
%         hold on
%         grid on
%         set(gca,'ydir','reverse')
%         set(gca,'FontSize',22)
%         
%         title([titl,' :Temperature Profiles in Argentine Basin'],'FontSize',22)
%         xlabel('Temperature [C]','FontSize',22)
%         ylabel('Pressure [dbars]','FontSize',22)
%         
%         xmin = 0;
%         xmax = (size(dd.data.lat)).^2;
%         dx = size(dd.data.lat);
%         depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
%         nlayers = size(depths,2);
%         cols = hsv(nlayers);
%         
%         for i = 1:size(dd.idat.T)
%             plot(dd.idat.T(i,:),dd.idat.P(i,:),'color',cols(i,:))
%         end
%         
%         if titl == '12700'
%             mar4th = plot(dd.idat.T(17,:),dd.idat.P(17,:),'k-','LineWidth',2);
%         elseif titl == '12778'
%             dec20th = plot(dd.idat.T(6,:),dd.idat.P(6,:),'k-','LineWidth',2);
%         end
% %         dec20th = plot(dd.idat.DO(106,:),dd.idat.P(106,:),'r-','LineWidth',2)
% %         
% %         mar2nd = plot(dd.idat.DO(113,:),dd.idat.P(113,:),'k-','LineWidth',1.5)
% %         legend([dec20th mar2nd],'Profile taken on Dec 20th','Profile March 2nd');
%         
%         hold off
%         
%         pos1 = [0.05 0.1 0.35 0.85];
%         subplot('Position',pos1)
%         
%          %   %closes all open figure windows 
%         %=========================================
%         Meta = MetaFile(cruise);
% 
%         hold on   %keep figure up throughout the script so that we can keep plotting awesome data over it!
%         title([titl,': Trajectory from ',datestr(dd.data.time(1,:)),' to ',datestr(dd.data.time(end,:))],'FontSize',15)
% 
%         Bath = load(Meta.BathyFile);         %loads Bathymetry file for cruise specified above
%         Bathf = fields(Bath);           %gets the fields from the loaded Bathymetry file
%         BathFile = Bath.([Bathf{:}]);   %creates final bathfile
% 
%         [numrows, numcols] = size(BathFile.z);  %finds the size of the depth part of the Bathymetry file and saves the number of rows/collumns
%         for r = 1:numrows        %double for loop goes through all of the data in the bath file
%             for c = 1:numcols
%                 if (BathFile.z(r,c)) >= 0       %if the elevation is greater or equal to zero
%                     BathFile.z(r,c) = 51;      %set the elevation to 51 (the colorscale only colors it in as land if it is greater than 50 so this makes all the land colored in)
%                 end
%             end
%         end
% 
% 
%         [ahomap,cblev] = AHO_gebco_matlab;      %makes the colormap the gebco map and assigns rgb color code values to ahomap and assigns depths to cblev
%         caxis([cblev(1) cblev(end)])            %makes the axis go from 100 (tan land color) to -6000 (dark blue)
%         cba = colorbar;                         %cba is our colorbar variable
%         % se added comments cba.Label.String = 'Depth (meters)';
%         set(cba,'YTick', cblev(1:end-1))        %set tick marks on cba to show depths from cblev (all except the last the land one (100): 0, -200, -500, -1000, -2000, -3000, -4000, -5000, -6000), these ticks correspond to changing color levels
%         colormap(ahomap)                        %makes the colormap use the ahomap rgb codes
% 
% 
%         [X,Y] = meshgrid(BathFile.x,BathFile.y);      %makes grid using longitudes and latitudes from the Bathymetry file
%         imagesc(BathFile.x,BathFile.y,BathFile.z');   %plots image using corresponding lon/lat and the depths from the Bathymetry file (rotated so that they match up correctly)
%         set(gca,'YDir','normal');        %gets handle of current axis (gca) and sets it to a normal y-direction
%         set(gca,'YLim',[Meta.LatMin Meta.LatMax]); % sets map y-limits to the Latitudes specified for your area
%         set(gca,'XLim',[Meta.LonMin Meta.LonMax]); % sets map x-limits to the Longitudes specified for you area
%         set(gca,'FontSize',18)
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %==================================================================
%         d = load([titl,'.mat']);
% 
%         Latitude = dd.data.lat;   %get Latitude variable
%         Longitude = dd.data.lon;
%         Cruise = line(Longitude(1:end),Latitude(1:end),'LineWidth',2,'Color',[0 0 0]);
%         
%         xmin = 0;
%         xmax = (size(dd.data.lat)).^2;
%         dx = size(dd.data.lat);
%         depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
%         nlayers = size(depths,2);
%         cols = hsv(nlayers);
%         
%         for row = 1:size(dd.data.lat)    %loop through the positions in the file
%              %get Longitude variable
%             Flo = plot(Longitude(row),Latitude(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(row,:)); 
%             if titl == '12700'
%                 mar4th = plot(Longitude(17,:),Latitude(17,:),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
%             elseif titl == '12778'
%                 dec20th = plot(Longitude(6,:),Latitude(6,:),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);;
%             end  
%         end
%         
%         %For 12700:
% %         for row = 17 %length(d.data.lat)    %loop through the positions in the file
% %              %get Longitude variable
% %             mar2nd = plot(Longitude(row),Latitude(row),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
% %             %WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);
% %         end
%         
%         
% %         for row = 113 %length(d.data.lat)    %loop through the positions in the file
% %              %get Longitude variable
% %             mar2nd = plot(Longitude(row),Latitude(row),'.','MarkerSize',2.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]); 
% %             %WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]); 
% %         end
% %         
% %         for row = 106 %length(d.data.lat)    %loop through the positions in the file
% %              %get Longitude variable
% %             dec20th = plot(Longitude(row),Latitude(row),'.','MarkerSize',2.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]); 
% %             %WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]); 
% %         end
% 
% 
%         %legend([WeirdFlo WeirdFlo2],'Profile 17 of 12700','Profile 6 of 12788')
% 
%         %==================================================================
%         
%         ncfname   = fullfile(Meta.common_loc,Meta.SeaIceFile1); %gets ice file by calling fullfile(FolderName,'FileName') and referring to variables specified at the top
%         ncid = netcdf.open(ncfname,'NC_NOWRITE');     %opens ice file so that you can't edit the file and assigns file to "ncid"
%         [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);   %keeps track of netcdf variables
%         for i = 1:numvars                             %loops through variables
%             [varnames{i}, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,i-1);   %saves variables
%             Data.([varnames{i}]) = netcdf.getVar(ncid,i-1);   %finds what variables are and saves all of them into Data so that you can access them
%         end
%         netcdf.close(ncid);  %closes netcdf file
%         [numrows,numcols] = size(Data.longitude);  %finds size of longitude,latitude,and Ice Conc arrays (they will each be same size) and assigns those values to numrows and numcols
%         for r = 1:numrows  %double for loop goes through each row and collumn or Ice Concentration array
%             for c = 1:numcols  %goes through all the collumns
%                 if (isnan(Data.IceConc(r,c)) == 1) || (Data.IceConc(r,c) == 0) %checks for NaN data and puts 0 in instead
%                     Data.IceConc(r,c)=0;
%                 elseif ((Data.longitude(r,c) > 180) && (Data.longitude(r,c) < 360)) == 1   %if longitude values are between 180 and 360 in the ice file, subtract 360 from them to make the value match up to our -180 to 180 scale
%                     Ice = plot(Data.longitude(r,c)-360, Data.latitude(r,c),'o','MarkerSize',4.5,'LineWidth',1,'MarkerEdgeColor',[0 1 1]); %plot where the ice is at the correct lat lon with cyan circles (keep middle of cirles clear so that you can see bathymetry levels) and assign to "ice" so that we can use in the legend
%                 else   %if the ice location is between 0 and 180 degrees east
%                     Ice = plot(Data.longitude(r,c), Data.latitude(r,c),'o','MarkerSize',4.5,'LineWidth',1,'MarkerEdgeColor',[0 1 1]); %plot where the ice is at the correct lat lon with cyan circles (keep middle of circles clear so that you can see bathymetry) and assign to "ice" so that we can use in legend
%                 end
%             end
%         end
% 
%         ncfname   = fullfile(Meta.common_loc,Meta.SeaIceFile2); %gets second ice file by calling fullfile(FolderName,'FileName') and referring to variables specified above
%         ncid = netcdf.open(ncfname,'NC_NOWRITE');     %opens ice file so that you can't edit the file and assigns file to "ncid"
%         [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);  %keeps track of netcdf variables
%         for i = 1:numvars                                                  %loops through variables
%             [varnames{i}, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,i-1);    %saves variables
%             Data.([varnames{i}]) = netcdf.getVar(ncid,i-1);   %finds what variables are and saves all of them into Data so that you can access them
%         end
%         netcdf.close(ncid);
%         [numrows,numcols] = size(Data.longitude);  %finds size of longitude,latitude,and Ice Conc arrays (they will each be same size) and assigns those values to numrows and numcols
%         for r = 1:numrows  %double for loop goes through each row and each collumn
%             for c = 1:numcols
%                 if (isnan(Data.IceConc(r,c)) == 1) || (Data.IceConc(r,c) == 0) %checks for NaN data and puts 0 in instead
%                     Data.IceConc(r,c)=0;
%                 elseif ((Data.longitude(r,c) > 180) && (Data.longitude(r,c) < 360)) == 1    %if longitude values are between 180 and 360 in the ice file, subtract 360 from them to make the value match up to our -180 to 180 scale
%                     Feb = plot(Data.longitude(r,c)-360, Data.latitude(r,c),'^','MarkerSize',5,'LineWidth',.45,'MarkerEdgeColor',[0 0 1]); %plot where the ice is at the correct lat lon with a blue triangle (middle open so that you can see bathymetry) and assign to "Feb" so that we can use in legend
%                 else              %if longitude values for ice are between 0 and 180 east
%                     Feb = plot(Data.longitude(r,c), Data.latitude(r,c),'^','MarkerSize',5,'LineWidth',.45,'MarkerEdgeColor',[0 0 1]);    %plot where the ice is at the correct lat lon with blue triangle (middle open so that you can see bathymetry) and assign to "Feb so that we can use in legend
%                 end
%             end
%         end
% 
% 
% 
%         SAF = load(fullfile(Meta.common_loc,'saf.asc.txt'));                            %load Subantarctic Front file
%         plot(SAF(:,1),SAF(:,2),'-','Color',[.9 .87 .88],'LineWidth', 1.5)      %plot the longitude/latitude of this front in a tan color
%         % This code does not work:  see working code below commented lines
%         %for row = 1:80:length(SAF)                            %loop through length of SAF file so that you label every nth point (adjust n to make make more/less cluttered)
%         %    x = SAF(row,1);                                   %assign x to longitude of SAF
%         %    y = SAF(row,2);                                   %ssign y to latitude of SAF
%         %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))  %only add text if the text will stay inside the range of the map
%         %        text(x,y,'SAF','FontSize',8)                  %labels the fronts
%         %    end
%         %end
%         % use below instead
%         %    x = SAF(row,1);                                   %assign x to longitude of SAF
%         %    y = SAF(row,2);                                   %ssign y to latitude of SAF
%         ipos = find( (SAF(:,1) > (Meta.LonMin+1)) & (SAF(:,1) < (Meta.LonMax-1)) & (SAF(:,2) > (Meta.LatMin+1)) & (SAF(:,2) < (Meta.LatMax-1)));  %only add text if the text will stay inside the range of the map
%         text(SAF(ipos(1:40:end),1),SAF(ipos(1:40:end),2),'SAF','FontSize',8)  %label the front only once
%         clear ipos
% 
%         SBDY = load(fullfile(Meta.common_loc,'sbdy.asc.txt'));                           %load the southern boundary file
%         plot(SBDY(:,1),SBDY(:,2),'-','Color',[.9 .87 .88],'LineWidth', 1.5)     %plot the longitude/latitude of this in a tan color
%         %for row = 1:18:length(SBDY)                            %loop through length of SBDY file so that you label every nth point (adjust n to make make more/less cluttered)
%         %    x = SBDY(row,1);                                   %assign x to longitude of SBDY
%         %    y = SBDY(row,2);                                   %assign y to latitude of SBDY
%         %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))  %only add text if the text will stay inside the range of the map
%         %        text(x,y,'SBDY','FontSize',8)                  %add text
%         %    end
%         %end
%         ipos = find( (SBDY(:,1) > (Meta.LonMin+1)) & (SBDY(:,1) < (Meta.LonMax-1)) & (SBDY(:,2) > (Meta.LatMin+1)) & (SBDY(:,2) < (Meta.LatMax-1))) ; %only add text if the text will stay inside the range of the map
%         text(SBDY(ipos(1:40:end),1),SBDY(ipos(1:40:end),2),'SBDY','FontSize',8)  %label the front only once
%         clear ipos
% 
%         STF = load(fullfile(Meta.common_loc,'stf.asc.txt'));                              %loads Subtropical Front file (only appears on a few maps because it is farther north)
%         plot(STF(:,1),STF(:,2),'-','Color', [.9 .87 .88],'LineWidth', 1.5)       %plot the longitude/latitude of this in a tan color
%         %for row = 1:30:length(STF)                              %loop through length of STF file so that you label every nth point (adjust n to make make more/less cluttered)
%         %    x = STF(row,1);                                     %assign x to longitude of STF
%         %    y = STF(row,2);                                     %assign y to latitude of STF
%         %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))   %only add text if the text will stay inside the range of the map
%         %        text(x,y,'STF','FontSize',8)                    %add text
%         %    end
%         %end
%         ipos = find( (STF(:,1) > (Meta.LonMin+1)) & (STF(:,1) < (Meta.LonMax-1)) & (STF(:,2) > (Meta.LatMin+1)) & (STF(:,2) < (Meta.LatMax-1))) ; %only add text if the text will stay inside the range of the map
%         if ipos    
%            text(STF(ipos(1:40:end),1),STF(ipos(1:40:end),2),'STF','FontSize',8)  %label the front only once
%         end
%         clear ipos
% 
%         PF = load(fullfile(Meta.common_loc,'pf.asc.txt'));                               %loads Polar Front file
%         plot(PF(:,1),PF(:,2),'-','Color',[.9 .87 .88],'LineWidth', 1.5)         %plot the longitude/latitude of this in a tan color
%         %for row = 1:90:length(PF)                              %loop through length of PF file so that you label every nth point (adjust n to make make more/less cluttered)
%         %    x = PF(row,1);                                     %assign x to longitude of STF file
%         %    y = PF(row,2);                                     %assign y to latitude of STF file
%         %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))  %only add text if the text will stay inside the range of the map
%         %        text(x,y,'PF','FontSize',8)                    %add text
%         %    end
%         %end
%         ipos = find( (PF(:,1) > (Meta.LonMin+1)) & (PF(:,1) < (Meta.LonMax-1)) & (PF(:,2) > (Meta.LatMin+1)) & (PF(:,2) < (Meta.LatMax-1)));  %only add text if the text will stay inside the range of the map
%         text(PF(ipos(1:40:end),1),PF(ipos(1:40:end),2),'PF','FontSize',8)  %label the front 
%         clear ipos
% 
%         SACCF = load(fullfile(Meta.common_loc,'saccf.asc.txt'));                        %load the Southern Atlantic Circumpolar Current Front file
%         plot(SACCF(:,1),SACCF(:,2),'-','Color',[.9 .87 .88],'LineWidth', 1.5)  %plot the longitude/latitude of this front in a tan color
%         %for row = 1:25:length(SACCF)                          %loop through length of SACCF file so that you label every nth point (adjust n to make make more/less cluttered)
%         %    x = SACCF(row,1);                                 %assign x to longitude of SACCF
%         %    y = SACCF(row,2);                                 %assign y to latitude of SACCF
%         %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))     %only add text if the text will stay inside the range of the map
%         %        text(x,y,'SACCF','FontSize', 8)               %labels the fronts
%         %    end
%         %end
%         ipos = find( (SACCF(:,1) > (Meta.LonMin+1)) & (SACCF(:,1) < (Meta.LonMax-1)) & (SACCF(:,2) > (Meta.LatMin+1)) & (SACCF(:,2) < (Meta.LatMax-1))) ; %only add text if the text will stay inside the range of the map
%         text(SACCF(ipos(1:40:end),1),SACCF(ipos(1:40:end),2),'SACCF','FontSize',8)  %label the front only once
%         clear ipos
% 
%         contour(BathFile.x,BathFile.y,BathFile.z',[0 0],'Color','k')                  % outline the coasts with black contour by accessing Bathymetry file lon/lat/dep
%         [L,m] = contour(BathFile.x,BathFile.y,BathFile.z',[-4000 -3000 -2000 -1000],'Color','k'); % -2000m isobath (farthest floats go to) black countour and assign lon/lat positions to L/M
%         clabel(L,m,'LabelSpacing',400);                                               %label contour with "2000" every nth datapoint
% 
%         ylabel('Latitude(degrees)','FontSize',16)    %add y label for Latitude
%         xlabel('Longitude(degrees)','FontSize',16)   %add x label for Longitude
% 
% 
%         % if  strcmp(area,'OOI')      %if your area is OOI
%         %     legend([Ice Feb More Flo], 'Ice on September 20, 2006', 'Ice on February 20, 2006','Station Locations','Mooring','Proposed Floats','Location','Northwest')  %add legend
%         % elseif strcmp(area,'CSIRO') %if your area is CSIRO
%         %     legend([Ice Feb Flo], 'Ice on September 20, 2006', 'Ice on February 20, 2006','Proposed Floats','Location','Northwest')  %add legend
% 
%         hold off    %we are done plotting
%     end 
% 
%     function SA_w_traj(cruise,titl)    
%         figure(3)
%         
%         dd = load([titl,'.mat']);
%         
%         pos2 = [0.48 0.1 0.49 0.85];
%         %left bottom width height
%         subplot('Position',pos2)
%         
%         hold on
%         grid on
%         set(gca,'ydir','reverse')
%         set(gca,'FontSize',22)
%         
%         title([titl,': Absolute Salinity Profiles in Argentine Basin'],'FontSize',22)
%         xlabel('Abs Sal [C]','FontSize',22)
%         ylabel('Pressure [dbars]','FontSize',22)
%         
%         xmin = 0;
%         xmax = (size(dd.data.lat)).^2;
%         dx = size(dd.data.lat);
%         depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
%         nlayers = size(depths,2);
%         cols = hsv(nlayers);
%         
%         for i = 1:size(dd.idat.sa)
%             plot(dd.idat.sa(i,:),dd.idat.P(i,:),'color',cols(i,:))
%         end
%         
%         %For 12700:
%         
%         if titl == '12700'
%             mar4th = plot(dd.idat.sa(17,:),dd.idat.P(17,:),'k-','LineWidth',2);
%         elseif titl == '12778'
%             dec20th = plot(dd.idat.sa(6,:),dd.idat.P(6,:),'k-','LineWidth',2);
%         end
%         
%         %For 12778:
%        
% %         dec20th = plot(dd.idat.DO(106,:),dd.idat.P(106,:),'r-','LineWidth',2)
% %         
% %         mar2nd = plot(dd.idat.DO(113,:),dd.idat.P(113,:),'k-','LineWidth',1.5)
% %         legend([dec20th mar2nd],'Profile taken on Dec 20th','Profile March 2nd');
%         
%         hold off
%         
%         pos1 = [0.05 0.1 0.35 0.85];
%         subplot('Position',pos1)
%         
%          %   %closes all open figure windows 
%         %=========================================
%         Meta = MetaFile(cruise);
% 
%         hold on   %keep figure up throughout the script so that we can keep plotting awesome data over it!
%         title([titl,': Trajectory from ',datestr(dd.data.time(1,:)),' to ',datestr(dd.data.time(end,:))],'FontSize',15)
% 
%         Bath = load(Meta.BathyFile);         %loads Bathymetry file for cruise specified above
%         Bathf = fields(Bath);           %gets the fields from the loaded Bathymetry file
%         BathFile = Bath.([Bathf{:}]);   %creates final bathfile
% 
%         [numrows, numcols] = size(BathFile.z);  %finds the size of the depth part of the Bathymetry file and saves the number of rows/collumns
%         for r = 1:numrows        %double for loop goes through all of the data in the bath file
%             for c = 1:numcols
%                 if (BathFile.z(r,c)) >= 0       %if the elevation is greater or equal to zero
%                     BathFile.z(r,c) = 51;      %set the elevation to 51 (the colorscale only colors it in as land if it is greater than 50 so this makes all the land colored in)
%                 end
%             end
%         end
% 
% 
%         [ahomap,cblev] = AHO_gebco_matlab;      %makes the colormap the gebco map and assigns rgb color code values to ahomap and assigns depths to cblev
%         caxis([cblev(1) cblev(end)])            %makes the axis go from 100 (tan land color) to -6000 (dark blue)
%         cba = colorbar;                         %cba is our colorbar variable
%         % se added comments cba.Label.String = 'Depth (meters)';
%         set(cba,'YTick', cblev(1:end-1))        %set tick marks on cba to show depths from cblev (all except the last the land one (100): 0, -200, -500, -1000, -2000, -3000, -4000, -5000, -6000), these ticks correspond to changing color levels
%         colormap(ahomap)                        %makes the colormap use the ahomap rgb codes
% 
% 
%         [X,Y] = meshgrid(BathFile.x,BathFile.y);      %makes grid using longitudes and latitudes from the Bathymetry file
%         imagesc(BathFile.x,BathFile.y,BathFile.z');   %plots image using corresponding lon/lat and the depths from the Bathymetry file (rotated so that they match up correctly)
%         set(gca,'YDir','normal');        %gets handle of current axis (gca) and sets it to a normal y-direction
%         set(gca,'YLim',[Meta.LatMin Meta.LatMax]); % sets map y-limits to the Latitudes specified for your area
%         set(gca,'XLim',[Meta.LonMin Meta.LonMax]); % sets map x-limits to the Longitudes specified for you area
%         set(gca,'FontSize',18)
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %==================================================================
%         d = load([titl,'.mat']);
% 
%         Latitude = dd.data.lat;   %get Latitude variable
%         Longitude = dd.data.lon;
%         Cruise = line(Longitude(1:end),Latitude(1:end),'LineWidth',2,'Color',[0 0 0]);
%         
%         xmin = 0;
%         xmax = (size(dd.data.lat)).^2;
%         dx = size(dd.data.lat);
%         depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
%         nlayers = size(depths,2);
%         cols = hsv(nlayers);
%         
%         for row = 1:size(dd.data.lat)    %loop through the positions in the file
%              %get Longitude variable
%             Flo = plot(Longitude(row),Latitude(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(row,:)); 
%             if titl == '12700'
%                 mar4th = plot(Longitude(17,:),Latitude(17,:),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
%             elseif titl == '12778'
%                 dec20th = plot(Longitude(6,:),Latitude(6,:),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);;
%             end 
%         end
%         
%         %For 12700:
% %         for row = 17 %length(d.data.lat)    %loop through the positions in the file
% %              %get Longitude variable
% %             mar2nd = plot(Longitude(row),Latitude(row),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
% %             %WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);
% %         end
%         
%         
% %         for row = 113 %length(d.data.lat)    %loop through the positions in the file
% %              %get Longitude variable
% %             mar2nd = plot(Longitude(row),Latitude(row),'.','MarkerSize',2.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]); 
% %             %WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]); 
% %         end
% %         
% %         for row = 106 %length(d.data.lat)    %loop through the positions in the file
% %              %get Longitude variable
% %             dec20th = plot(Longitude(row),Latitude(row),'.','MarkerSize',2.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]); 
% %             %WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]); 
% %         end
% 
% 
%         %legend([WeirdFlo WeirdFlo2],'Profile 17 of 12700','Profile 6 of 12788')
% 
%         %==================================================================
%         
%         ncfname   = fullfile(Meta.common_loc,Meta.SeaIceFile1); %gets ice file by calling fullfile(FolderName,'FileName') and referring to variables specified at the top
%         ncid = netcdf.open(ncfname,'NC_NOWRITE');     %opens ice file so that you can't edit the file and assigns file to "ncid"
%         [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);   %keeps track of netcdf variables
%         for i = 1:numvars                             %loops through variables
%             [varnames{i}, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,i-1);   %saves variables
%             Data.([varnames{i}]) = netcdf.getVar(ncid,i-1);   %finds what variables are and saves all of them into Data so that you can access them
%         end
%         netcdf.close(ncid);  %closes netcdf file
%         [numrows,numcols] = size(Data.longitude);  %finds size of longitude,latitude,and Ice Conc arrays (they will each be same size) and assigns those values to numrows and numcols
%         for r = 1:numrows  %double for loop goes through each row and collumn or Ice Concentration array
%             for c = 1:numcols  %goes through all the collumns
%                 if (isnan(Data.IceConc(r,c)) == 1) || (Data.IceConc(r,c) == 0) %checks for NaN data and puts 0 in instead
%                     Data.IceConc(r,c)=0;
%                 elseif ((Data.longitude(r,c) > 180) && (Data.longitude(r,c) < 360)) == 1   %if longitude values are between 180 and 360 in the ice file, subtract 360 from them to make the value match up to our -180 to 180 scale
%                     Ice = plot(Data.longitude(r,c)-360, Data.latitude(r,c),'o','MarkerSize',4.5,'LineWidth',1,'MarkerEdgeColor',[0 1 1]); %plot where the ice is at the correct lat lon with cyan circles (keep middle of cirles clear so that you can see bathymetry levels) and assign to "ice" so that we can use in the legend
%                 else   %if the ice location is between 0 and 180 degrees east
%                     Ice = plot(Data.longitude(r,c), Data.latitude(r,c),'o','MarkerSize',4.5,'LineWidth',1,'MarkerEdgeColor',[0 1 1]); %plot where the ice is at the correct lat lon with cyan circles (keep middle of circles clear so that you can see bathymetry) and assign to "ice" so that we can use in legend
%                 end
%             end
%         end
% 
%         ncfname   = fullfile(Meta.common_loc,Meta.SeaIceFile2); %gets second ice file by calling fullfile(FolderName,'FileName') and referring to variables specified above
%         ncid = netcdf.open(ncfname,'NC_NOWRITE');     %opens ice file so that you can't edit the file and assigns file to "ncid"
%         [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);  %keeps track of netcdf variables
%         for i = 1:numvars                                                  %loops through variables
%             [varnames{i}, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,i-1);    %saves variables
%             Data.([varnames{i}]) = netcdf.getVar(ncid,i-1);   %finds what variables are and saves all of them into Data so that you can access them
%         end
%         netcdf.close(ncid);
%         [numrows,numcols] = size(Data.longitude);  %finds size of longitude,latitude,and Ice Conc arrays (they will each be same size) and assigns those values to numrows and numcols
%         for r = 1:numrows  %double for loop goes through each row and each collumn
%             for c = 1:numcols
%                 if (isnan(Data.IceConc(r,c)) == 1) || (Data.IceConc(r,c) == 0) %checks for NaN data and puts 0 in instead
%                     Data.IceConc(r,c)=0;
%                 elseif ((Data.longitude(r,c) > 180) && (Data.longitude(r,c) < 360)) == 1    %if longitude values are between 180 and 360 in the ice file, subtract 360 from them to make the value match up to our -180 to 180 scale
%                     Feb = plot(Data.longitude(r,c)-360, Data.latitude(r,c),'^','MarkerSize',5,'LineWidth',.45,'MarkerEdgeColor',[0 0 1]); %plot where the ice is at the correct lat lon with a blue triangle (middle open so that you can see bathymetry) and assign to "Feb" so that we can use in legend
%                 else              %if longitude values for ice are between 0 and 180 east
%                     Feb = plot(Data.longitude(r,c), Data.latitude(r,c),'^','MarkerSize',5,'LineWidth',.45,'MarkerEdgeColor',[0 0 1]);    %plot where the ice is at the correct lat lon with blue triangle (middle open so that you can see bathymetry) and assign to "Feb so that we can use in legend
%                 end
%             end
%         end
% 
% 
% 
%         SAF = load(fullfile(Meta.common_loc,'saf.asc.txt'));                            %load Subantarctic Front file
%         plot(SAF(:,1),SAF(:,2),'-','Color',[.9 .87 .88],'LineWidth', 1.5)      %plot the longitude/latitude of this front in a tan color
%         % This code does not work:  see working code below commented lines
%         %for row = 1:80:length(SAF)                            %loop through length of SAF file so that you label every nth point (adjust n to make make more/less cluttered)
%         %    x = SAF(row,1);                                   %assign x to longitude of SAF
%         %    y = SAF(row,2);                                   %ssign y to latitude of SAF
%         %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))  %only add text if the text will stay inside the range of the map
%         %        text(x,y,'SAF','FontSize',8)                  %labels the fronts
%         %    end
%         %end
%         % use below instead
%         %    x = SAF(row,1);                                   %assign x to longitude of SAF
%         %    y = SAF(row,2);                                   %ssign y to latitude of SAF
%         ipos = find( (SAF(:,1) > (Meta.LonMin+1)) & (SAF(:,1) < (Meta.LonMax-1)) & (SAF(:,2) > (Meta.LatMin+1)) & (SAF(:,2) < (Meta.LatMax-1)));  %only add text if the text will stay inside the range of the map
%         text(SAF(ipos(1:40:end),1),SAF(ipos(1:40:end),2),'SAF','FontSize',8)  %label the front only once
%         clear ipos
% 
%         SBDY = load(fullfile(Meta.common_loc,'sbdy.asc.txt'));                           %load the southern boundary file
%         plot(SBDY(:,1),SBDY(:,2),'-','Color',[.9 .87 .88],'LineWidth', 1.5)     %plot the longitude/latitude of this in a tan color
%         %for row = 1:18:length(SBDY)                            %loop through length of SBDY file so that you label every nth point (adjust n to make make more/less cluttered)
%         %    x = SBDY(row,1);                                   %assign x to longitude of SBDY
%         %    y = SBDY(row,2);                                   %assign y to latitude of SBDY
%         %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))  %only add text if the text will stay inside the range of the map
%         %        text(x,y,'SBDY','FontSize',8)                  %add text
%         %    end
%         %end
%         ipos = find( (SBDY(:,1) > (Meta.LonMin+1)) & (SBDY(:,1) < (Meta.LonMax-1)) & (SBDY(:,2) > (Meta.LatMin+1)) & (SBDY(:,2) < (Meta.LatMax-1))) ; %only add text if the text will stay inside the range of the map
%         text(SBDY(ipos(1:40:end),1),SBDY(ipos(1:40:end),2),'SBDY','FontSize',8)  %label the front only once
%         clear ipos
% 
%         STF = load(fullfile(Meta.common_loc,'stf.asc.txt'));                              %loads Subtropical Front file (only appears on a few maps because it is farther north)
%         plot(STF(:,1),STF(:,2),'-','Color', [.9 .87 .88],'LineWidth', 1.5)       %plot the longitude/latitude of this in a tan color
%         %for row = 1:30:length(STF)                              %loop through length of STF file so that you label every nth point (adjust n to make make more/less cluttered)
%         %    x = STF(row,1);                                     %assign x to longitude of STF
%         %    y = STF(row,2);                                     %assign y to latitude of STF
%         %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))   %only add text if the text will stay inside the range of the map
%         %        text(x,y,'STF','FontSize',8)                    %add text
%         %    end
%         %end
%         ipos = find( (STF(:,1) > (Meta.LonMin+1)) & (STF(:,1) < (Meta.LonMax-1)) & (STF(:,2) > (Meta.LatMin+1)) & (STF(:,2) < (Meta.LatMax-1))) ; %only add text if the text will stay inside the range of the map
%         if ipos    
%            text(STF(ipos(1:40:end),1),STF(ipos(1:40:end),2),'STF','FontSize',8)  %label the front only once
%         end
%         clear ipos
% 
%         PF = load(fullfile(Meta.common_loc,'pf.asc.txt'));                               %loads Polar Front file
%         plot(PF(:,1),PF(:,2),'-','Color',[.9 .87 .88],'LineWidth', 1.5)         %plot the longitude/latitude of this in a tan color
%         %for row = 1:90:length(PF)                              %loop through length of PF file so that you label every nth point (adjust n to make make more/less cluttered)
%         %    x = PF(row,1);                                     %assign x to longitude of STF file
%         %    y = PF(row,2);                                     %assign y to latitude of STF file
%         %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))  %only add text if the text will stay inside the range of the map
%         %        text(x,y,'PF','FontSize',8)                    %add text
%         %    end
%         %end
%         ipos = find( (PF(:,1) > (Meta.LonMin+1)) & (PF(:,1) < (Meta.LonMax-1)) & (PF(:,2) > (Meta.LatMin+1)) & (PF(:,2) < (Meta.LatMax-1)));  %only add text if the text will stay inside the range of the map
%         text(PF(ipos(1:40:end),1),PF(ipos(1:40:end),2),'PF','FontSize',8)  %label the front 
%         clear ipos
% 
%         SACCF = load(fullfile(Meta.common_loc,'saccf.asc.txt'));                        %load the Southern Atlantic Circumpolar Current Front file
%         plot(SACCF(:,1),SACCF(:,2),'-','Color',[.9 .87 .88],'LineWidth', 1.5)  %plot the longitude/latitude of this front in a tan color
%         %for row = 1:25:length(SACCF)                          %loop through length of SACCF file so that you label every nth point (adjust n to make make more/less cluttered)
%         %    x = SACCF(row,1);                                 %assign x to longitude of SACCF
%         %    y = SACCF(row,2);                                 %assign y to latitude of SACCF
%         %    if (x > (Meta.LonMin+1) && x < (Meta.LonMax-1)) && (y > (Meta.LatMin+1) && y < (Meta.LatMax-1))     %only add text if the text will stay inside the range of the map
%         %        text(x,y,'SACCF','FontSize', 8)               %labels the fronts
%         %    end
%         %end
%         ipos = find( (SACCF(:,1) > (Meta.LonMin+1)) & (SACCF(:,1) < (Meta.LonMax-1)) & (SACCF(:,2) > (Meta.LatMin+1)) & (SACCF(:,2) < (Meta.LatMax-1))) ; %only add text if the text will stay inside the range of the map
%         text(SACCF(ipos(1:40:end),1),SACCF(ipos(1:40:end),2),'SACCF','FontSize',8)  %label the front only once
%         clear ipos
% 
%         contour(BathFile.x,BathFile.y,BathFile.z',[0 0],'Color','k')                  % outline the coasts with black contour by accessing Bathymetry file lon/lat/dep
%         [L,m] = contour(BathFile.x,BathFile.y,BathFile.z',[-4000 -3000 -2000 -1000],'Color','k'); % -2000m isobath (farthest floats go to) black countour and assign lon/lat positions to L/M
%         clabel(L,m,'LabelSpacing',400);                                               %label contour with "2000" every nth datapoint
% 
%         ylabel('Latitude(degrees)','FontSize',16)    %add y label for Latitude
%         xlabel('Longitude(degrees)','FontSize',16)   %add x label for Longitude
% 
% 
%         % if  strcmp(area,'OOI')      %if your area is OOI
%         %     legend([Ice Feb More Flo], 'Ice on September 20, 2006', 'Ice on February 20, 2006','Station Locations','Mooring','Proposed Floats','Location','Northwest')  %add legend
%         % elseif strcmp(area,'CSIRO') %if your area is CSIRO
%         %     legend([Ice Feb Flo], 'Ice on September 20, 2006', 'Ice on February 20, 2006','Proposed Floats','Location','Northwest')  %add legend
% 
%         hold off    %we are done plotting
    %end 
end 