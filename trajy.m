function trajy(cruise,titl)


clf 

        %   %closes all open figure windows 
        %=========================================
        Meta = MetaFile(cruise);
        (genpath('.../Research/Profile17/09624'));
        (genpath('.../Research/Profile17/09750'));
        (genpath('.../Research/Profile17/12379'));
        (genpath('.../Research/Profile17/12386'));
        (genpath('.../Research/Profile17/12542'));
        (genpath('.../Research/Profile17/12543'));
        (genpath('.../Research/Profile17/12545'));
        
        (genpath('.../Research/09646'));
        (genpath('.../Research/12573'));
        (genpath('.../Research/12575'));
        (genpath('.../Research/12696'));
        (genpath('.../Research/12700'));
        (genpath('.../Research/12778'));
        (genpath('.../Research/12881'));
        (genpath('.../Research/12883'));
        
%        Meta = MetaInfo(cruise);
        % % these may be needed for previous cruises
        % parts = strread(cruise,'%s','delimiter','_');
        % year = parts{1};
        % area = parts{2};

        %below now in MetaInfo.m
        %common_loc = '~/work/sogle_plots/common_files/';     % path to common soccom files
        %SeaIceFile1 = 'nt_20060920_f13_v01_s.nc';       %greatest ice extent (change this depending on what sea ice file you want (this is September 20, 2006))
        %SeaIceFile2 = 'nt_20060220_f13_v01_s.nc';       %lower ice extent (change this depending on what sea ice file you want (this is February 20, 2006))
        %=========================================

        %EDIT variables below
        %If you want to add a new area add an elseif with this code:
        %elseif strcmp(area,'your new area') 
        %    Station = 'name of your station file'  %file must have 5 collums with
        %        station number in the first collumn, latitude in the 3rd collumn,
        %        and longitude in the 4th column (lat/lon in decimal form -180 to
        %        180 form) 
        %        text file
        %        if there isn't a station file set Station to "false"
        %    Floats = 'name of float file'       %file must have  at least 3 collumns with
        %        float number in the first collumn,latitude in the 2nd collumn and
        %        longitude in the 3rd collumn (lat/lon in decimal -180 to 180 form) 
        %        text file
        %        if you don't have a float file, then assign Floats to "false"
        %    BathyFile = 'name of bathymetry file'   %file must be for the
        %        specific location in .mat format with ocean depths as negative #s
        %    LatMin = #; LatMax = #; LonMin = #; LonMax = #; (sets boundary for
        %        graph (lat/lon must be in -180 to 180 scale/decimal form)
        %    title('input any title you want')
        %now call the function by typing LatLonDep('area') in your command window

        hold on   %keep figure up throughout the script so that we can keep plotting awesome data over it!
        title(['Trajectories of Argentine Basin with DO'],'FontSize',15)

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
        set(gca,'FontSize',16)
        
        
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

        contour(BathFile.x,BathFile.y,BathFile.z',[0 0],'Color','k')                  % outline the coasts with black contour by accessing Bathymetry file lon/lat/dep
        [L,m] = contour(BathFile.x,BathFile.y,BathFile.z',[-4000 -3000 -2000 -1000],'Color','k'); % -2000m isobath (farthest floats go to) black countour and assign lon/lat positions to L/M
        clabel(L,m,'LabelSpacing',400); 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %==================================================================        
        d = load('09642.mat');
        
        xmin = 0;
        xmax = 32;
        dx = 2;
        depths = fliplr(-[xmin:dx:xmax-dx; xmin+dx:dx:xmax]);
        nlayers = size(depths,2);
        cols = hsv(nlayers);
        
        Latitude = d.data.lat;   %get Latitude variable
        Longitude = d.data.lon;  %get Longitude variable
        
        Cruise = line(Longitude,Latitude,'LineWidth',2,'Color',cols(1,:));

        for row = 1:1:length(d.data.lat)    %loop through the positions in the file
            Flo = plot(Longitude(row),Latitude(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(16,:)); 
            %WeirdFlo2 = plot(Longitude(6),Latitude(6),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d1 = load('09750.mat');
        
        Latitude1 = d1.data.lat;   %get Latitude variable
        Longitude1 = d1.data.lon;  %get Longitude variable
        
        Cruise1 = line(Longitude1,Latitude1,'LineWidth',2,'Color',cols(2,:));

        for row = 1:1:length(d1.data.lat)    %loop through the positions in the file
            Flo1 = plot(Longitude1(row),Latitude1(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(2,:)); 
            %WeirdFlo2 = plot(Longitude2(6),Latitude2(6),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end

        
        d2 = load('12379.mat');
        
        Latitude2 = d2.data.lat;   %get Latitude variable
        Longitude2 = d2.data.lon;  %get Longitude variable
        Cruise2 = line(Longitude2,Latitude2,'LineWidth',2,'Color',cols(3,:));

        for row = 1:1:length(d2.data.lat)    %loop through the positions in the file
           Flo2 = plot(Longitude2(row),Latitude2(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(3,:)); 
           %WeirdFlo3 = plot(Longitude3(17),Latitude3(17),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d3 = load('12386.mat');
        
        Latitude3 = d3.data.lat;   %get Latitude variable
        Longitude3 = d3.data.lon;  %get Longitude variable
        Cruise3 = line(Longitude3,Latitude3,'LineWidth',2,'Color',cols(4,:));

        for row = 1:1:length(d3.data.lat)    %loop through the positions in the file
           Flo3 = plot(Longitude3(row),Latitude3(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(4,:)); 
           %WeirdFlo3 = plot(Longitude3(17),Latitude3(17),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d4 = load('12542.mat');
        
        Latitude4 = d4.data.lat;   %get Latitude variable
        Longitude4 = d4.data.lon;  %get Longitude variable
        Cruise4 = line(Longitude4,Latitude4,'LineWidth',2,'Color',cols(5,:));

        for row = 1:1:length(d4.data.lat)    %loop through the positions in the file
           Flo4 = plot(Longitude4(row),Latitude4(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(5,:)); 
           %WeirdFlo3 = plot(Longitude3(17),Latitude3(17),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d5 = load('12543.mat');
        
        Latitude5 = d5.data.lat;   %get Latitude variable
        Longitude5 = d5.data.lon;  %get Longitude variable
        Cruise5 = line(Longitude5,Latitude5,'LineWidth',2,'Color',cols(6,:));

        for row = 1:1:length(d5.data.lat)    %loop through the positions in the file
           Flo5 = plot(Longitude5(row),Latitude5(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(6,:)); 
           %WeirdFlo3 = plot(Longitude3(17),Latitude3(17),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d6 = load('12545.mat');
        
        Latitude6 = d6.data.lat;   %get Latitude variable
        Longitude6 = d6.data.lon;  %get Longitude variable
        Cruise6 = line(Longitude6,Latitude6,'LineWidth',2,'Color',cols(7,:));

        for row = 1:1:length(d6.data.lat)    %loop through the positions in the file
           Flo6 = plot(Longitude6(row),Latitude6(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(7,:)); 
           %WeirdFlo3 = plot(Longitude3(17),Latitude3(17),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d7 = load('09668.mat');
        
        Latitude7 = d7.data.lat;   %get Latitude variable
        Longitude7 = d7.data.lon;  %get Longitude variable
        Cruise7 = line(Longitude7,Latitude7,'LineWidth',2,'Color',cols(8,:));

        for row = 1:1:length(d7.data.lat)    %loop through the positions in the file
           Flo7 = plot(Longitude7(row),Latitude7(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(8,:)); 
           %WeirdFlo3 = plot(Longitude3(17),Latitude3(17),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d8 = load('12369.mat');
        
        Latitude8 = d8.data.lat;   %get Latitude variable
        Longitude8 = d8.data.lon;  %get Longitude variable
        Cruise8 = line(Longitude8,Latitude8,'LineWidth',2,'Color',cols(9,:));

        for row = 1:1:length(d8.data.lat)    %loop through the positions in the file
           Flo8 = plot(Longitude8(row),Latitude8(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(9,:)); 
           %WeirdFlo3 = plot(Longitude3(17),Latitude3(17),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d9 = load('12540.mat');
        
        Latitude9 = d9.data.lat;   %get Latitude variable
        Longitude9 = d9.data.lon;  %get Longitude variable
        Cruise9 = line(Longitude9,Latitude9,'LineWidth',2,'Color',cols(10,:));

        for row = 1:1:length(d9.data.lat)    %loop through the positions in the file
           Flo9 = plot(Longitude9(row),Latitude9(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(10,:)); 
           %WeirdFlo3 = plot(Longitude3(17),Latitude3(17),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d10 = load('12573.mat');
        
        Latitude10 = d10.data.lat;   %get Latitude variable
        Longitude10 = d10.data.lon;  %get Longitude variable
        Cruise10 = line(Longitude10,Latitude10,'LineWidth',2,'Color',cols(11,:));

        for row = 1:1:length(d10.data.lat)    %loop through the positions in the file
           Flo10 = plot(Longitude10(row),Latitude10(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(11,:)); 
           %WeirdFlo3 = plot(Longitude3(17),Latitude3(17),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d11 = load('12575.mat');
        
        Latitude11 = d11.data.lat;   %get Latitude variable
        Longitude11 = d11.data.lon;  %get Longitude variable
        Cruise11 = line(Longitude11(2:86),Latitude11(2:86),'LineWidth',2,'Color',cols(12,:));

        for row = 1:1:length(d11.data.lat)    %loop through the positions in the file
           Flo11 = plot(Longitude11(row),Latitude11(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(12,:)); 
           %WeirdFlo3 = plot(Longitude3(17),Latitude3(17),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d12 = load('12700.mat');
        
        Latitude12 = d12.data.lat;   %get Latitude variable
        Longitude12 = d12.data.lon;  %get Longitude variable
        Cruise12 = line(Longitude12,Latitude12,'LineWidth',2,'Color',cols(13,:));

        for row = 1:1:length(d12.data.lat)    %loop through the positions in the file
           Flo12 = plot(Longitude12(row),Latitude12(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(13,:)); 
           WeirdFlo12 = plot(Longitude12(17),Latitude12(17),'.','MarkerSize',2.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d13 = load('12778.mat');
        
        Latitude13 = d13.data.lat;   %get Latitude variable
        Longitude13 = d13.data.lon;  %get Longitude variable
        Cruise13 = line(Longitude13,Latitude13,'LineWidth',2,'Color',cols(14,:));

        for row = 1:1:length(d13.data.lat)    %loop through the positions in the file
           Flo13 = plot(Longitude13(row),Latitude13(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(14,:)); 
           WeirdFlo13 = plot(Longitude13(6),Latitude13(6),'.','MarkerSize',2.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d14 = load('12881.mat');
        
        Latitude14 = d14.data.lat;   %get Latitude variable
        Longitude14 = d14.data.lon;  %get Longitude variable
        Cruise14 = line(Longitude14,Latitude14,'LineWidth',2,'Color',cols(15,:));

        for row = 1:1:length(d14.data.lat)    %loop through the positions in the file
           Flo14 = plot(Longitude14(row),Latitude14(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',cols(15,:)); 
           %WeirdFlo3 = plot(Longitude3(17),Latitude3(17),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        d15 = load('09646.mat');
        
        Latitude15 = d15.data.lat;   %get Latitude variable
        Longitude15 = d15.data.lon;  %get Longitude variable
        Cruise15 = line(Longitude15,Latitude15,'LineWidth',2,'Color',[0 0 0]);

        for row = 1:1:length(d15.data.lat)    %loop through the positions in the file
           Flo15 = plot(Longitude15(row),Latitude15(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]); 
           %WeirdFlo3 = plot(Longitude3(17),Latitude3(17),'.','MarkerSize',1.5*Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0 0 0]);
        end
        
        legend([Flo Flo2 Flo3 Flo4 Flo5 Flo6 Flo7 Flo8 Flo9 Flo9 Flo10 Flo11 Flo12 Flo13 Flo14 Flo15],'9642','9750','12379','12386',...
            '12542','12543','12545','9668','12369','9646','12573','12575','12700','12778','12881','12540')

        %==================================================================
        %label contour with "2000" every nth datapoint

        ylabel('Latitude(degrees)','FontSize',16)    %add y label for Latitude
        xlabel('Longitude(degrees)','FontSize',16)   %add x label for Longitude


        % if  strcmp(area,'OOI')      %if your area is OOI
        %     legend([Ice Feb More Flo], 'Ice on September 20, 2006', 'Ice on February 20, 2006','Station Locations','Mooring','Proposed Floats','Location','Northwest')  %add legend
        % elseif strcmp(area,'CSIRO') %if your area is CSIRO
        %     legend([Ice Feb Flo], 'Ice on September 20, 2006', 'Ice on February 20, 2006','Proposed Floats','Location','Northwest')  %add legend

        hold off    %we are done plotting
end