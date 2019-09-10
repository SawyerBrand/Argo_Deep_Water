function Alty(cruise,titl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was made by Sawyer Brand, SIO/SOCCOM            %
% during research into the Argentine Basin mesoscale activities %
% and CO2/Heat Flux in that region                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Meta = MetaFile(cruise);
%meta = MetaInfo(cruise);

A = load(Meta.AltimFile);  %links to the AltimFile within the MetaInfo file 

parts = strread(cruise,'%s','delimiter','_');
year = parts{1};
area = parts{2};

clf;

figure(1)

hold on 

title(['Mean Aviso Altimetry for Drake Passage'],'FontSize',20)%,titl],'FontSize',16) %creates title
ylabel('Latitude (degrees)','FontSize',16)    %add y label for Latitude
xlabel('Longitude (degrees)','FontSize',16)   %add x label for Longitude
imagesc( A.lon,A.lat,A.avdepth) % creates the basic image we want to edit
set(gca,'YLim',[Meta.LatMin Meta.LatMax]); % sets map y-limits to the Latitudes specified for your area
set(gca,'XLim',[Meta.LonMin Meta.LonMax]); % sets map x-limits to the Latitudes specified for your area
set(gca,'FontSize',18)

contour( A.lon,A.lat,A.avdepth,[0 0],'Color','k')                  % outline the coasts with black contour by accessing Bathymetry file lon/lat/dep
[L,m] = contour( A.lon,A.lat,A.avdepth,[-0.4:0.1:0.8],'Color','k'); % -2000m isobath (farthest floats go to) black countour and assign lon/lat positions to L/M
clabel(L,m,'LabelSpacing',200); %labels the contours

% %==================================================================
% d = load('12700.mat');
% 
% Latitude = d.data.lat;   %get Latitude variable
% Longitude = d.data.lon;  %get Longitude variable
% Cruise = line(Longitude,Latitude,'LineWidth',2,'Color',[0 0 1]);
% 
% for row = 1:1:size(d.data.lat,1)    %loop through the positions in the file
%     Flo = plot(Longitude(row),Latitude(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0.5 0 0]); 
% end
% 
% WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);


% d2 = load('12881dat mat');
% 
% for row2 = 1:1:size(d2.dat lat,1)    %loop through the positions in the file
%     Latitude2 = d2.dat lat;   %get Latitude variable
%     Longitude2 = d2.dat lon;  %get Longitude variable
%     Flo2 = plot(Longitude2(row2),Latitude2(row2),'.','MarkerSize',met MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0.5 0 0]); 
%     WeirdFlo = plot(Longitude2(17),Latitude2(17),'.','MarkerSize',met MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);
% end
% 
% Cruise2 = line(Longitude2,Latitude2,'LineWidth',2,'Color',[0 0 0.5]);

%==================================================================

coast = load('coast'); %loads coast
geoshow ('landareas.shp','FaceColor','black'); %creates continent on map

h = colorbar; %creates colorbar
set(get(h,'label'),'string','[m]','FontSize',16) %labels colorbar


%legend([Cruise Flo WeirdFlo],'Track 12700','Float Locations 12700','Profile 17 of 12700','location','northwest') %MUST CHANGE PER CRUISE

hold off

figure(2)

hold on 

title(['Aviso Standard Deviation (2017-2018) for ',titl],'FontSize',16) %creates title
ylabel('Latitude (degrees)','FontSize',16)    %add y label for Latitude
xlabel('Longitude (degrees)','FontSize',16)   %add x label for Longitude
imagesc( A.lon,A.lat,A.stddev) % creates the basic image we want to edit
set(gca,'YLim',[Meta.LatMin Meta.LatMax]); % sets map y-limits to the Latitudes specified for your area
set(gca,'XLim',[Meta.LonMin Meta.LonMax]); % sets map x-limits to the Latitudes specified for your area
set(gca,'FontSize',18)

contour( A.lon,A.lat,A.stddev,[0 0],'Color','k')                  % outline the coasts with black contour by accessing Bathymetry file lon/lat/dep
[L,m] = contour( A.lon,A.lat,A.stddev,[-0.4:0.05:0.5],'Color','k'); % -2000m isobath (farthest floats go to) black countour and assign lon/lat positions to L/M
clabel(L,m,'LabelSpacing',200); %labels the contours
%==================================================================
% d = load('12700.mat');
% 
% Latitude = d.data.lat;   %get Latitude variable
% Longitude = d.data.lon;  %get Longitude variable
% Cruise = line(Longitude,Latitude,'LineWidth',2,'Color',[0 0 1]);
% 
% for row = 1:1:size(d.data.lat,1)    %loop through the positions in the file
%     Flo = plot(Longitude(row),Latitude(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0.5 0 0]); 
% end
% 
% WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);


% d2 = load('12881dat mat');
% 
% for row2 = 1:1:size(d2.dat lat,1)    %loop through the positions in the file
%     Latitude2 = d2.dat lat;   %get Latitude variable
%     Longitude2 = d2.dat lon;  %get Longitude variable
%     Flo2 = plot(Longitude2(row2),Latitude2(row2),'.','MarkerSize',met MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0.5 0 0]); 
%     WeirdFlo = plot(Longitude2(17),Latitude2(17),'.','MarkerSize',met MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);
% end
% 
% Cruise2 = line(Longitude2,Latitude2,'LineWidth',2,'Color',[0 0 0.5]);

%==================================================================

coast = load('coast'); %loads coast
geoshow ('landareas.shp','FaceColor','black'); %creates continent on map

h = colorbar; %creates colorbar
set(get(h,'label'),'string','[m]','FontSize',16) %labels colorbar

%legend([Cruise Flo WeirdFlo],'Track 12700','Float Locations 12700','Profile 17 of 12700','location','northwest') %MUST CHANGE PER CRUISE

hold off

% figure(3)
% 
% hold on 
% 
% title(['Aviso on March 1st 2019 Divided by Standard Deviation (2017-2018) for ',titl],'FontSize',16) %creates title
% ylabel('Latitude (degrees)','FontSize',16)    %add y label for Latitude
% xlabel('Longitude (degrees)','FontSize',16)   %add x label for Longitude
% imagesc( lon, lat, depth./ stddev) % creates the basic image we want to edit
% set(gca,'YLim',[Met LatMin Met LatMax]); % sets map y-limits to the Latitudes specified for your area
% set(gca,'XLim',[Met LonMin Met LonMax]); % sets map x-limits to the Latitudes specified for your area
% set(gca,'FontSize',18)
% 
% contour( lon, lat, depth./ stddev,[0 0],'Color','k')                  % outline the coasts with black contour by accessing Bathymetry file lon/lat/dep
% [L,m] = contour( lon, lat, depth./ stddev,[-26:2:26],'Color','k'); % -2000m isobath (farthest floats go to) black countour and assign lon/lat positions to L/M
% clabel(L,m,'LabelSpacing',200); %labels the contours
% 
% %==================================================================
% d = load('12700.mat');
% 
% Latitude = d.dat lat;   %get Latitude variable
% Longitude = d.dat lon;  %get Longitude variable
% Cruise = line(Longitude,Latitude,'LineWidth',2,'Color',[0 0 1]);
% 
% for row = 1:1:size(d.dat lat,1)    %loop through the positions in the file
%     Flo = plot(Longitude(row),Latitude(row),'.','MarkerSize',met MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0.5 0 0]); 
% end
% 
% WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',met MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);
% 
% 
% % d2 = load('12881dat mat');
% % 
% % for row2 = 1:1:size(d2.dat lat,1)    %loop through the positions in the file
% %     Latitude2 = d2.dat lat;   %get Latitude variable
% %     Longitude2 = d2.dat lon;  %get Longitude variable
% %     Flo2 = plot(Longitude2(row2),Latitude2(row2),'.','MarkerSize',met MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0.5 0 0]); 
% %     WeirdFlo = plot(Longitude2(17),Latitude2(17),'.','MarkerSize',met MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);
% % end
% % 
% % Cruise2 = line(Longitude2,Latitude2,'LineWidth',2,'Color',[0 0 0.5]);
% 
% %==================================================================
% 
% coast = load('coast'); %loads coast
% geoshow ('landareas.shp','FaceColor','black'); %creates continent on map
% 
% h = colorbar; %creates colorbar
% set(get(h,'label'),'string','[m]','FontSize',16) %labels colorbar
% 
% legend([Cruise Flo WeirdFlo],'Track 12700','Float Locations 12700','Profile 17 of 12700','location','northwest') %MUST CHANGE PER CRUISE
% 
% hold off

%%%%%%%%%%%

% meta = MetaFile(cruise);
% 
% % vidObj = VideoWriter('SSH.avi');
% vidObj = VideoWriter('SSH.mp4','MPEG-4');	% on linux, MATLAB can only make .avi movies
% vidObj.Quality = 100;				% I think this should work on a mac
% vidObj.FrameRate = 6;				% this sets the frame rate
% open(vidObj);
% writeVideo(vidObj, getframe(gcf));
% 
% for i = 1:1:size( dsave)
%     LatInd = find(Alt.ysave <= met LatMax & Alt.ysave >= met LatMin); %finds lat indices within max/min bounds
%     LonInd = find(Alt.xsave <= met LonMax & Alt.xsave >= met LonMin); %finds lon indices within max/min bounds
%     
%     latt = Alt.ysave(LatInd);
%     lonn = Alt.xsave(LonInd);
%     
%     depth = squeeze(Alt.dsave(i,LatInd,LonInd))
%     imagesc( lon, lat, depth);
%     c = colorbar('eastoutside');
%     hold on
%     caxis([lb ub])
%     title('SSH anomaly')
%     xlabel('Latitude')
%     ylabel('Longitude')
%     xtickformat('degrees')
%     ytickformat('degrees')
%     acc_movie
%     hold off
% 
%     drawnow()				% these two commands save the current plot as the next
%     writeVideo(vidObj, getframe(gcf));	% frame in the movie
% end
% 
% close(vidObj);

figure(4)

hold on 

title(['Aviso SLA on March 1st 2019 for ',titl],'FontSize',16) %creates title
ylabel('Latitude (degrees)','FontSize',16)    %add y label for Latitude
xlabel('Longitude (degrees)','FontSize',16)   %add x label for Longitude
imagesc( A.lon,A.lat,A.anom) % creates the basic image we want to edit
set(gca,'YLim',[Meta.LatMin Meta.LatMax]); % sets map y-limits to the Latitudes specified for your area
set(gca,'XLim',[Meta.LonMin Meta.LonMax]); % sets map x-limits to the Latitudes specified for your area
set(gca,'FontSize',18)

contour( A.lon,A.lat,A.anom,[0 0],'Color','k')                  % outline the coasts with black contour by accessing Bathymetry file lon/lat/dep
[L,m] = contour( A.lon,A.lat,A.anom,[-0.4:0.1:0.8],'Color','k'); % -2000m isobath (farthest floats go to) black countour and assign lon/lat positions to L/M
clabel(L,m,'LabelSpacing',200); %labels the contours

%==================================================================
% d = load('12700.mat');
% 
% Latitude = d.data.lat;   %get Latitude variable
% Longitude = d.data.lon;  %get Longitude variable
% Cruise = line(Longitude,Latitude,'LineWidth',2,'Color',[0 0 1]);
% 
% for row = 1:1:size(d.data.lat,1)    %loop through the positions in the file
%     Flo = plot(Longitude(row),Latitude(row),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0.5 0 0]); 
% end
% 
% WeirdFlo = plot(Longitude(17),Latitude(17),'.','MarkerSize',Meta.MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);
% 

% d2 = load('12881dat mat');
% 
% for row2 = 1:1:size(d2.dat lat,1)    %loop through the positions in the file
%     Latitude2 = d2.dat lat;   %get Latitude variable
%     Longitude2 = d2.dat lon;  %get Longitude variable
%     Flo2 = plot(Longitude2(row2),Latitude2(row2),'.','MarkerSize',met MarkerSize,'LineWidth',4,'MarkerEdgeColor',[0.5 0 0]); 
%     WeirdFlo = plot(Longitude2(17),Latitude2(17),'.','MarkerSize',met MarkerSize,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);
% end
% 
% Cruise2 = line(Longitude2,Latitude2,'LineWidth',2,'Color',[0 0 0.5]);

%==================================================================

coast = load('coast'); %loads coast
geoshow ('landareas.shp','FaceColor','black'); %creates continent on map

h = colorbar; %creates colorbar
set(get(h,'label'),'string','[m]','FontSize',16) %labels colorbar

%legend([Cruise Flo WeirdFlo],'Track 12700','Float Locations 12700','Profile 17 of 12700','location','northwest') %MUST CHANGE PER CRUISE

hold off


end