%% get_aviso info
    function get_Matt_alt(cruise,num)  %sets as function

    meta = MetaFile(cruise);    %locates info in metainfo file 

    Alt = load('nrt_Matt_Alt.mat');
    %outfile = ['Alty_matt_',cruise,'_',num,'.mat'];
    
    allLon1 = Alt.xsave(723:1440)-360; %this loads the first collumn (the longitude) and adjusts from 180:360 to be 0:180
    allLon2 = Alt.xsave(1:722)-360; %adjusts from 1:180 to -180:0
    uLon2 = unique(allLon1); %finds the unique longitude values each time - the values
    uLon1 = unique(allLon2);

    Alt.xsave = [uLon1; uLon2];

    %create new indexes for the correct cruise lat/lon to be called later 
    LatInd = find(Alt.ysave <= meta.LatMax & Alt.ysave >= meta.LatMin); %finds lat indices within max/min bounds
    LonInd = find(Alt.xsave <= meta.LonMax & Alt.xsave >= meta.LonMin); %finds lon indices within max/min bounds

    %find the curl data for those correct indexes
    lat = Alt.ysave(LatInd);
    lon = Alt.xsave(LonInd);
    
    % vidObj = VideoWriter('SSH.avi');
    vidObj = VideoWriter('P17_Aviso.mp4','MPEG-4');	% on linux, MATLAB can only make .avi movies
    vidObj.Quality = 100;				% I think this should work on a mac
    vidObj.FrameRate = 6;				% this sets the frame rate
    open(vidObj);
    writeVideo(vidObj, getframe(gcf));
    
    %depth = squeeze(Alt.dsave(348,LatInd,LonInd)); %for 12778
    for i = 340:470
        f = load('12700.mat');
        depth = squeeze(Alt.dsave(i,LatInd,LonInd)); %for 12700
        %depppp = Alt.dsave;
        %avdepth = squeeze(mean(Alt.dsave(:,LatInd,LonInd)));
        %anom = depth-avdepth;
        time = datestr(Alt.tsave(i,:,:));
        %stddev = squeeze(std(Alt.dsave(:,LatInd,LonInd)));
        
        %outfile = ['Alty_matt_',num2str(i),'.mat'];
        %save(outfile)
        
        hold on
        
        title(['Aviso: ',num2str(time)],'FontSize',20)%,titl],'FontSize',16) %creates title
        ylabel('Latitude (degrees)','FontSize',16)    %add y label for Latitude
        xlabel('Longitude (degrees)','FontSize',16)   %add x label for Longitude
        imagesc(lon,lat,(depth)) % creates the basic image we want to edit
         set(gca,'YLim',[meta.LatMin meta.LatMax]); % sets map y-limits to the Latitudes specified for your area
         set(gca,'XLim',[meta.LonMin meta.LonMax]); % sets map x-limits to the Latitudes specified for your area
        set(gca,'FontSize',18)
        
        coast = load('coast'); %loads coast
        geoshow ('landareas.shp','FaceColor','black'); %creates continent on map

        h = colorbar; %creates colorbar
        set(get(h,'label'),'string','[m]','FontSize',16) %labels colorbar
        
        contour(lon,lat,(depth),[0 0],'Color','k')                  % outline the coasts with black contour by accessing Bathymetry file lon/lat/dep
        [L,m] = contour(lon,lat,(depth),[-0.4:0.1:0.8],'Color','k'); % -2000m isobath (farthest floats go to) black countour and assign lon/lat positions to L/M
        clabel(L,m,'LabelSpacing',200); %labels the contours
        
        plot(-41.875,-33.875,'r.','markersize',20)
        
        %saveas( gcf, ['Aviso_',num2str(i)], 'jpg' );
        hold off
        
        drawnow()				% these two commands save the current plot as the next
        open(vidObj);
        writeVideo(vidObj, getframe(gcf));	% frame in the movie
        
    end 
        
%     %save as an outfile
%     save(outfile)

    end