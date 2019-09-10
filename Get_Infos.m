function Get_Infos(cruise) 

get_aviso(cruise)
get_area_from_etopo2(cruise) 
% get_buoyancyflux(cruise)
% get_chloro1(cruise)
% get_chloro2(cruise)
% get_curl(cruise)
% get_curl_vectors(cruise)
% get_fresh_water(cruise)
% get_heat(cruise)
% get_sst(cruise)

%% get_aviso info
    function get_aviso(cruise)  %sets as function

    Meta = MetaFile(cruise);    %locates info in MetaFile file 

    Alt = load('Aviso_Data.mat');
    outfile = ['Aviso_',cruise,'.mat'];

    %create new indexes for the correct cruise lat/lon to be called later 
    LatInd = find(Alt.Aviso.lat <= Meta.LatMax & Alt.Aviso.lat >= Meta.LatMin); %finds lat indices within max/min bounds
    LonInd = find(Alt.Aviso.lon <= Meta.LonMax & Alt.Aviso.lon >= Meta.LonMin); %finds lon indices within max/min bounds

    %find the curl data for those correct indexes
    lat = Alt.Aviso.lat(LatInd);
    lon = Alt.Aviso.lon(LonInd);
    sla = Alt.Aviso.sla(LatInd,LonInd);
    adt = Alt.Aviso.adt(LatInd,LonInd);

    %save as an outfile
    save(outfile)

    end
%% bathymetry info
    function get_area_from_etopo2(cruise) 

    Meta = MetaFile(cruise);

    outfile = ['ETOPO2v2g_f4_',cruise,'.mat']; %outfile we want to create

    % set variables needed for cruise of interest
    Meta = MetaFile(cruise);


    % etopo_loc = '/dat2/etopo/etopo2';  % location of etopo file % now in MetaFile.m

    F_etopo   = fullfile(Meta.etopo_loc,'ETOPO2v2g_f4.nc'); % name of etopo1 file

    % Read the variables in the netcdf (.nc) etopo file
    varnames = {};
    ncid = netcdf.open(F_etopo,'NC_NOWRITE');
    [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
    for i = 1:numvars
       [varnames{i}, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,i-1);
       Data.([varnames{i}]) = netcdf.getVar(ncid,i-1);
    end
    netcdf.close(ncid);


    % finds longitude within the right area for cruise
    ix =  find(Data.x <= Meta.LonMax  & Data.x >= Meta.LonMin);
    Bathy.x = Data.x(ix);

    % finds latitude within the right area
    iy =  find(Data.y <= Meta.LatMax  & Data.y >= Meta.LatMin);
    Bathy.y = Data.y(iy);

    %reads the bathy data and links it to the lat/lon indexes
    Bathy.z = Data.z(ix,iy);

    % dont do this anymore
    %etopo_P15S.z(find(etopo_P15S.z >= 0)) = 0; % make all bathy>0 be 0


    save(outfile,'Bathy')

    end
%% get_buoyancy
    function get_buoyancyflux(cruise);

    Meta = MetaFile(cruise);

    Buoyancy = load('BuoyancyData.mat'); %loads HeatFlux file 
    outfile = ['BuoyancyFlux_',cruise,'.mat'];

    %create new indexes for the correct cruise lat/lon to be called later 
    LatInd = find(Buoyancy.Buoyancy.lat <= Meta.LatMax & Buoyancy.Buoyancy.lat >= Meta.LatMin);
    LonInd = find(Buoyancy.Buoyancy.lon <= Meta.LonMax & Buoyancy.Buoyancy.lon >= Meta.LonMin);

    %find the curl data for those correct indexes
    BuoyancyLat = Buoyancy.Buoyancy.lat(LatInd);
    BuoyancyLon = Buoyancy.Buoyancy.lon(LonInd);
    BuoyancyImage = Buoyancy.Buoyancy.data(LatInd,LonInd);

    %save as an outfile
    save(outfile)

    end
%% 
    function get_chloro1(cruise)

    Meta = MetaFile(cruise);

    outfile = ['Chloro_data_' , cruise ,'.mat'];

    Chl = load('Chloro_data_.mat');

    ix = find(Chl.Chlor.lon >= Meta.LonMin & Chl.Chlor.lon <= Meta.LonMax);
    Chloro.Lon = Chl.Chlor.lon(ix); 

    iy = find(Chl.Chlor.lat >= Meta.LatMin & Chl.Chlor.lat <= Meta.LatMax);
    Chloro.Lat = Chl.Chlor.lat(iy); 

    Chloro.chlor = Chl.Chlor.Chlor(ix,iy);
    Chloro.chlor = Chloro.chlor(Chloro.chlor <= 1);

    save(outfile)
    
    end

%%
    function get_chloro2(cruise)

    Meta = MetaFile(cruise);

    outfile = ['Chloro2_data_' , cruise ,'.mat'];

    Chl = load('Chloro_data_.mat');

    ix = find(Chl.Chlor.lon >= Meta.LonMin & Chl.Chlor.lon <= Meta.LonMax);
    Chloro.Lon = Chl.Chlor.lon(ix); 

    iy = find(Chl.Chlor.lat >= Meta.LatMin & Chl.Chlor.lat <= Meta.LatMax);
    Chloro.Lat = Chl.Chlor.lat(iy); 

    Chloro.chlor = Chl.Chlor.Chlor(ix,iy);
    Chloro.chlor = Chloro.chlor(Chloro.chlor <= 2);

    save(outfile)
    
        
    end
%% 
    function get_curl(cruise)  %sets up code as a function that can be called 

    curl = load('Curl.mat');  %loads the curl file from common files
    Meta = MetaFile(cruise); %links us to the Meta info so we can call what we need to call
    outfile = ['CurlData_',cruise,'.mat']; %creates name of outfile that will be created by the data

    %create new indexes for the correct cruise lat/lon to be called later 
    newLat = find(curl.Curl.lat <= Meta.LatMax & curl.Curl.lat >= Meta.LatMin);
    newLon = find(curl.Curl.lon <= Meta.LonMax & curl.Curl.lon >= Meta.LonMin);
    %below is an attempt at isolating and changing color of NaNs that hasn't
    %been successful yet
    newNaN = find(isnan(curl.Curl.curl)); 

    %find the curl data for those correct indexes
    CurlLat = curl.Curl.lat(newLat);
    CurlLon = curl.Curl.lon(newLon);
    CurlImage = curl.Curl.curl(newLat,newLon); %calls indices of lat and lon to match them up
    %below is an attempt at isolating and changing color of NaNs that hasn't
    %been successful yet
    CurlNaN = curl.Curl.curl(newNaN);

    %save as an outfile
    save(outfile)
    
    end
%% 
    function get_curl_vectors(cruise)  %sets up code as a function that can be called 

    vector = load('CurlVector.mat');  %loads the information code we need
    Meta = MetaFile(cruise); %links us to the Meta info so we can call what we need to call
    outfile = ['VectorCurlData_',cruise,'.mat'];

    %create new indexes for the correct cruise lat/lon to be called later 
    newLat = find(vector.Vector.lat <= Meta.LatMax & vector.Vector.lat >= Meta.LatMin);
    newLon = find(vector.Vector.lon <= Meta.LonMax & vector.Vector.lon >= Meta.LonMin);

    %find the curl data for those correct indexes
    VLat = vector.Vector.lat(newLat);
    VLon = vector.Vector.lon(newLon);
    VAngle = vector.Vector.Ang(newLat,newLon)*0.0174533;
    VSpeed = vector.Vector.Speed(newLat,newLon);

    Lat2 = (VSpeed.*sin(VAngle)); %creates the x-component of the vector we want to map
    Lon2 = (VSpeed.*cos(VAngle)); %creates the y-component of the vector we want to map

    %save as an outfile
    save(outfile)
    
    end
%% 
    function get_fresh_water(cruise);

    Meta = MetaFile(cruise);

    fwater = load('FreshWater.mat'); %loads FreshWaterFlux file 
    outfile = ['FreshWater_',cruise,'.mat'];

    %create new indexes for the correct cruise lat/lon to be called later 
    LatInd = find(fwater.FreshWater.lat <= Meta.LatMax & fwater.FreshWater.lat >= Meta.LatMin);
    LonInd = find(fwater.FreshWater.lon <= Meta.LonMax & fwater.FreshWater.lon >= Meta.LonMin);

    %find the curl data for those correct indexes
    FreshWaterLat = fwater.FreshWater.lat(LatInd);
    FreshWaterLon = fwater.FreshWater.lon(LonInd);
    FreshWaterImage = fwater.FreshWater.water(LatInd,LonInd);

    %save as an outfile
    save(outfile)

    end
%% 
    function get_heat(cruise);

    Meta = MetaFile(cruise);

    HeatFlux = load('HeatFlux.mat'); %loads HeatFlux file 
    outfile = ['HeatFlux_',cruise,'.mat'];

    %create new indexes for the correct cruise lat/lon to be called later 
    LatInd = find(HeatFlux.Heat.lat <= Meta.LatMax & HeatFlux.Heat.lat >= Meta.LatMin);
    LonInd = find(HeatFlux.Heat.lon <= Meta.LonMax & HeatFlux.Heat.lon >= Meta.LonMin);

    %find the curl data for those correct indexes
    HeatLat = HeatFlux.Heat.lat(LatInd);
    HeatLon = HeatFlux.Heat.lon(LonInd);
    HeatImage = HeatFlux.Heat.heat(LatInd,LonInd);

    %save as an outfile
    save(outfile)
    
    end
%% 
    function get_sst(cruise);

    Meta = MetaFile(cruise);

    SSTFlux = load('SST_Data.mat'); %loads SSTFlux file 
    outfile = ['SSTFlux_',cruise,'.mat'];

    %create new indexes for the correct cruise lat/lon to be called later 
    LatInd = find(SSTFlux.SST.lat <= Meta.LatMax & SSTFlux.SST.lat >= Meta.LatMin);
    LonInd = find(SSTFlux.SST.lon <= Meta.LonMax & SSTFlux.SST.lon >= Meta.LonMin);

    %find the curl data for those correct indexes
    SSTLat = SSTFlux.SST.lat(LatInd);
    SSTLon = SSTFlux.SST.lon(LonInd);
    SSTImage = SSTFlux.SST.sst(LonInd,LatInd,438);


    %save as an outfile
    save(outfile)
    
    end

end