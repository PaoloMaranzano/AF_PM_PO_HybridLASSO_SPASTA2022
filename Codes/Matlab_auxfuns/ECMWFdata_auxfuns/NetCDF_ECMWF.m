function OutputDataset = NetCDF_ECMWF(Varargin,debug)

%% Open stations registry
stats = readtable('Data/ARPA_data/ARPA_stats_registry.csv');

%% put 1 for debugging or put 0 for function usage
if 0
    debug = 1
end

if debug == 1
    v = 1;
    vars = {'Temperature'};
    r = find(~ismember(vars, {'RelHumid','Humidex'}));
    vars = [vars(r)];
else
    % Exclude from the list 'RelHumid' and 'Humidex' because they
    % are computed combining 'Temperature' and 'Dewpoint'.
    % They do not have a NetCDF.
    vars = Varargin;
    r = find(~ismember(vars, {'RelHumid','Humidex'}));
    vars = [vars(r)];
end
output = cell(length(vars),1);

for v = 1:length(vars)
    
    %% Open NetCDF
    var = vars{v};
    switch var
        case 'Temperature'
            ncid_name = 'Temp2m_Lom_20172020.nc';
            vname = 't2m';
        case 'Dewpoint'
            ncid_name = 'Temp2m_Lom_20172020.nc';
            vname = 'd2m';
        case 'Rainfall'
            ncid_name = 'SoilPrec_Lom_20172020.nc';
            vname = 'tp';
        case 'SoilTemp'
            ncid_name = 'SoilPrec_Lom_20172020.nc';
            vname = 'stl1';
        case 'Pressure'
            ncid_name = 'PressRadiation_Lom_20172020.nc';
            vname = 'sp';
        case 'WindU'
            ncid_name = 'Wind_Lom_20172020.nc';
            vname = 'u10';
        case 'WindV'
            ncid_name = 'Wind_Lom_20172020.nc';
            vname = 'v10';
        case 'VegetationLow'
            ncid_name = 'Vegetation_Lom_20172020.nc';
            vname = 'lai_lv';
        case 'VegetationHigh'
            ncid_name = 'Vegetation_Lom_20172020.nc';
            vname = 'lai_hv';
        case 'GeopotHeight'
            ncid_name = 'Geopotential.nc';
            vname = 'z';
        case 'SolarRadiation'
            ncid_name = 'PressRadiation_Lom_20172020.nc';
            vname = 'ssr';
        otherwise
            disp('other value')
    end

    ncid = netcdf.open(ncid_name);
    varid = netcdf.inqVarID(ncid,vname);
    
    %% Extract information from NetCDf file
    % info = ncinfo(ncid_name);
    % Dimensions
    % info.Dimensions.Length
    % Dimensions names
    % info.Dimensions.Name
    % Attributes
    % a = info.Variables.Attributes;
    % info.Variables.Name
    % Variables
    % b = info.Variables;
    
    %% Extract information from NetCDf file
    % Coordinates (longitude and latitude
    longid = netcdf.inqDimID(ncid,'longitude');
    long = netcdf.getVar(ncid,longid);
    latid = netcdf.inqDimID(ncid,'latitude');
    lat = netcdf.getVar(ncid,latid);
    % Time
    timeid = netcdf.inqDimID(ncid,'time');
    date = netcdf.getVar(ncid,timeid);
    t = double(date)/24 + datenum('1900-01-01 00:00:00');
    t2 = datetime(datestr(t,'dd-mm-yyyy HH:MM:SS'),'InputFormat','dd-MM-yyyy HH:mm:ss');
    t2 = datetime(t2,'InputFormat','yyyy-MM-dd HH:mm:ss');
    t2 = array2table(t2);
    t2.Properties.VariableNames = {'DateTime'};
    
    %% Data extraction
    dataset = netcdf.getVar(ncid,varid);
    dataset = double(dataset);
    % (34321/24)/365  %%%% ~ 4 years
    
    %% Scaling and missing data management
    [noFillMode,fillValue] = netcdf.inqVarFill(ncid,varid);
    idx_nan = dataset == fillValue;
    dataset(idx_nan) = nan;
    
    % Geopotential height does not have offset and scaling factor
    if strcmp(var,'GeopotHeight') == 0
        add_offset = netcdf.getAtt(ncid,varid,'add_offset');
        scale_factor = netcdf.getAtt(ncid,varid,'scale_factor');
    end
    
    % Measure units
    switch var
        case {'Temperature','Dewpoint','SoilTemp'}
            % From Kelvin to Celsius
            dataset = dataset.*scale_factor + add_offset - 273.15;
        case 'Pressure'
            % From Pascal to HettoPascal
            dataset = dataset.*scale_factor + add_offset;
            dataset = dataset/100;
        case 'Rainfall'
            % From metres to millimetres
            dataset = dataset.*scale_factor + add_offset;
            dataset = dataset.*1000;
        case 'GeopotHeight'
            % From geopotenzial to geopotential height 
            % Earth's gravitational acceleration, g (=9.80665 m s-2)
            dataset = dataset./9.80665;
        otherwise
            dataset = dataset.*scale_factor + add_offset;
    end
    
%     plot(squeeze(dataset(25,25,:)))
%     hist(squeeze(dataset(25,25,:)))
%     sum(squeeze(dataset(25,25,:)) == 0)
    
    %% Matching stations and gridded data
    y = zeros(length(date),size(stats,1));
    for i = 1:size(stats,1)
        coords = [stats(i,:).Latitude,stats(i,:).Longitude];
        % lat_south = lat(min(find(lat<=coords(1))));
        lat_north = lat(max(find(lat>=coords(1))));
        long_west = long(max(find(long<=coords(2))));
        % long_east = long(min(find(long>=coords(2))));
        NW = [find(long == long_west) , find(lat == lat_north)];
        y(:,i) = squeeze(dataset(NW(1),NW(2),:));
    end
    
    % coords_mins = [nan , nan];
    % coords_mins(1) = degreesDecimalToMinutes(coords(1))
    % degreesMinutesToDecimal([45,6,0])
    % degreesMinutesToDecimal([9,7,0])
    
    %% Stacking data as table
    y = array2table(y);
    y.Properties.VariableNames = stats.New_cod_stz;
    y = [t2 , y];
    y3 = stack(y,2:size(y,2));
    
    %% Variable name
    switch var
        case 'Temperature'
            y3.Properties.VariableNames = {'Date','IDStat','Temperature'};
        case 'Dewpoint'
            y3.Properties.VariableNames = {'Date','IDStat','Dewpoint'};
        case 'SoilTemp'
            y3.Properties.VariableNames = {'Date','IDStat','SoilTemp'};
        case 'Rainfall'
            y3.Properties.VariableNames = {'Date','IDStat','Rainfall'};
        case 'Pressure'
            y3.Properties.VariableNames = {'Date','IDStat','Pressure'};
        case 'WindU'
            y3.Properties.VariableNames = {'Date','IDStat','WindU'};
        case 'WindV'
            y3.Properties.VariableNames = {'Date','IDStat','WindV'};
        case 'VegetationHigh'
            y3.Properties.VariableNames = {'Date','IDStat','VegetationHigh'};
        case 'VegetationLow'
            y3.Properties.VariableNames = {'Date','IDStat','VegetationLow'};
        case 'GeopotHeight'
            y3.Properties.VariableNames = {'Date','IDStat','GeopotHeight'};
        case 'SolarRadiation'
            y3.Properties.VariableNames = {'Date','IDStat','SolarRadiation'};
        otherwise
            y3.Properties.VariableNames = {'Date','IDStat','Other'};
    end
    
    %% Missing data replacement
    nan_vals = isnan(y3{:,3});
    y3{nan_vals,3} = nanmean(y3{:,3});
    
    %% Output for the single variable
    output{v} = y3;
    netcdf.close(ncid);
    
end

%% Checking for the presence of Geopotenzial NetCDF
NonGeopot_idx = find(~strcmp(vars,'GeopotHeight'));
Geopot_idx = find(strcmp(vars,'GeopotHeight'));

%% Datasets merge
if sum(strcmp(vars,'GeopotHeight')) == 0
    T = MultiJoin({output{NonGeopot_idx}},{'Date','IDStat'});
else
    T = MultiJoin({output{NonGeopot_idx}},{'Date','IDStat'});
    output{Geopot_idx}.Date = [];
    T = MultiJoin({T,output{Geopot_idx}},{'IDStat'});
end

%% Relative humidity and Humidex
if sum(strcmp(T.Properties.VariableNames,{'Temperature'})) + ...
        sum(strcmp(T.Properties.VariableNames,{'Dewpoint'})) == 2
    T.RelHumid = 100*exp(17.625*T.Dewpoint./(243.04+T.Dewpoint))./exp(17.625*T.Temperature./(243.04+T.Temperature));
    T.Humidex = T.Temperature + 5/9*(6.112*10.^(7.5*T.Temperature./(237.7+T.Temperature)).*T.RelHumid./100-10);
end

%% Dataset filtering according to the date
d_start = datetime(2017,01,01,0,0,0);
d_end = datetime(2020,11,30,23,0,0);
idx_date = T.Date <= d_end & T.Date >= d_start;

%% Output dataset
OutputDataset = T(idx_date,:);

end