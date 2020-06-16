function [dataStruct] = getSwarmData(swarmPath,nameSWARM,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd)

tStart  = datetime(yearStart,monthStart,dayStart);
tEnd    = datetime(yearEnd,monthEnd,dayEnd);

longitudes = [];
latitudes = [];
altitudes = [];
densities = [];
densitiesOrbAvg = [];
localtimes = [];
dates     = [];
for t=tStart:tEnd
    yearStr = num2str(t.Year,'%04d');
    monthStr = num2str(t.Month,'%02d');
    dayStr = num2str(t.Day,'%02d');
    dateStr = [yearStr, monthStr, dayStr];
    
    % Construct filename
    % Example: SW_OPER_DNSAPOD_2__20200102T000000_20200102T235930_0201
    filename1 = 'SW_OPER_DNS';
    filename2 = 'POD_2__';
    filename3 = 'T000000_';
    filename4 = 'T235930_0201';
    filename = [filename1 nameSWARM filename2 dateStr filename3 dateStr filename4];
    
    filePath = fullfile(swarmPath,filename,[filename '.cdf']);
    try
        % SWARMdensityData
        data = cdfread(filePath);
        
        density = cell2mat(data(:,2)); % Actual density
        densityOrbAvg = cell2mat(data(:,3)); % Orbit-averaged density
        altitude = cell2mat(data(:,5)) / 1000; % Altitude in km
        latitude = cell2mat(data(:,6)); % Geodetic latitude
        longitude = cell2mat(data(:,7)); % Geodetic longitude
        localtime = cell2mat(data(:,8));
        
        date = zeros(length(density),1);
        for i=1:length(density)
            date(i) = todatenum(data{i,1}); % Datenum
        end
        
        validity = ~cell2mat(data(:,4));
        
        longitudes  = [longitudes; longitude(validity)];
        latitudes   = [latitudes; latitude(validity)];
        altitudes   = [altitudes; altitude(validity)];
        densities   = [densities; density(validity)];
        densitiesOrbAvg = [densitiesOrbAvg; densityOrbAvg(validity)];
        localtimes  = [localtimes; localtime(validity)];
        dates = [dates; date(validity)];
%         years       = [years; year];
%         doys        = [doys; doy_];
%         seconds     = [seconds; second];
    catch
        disp(filePath);
        disp(' not found');
    end
end

dataStruct = struct;
dataStruct.longitudes = longitudes;
dataStruct.latitudes = latitudes;
dataStruct.altitudes = altitudes;
dataStruct.densities = densities;
dataStruct.densities = densitiesOrbAvg;
dataStruct.localtimes = localtimes;
dataStruct.dates = dates;

end