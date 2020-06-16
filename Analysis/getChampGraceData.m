function [dataStruct] = getChampGraceData(dataPath,filenameStart,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd)
% Example:   [xx_GPS_J2000,et,jdate] = getGPSdata('0E0E',2018,1,1,2018,1,5);

tStart  = datetime(yearStart,monthStart,dayStart);
tEnd    = datetime(yearEnd,monthEnd,dayEnd);

longitudes = [];
latitudes = [];
altitudes = [];
densities = [];
localtimes = [];
years     = [];
doys      = [];
seconds   = [];
for t=tStart:tEnd
    yearStr = num2str(t.Year,'%04d');
    doy = dayofyear(t.Year,t.Month,t.Day);
    doyStr = sprintf('%03d',doy);

    filename = [filenameStart yearStr(3:4) '_' doyStr '.mat'];
    
    filePath = fullfile(dataPath,yearStr,filename);
    try
    data = load(filePath);
    
    longitude = data.data.Lon.data;
    latitude = data.data.Lat.data;
    altitude = data.data.Height.data;
    density = data.data.Density.data;
    localtime = data.data.LocTim.data;
    year = data.data.Year.data*ones(length(longitude),1);
    doy_ = data.data.Doy.data*ones(length(longitude),1);
    second = data.data.Sec.data;
    
    longitudes  = [longitudes; longitude];
    latitudes   = [latitudes; latitude];
    altitudes   = [altitudes; altitude];
    densities   = [densities; density];
    localtimes  = [localtimes; localtime];
    years       = [years; year];
    doys        = [doys; doy_];
    seconds     = [seconds; second];
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
dataStruct.localtimes = localtimes;
dataStruct.years = years;
dataStruct.doys = doys;
dataStruct.seconds = seconds;

end