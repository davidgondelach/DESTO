function [gpsData] = getGPSdataFromJsonFile(GPSdataPath,objectPlanetIDstr,objectNORADID,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd)
% Example:   [xx_GPS_J2000,et,jdate] = getGPSdata('0E0E',2018,1,1,2018,1,5);

% % Download GPS files
% for i=1:31
% websave(['gps_data_201801',num2str(i,'%02d'),'.zip'],['http://ephemerides.planet-labs.com/gps_data_201801',num2str(i,'%02d'),'.zip']);
% end

tStart  = datetime(yearStart,monthStart,dayStart);
tEnd    = datetime(yearEnd,monthEnd,dayEnd);

gpsData = struct('gnss_sample_gdop',{},'gnss_sample_prns_acquired',{},'gnss_sample_rx_ecef',{},'gnss_sample_ry_ecef',{},'gnss_sample_rz_ecef',{}, ...
                 'gnss_sample_t_j2000',{},'gnss_sample_vx_ecef',{},'gnss_sample_vy_ecef',{},'gnss_sample_vz_ecef',{},'sat_id',{},'t',{},'tET',{},'xx_J2000',{});

for t=tStart:tEnd
    yearStr = num2str(t.Year);
    monthStr = num2str(t.Month,'%02d');
    dayStr = num2str(t.Day,'%02d');
    dateStr = [yearStr,monthStr,dayStr];
    
    gpsfilepath = fullfile(GPSdataPath,['gps_data_' dateStr],['gps_' objectPlanetIDstr '_' num2str(objectNORADID) '_' dateStr '.json']);
    gpsDataSingleDay = readGPSdataFromJsonFile(gpsfilepath);

    gpsData = [gpsData; gpsDataSingleDay];
end

end