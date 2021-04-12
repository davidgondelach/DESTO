function [gpsData,gpsfilepath] = getGPSdataFromJsonFile(GPSdataPath,objectPlanetIDstr,objectNORADID,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd)
% Example:   [xx_GPS_J2000,et,jdate] = getGPSdata('0E0E',2018,1,1,2018,1,5);
%
% Copyright (C) 2021 by David Gondelach
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Aug 2020; Last revision: 31-Aug-2020

%------------- BEGIN CODE --------------

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

%------------- END OF CODE --------------