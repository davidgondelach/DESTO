function [gpsData] = readGPSdataFromJsonFile(filepath)
%readGPSmessage - Read GPS message and convert state to J2000 reference
%frame and time to ephemeris time.
%   WGS84 is assumed to be equal to the International Terrestrial Reference Frame (ITRF) with an accuracy of a few centimeters (Montenbruck and Gill)


% Read measurements
gpsData = jsondecode(fileread(filepath));

emptyStates = zeros(1,length(gpsData));
for i=1:length(gpsData)
    et = cspice_unitim( gpsData(i).gnss_sample_t_j2000, 'TDT', 'ET' );
    gpsData(i).tET = et;
    
    xx_ECEF(1,1) = gpsData(i).gnss_sample_rx_ecef;
    xx_ECEF(2,1) = gpsData(i).gnss_sample_ry_ecef;
    xx_ECEF(3,1) = gpsData(i).gnss_sample_rz_ecef;
    xx_ECEF(4,1) = gpsData(i).gnss_sample_vx_ecef;
    xx_ECEF(5,1) = gpsData(i).gnss_sample_vy_ecef;
    xx_ECEF(6,1) = gpsData(i).gnss_sample_vz_ecef;

%     % Convert from ECEF to J2000 and from meters to kilometers
%     % using Vallado code
%     jdatestr    = cspice_et2utc( et, 'J', 12 );
%     jdate       = str2double(jdatestr(4:end)); % Cut leading 'JD ' off from string
%     [rr_J2000, vv_J2000] = convertECEFtoJ2000(xx_ECEF(1:3), xx_ECEF(4:6), jdate);
%     gpsData(i).xx_J2000 = [rr_J2000; vv_J2000] / 1000;
    
    % Convert from ECEF to J2000 using SPICE
    % and convert from meters to kilometers
    xform = cspice_sxform('ITRF93', 'J2000', et );
    gpsData(i).xx_J2000 = xform * xx_ECEF / 1000;
    
    if norm(xx_ECEF(1:3)) == 0
        emptyStates(i) = 1;
    end
    
end
gpsData = gpsData(~emptyStates);

end