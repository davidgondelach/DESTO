function [gpsData] = readGPSdataFromJsonFile(filepath)
%readGPSmessage - Read GPS message and convert state to J2000 reference
%frame and time to ephemeris time.
%   WGS84 is assumed to be equal to the International Terrestrial Reference Frame (ITRF) with an accuracy of a few centimeters (Montenbruck and Gill)

try % Try to read a ready Struct from the hard drive to avoid time-consuming parsing.
    matFile = strcat(filepath(1:end-5),'.mat');
    gpsData = load( matFile, 'gpsData' ); % We already have a Struct saved there.
    gpsData = gpsData.gpsData;
    
catch % This file doesn't exist yet, parse it.
    
    % Read measurements
    gpsData = jsondecode(fileread(filepath));
    
    corruptStates = zeros(1,length(gpsData));
    for i=1:length(gpsData)
        et = cspice_unitim( gpsData(i).gnss_sample_t_j2000, 'TDT', 'ET' );
        gpsData(i).tET = et;
        
        xx_ECEF(1,1) = gpsData(i).gnss_sample_rx_ecef;
        xx_ECEF(2,1) = gpsData(i).gnss_sample_ry_ecef;
        xx_ECEF(3,1) = gpsData(i).gnss_sample_rz_ecef;
        xx_ECEF(4,1) = gpsData(i).gnss_sample_vx_ecef;
        xx_ECEF(5,1) = gpsData(i).gnss_sample_vy_ecef;
        xx_ECEF(6,1) = gpsData(i).gnss_sample_vz_ecef;
        
        % Convert from ECEF to J2000 using SPICE
        % and convert from meters to kilometers
        xform = cspice_sxform('ITRF93', 'J2000', et );
        gpsData(i).xx_J2000 = xform * xx_ECEF / 1000;
        
        if norm(xx_ECEF(1:3)) < 6000
            corruptStates(i) = 1;
        end
                
    end
    gpsData = gpsData(~corruptStates);
    
%     save(matFile,'gpsData'); % Save this file for later.
end

end