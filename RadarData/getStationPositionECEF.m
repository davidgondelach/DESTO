function [stationsPosECEF] = getStationPositionECEF(radarStations)

% Stations' longitude, latitude and altitude
stationsLonLatAlt = [   deg2rad([radarStations.longitude]); 
                        deg2rad([radarStations.latitude]); 
                        [radarStations.altitude] / 1e3];

% Convert to ECEF postion
% https://confluence.qps.nl/qinsy/9.1/en/world-geodetic-system-1984-wgs84-182618391.html
re = 6378.137; % WGS84: a (equatorial radius): 6378137.0 m
f = 1 / 298.257223563; % WGS84: 1/f (inverse flattening): 298.257223563
stationsPosECEF = cspice_georec( stationsLonLatAlt(1,:), stationsLonLatAlt(2,:), stationsLonLatAlt(3,:), re, f);

end

