function [stationsPosECEF] = getStationPositionECEF(radarStations)
%getStationPositionECEF - Compute ground station positions in
%Earth-centered Earth-fixed frame (WGS84) from longitude, latitude and
%altitude.
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

% Stations' longitude, latitude and altitude
stationsLonLatAlt = [   deg2rad([radarStations.longitude]); 
                        deg2rad([radarStations.latitude]); 
                        [radarStations.altitude] / 1e3];

% Convert to ECEF postion
re = 6378.137; % WGS84: a (equatorial radius): 6378137.0 m
f = 1 / 298.257223563; % WGS84: 1/f (inverse flattening): 298.257223563
stationsPosECEF = cspice_georec( stationsLonLatAlt(1,:), stationsLonLatAlt(2,:), stationsLonLatAlt(3,:), re, f);

end

%------------- END OF CODE --------------