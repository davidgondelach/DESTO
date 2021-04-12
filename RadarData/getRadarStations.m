function [radarStations] = getRadarStations(filepath)
%getRadarStations - Read radar stations data from json file
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

% Read station info
fname = fullfile(filepath,'instruments.json');
stationData = jsondecode(fileread(fname));
radarStations = stationData.instruments;

end

%------------- END OF CODE --------------