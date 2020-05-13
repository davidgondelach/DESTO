function [radarStations] = getRadarStations(path)

% Read station info
fname = fullfile(path,'instruments.json');
stationData = jsondecode(fileread(fname));
radarStations = stationData.instruments;

end