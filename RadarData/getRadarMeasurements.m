function [radarMeasurementsObjects] = getRadarMeasurements(path,objectIDs)
%getRadarMeasurements - Load Leolabs radar measurement data for objects
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

radarMeasurementsObjects = struct('noradID',{},'measurements',{});
for i=1:length(objectIDs)
    objectID = objectIDs(i);
    [measurements] = getRadarMeasurementsForObject(path,objectID);
    
    radarMeasurementsObjects(i).noradID = objectID;
    radarMeasurementsObjects(i).measurements = measurements.measurements;
end

%------------- END OF CODE --------------