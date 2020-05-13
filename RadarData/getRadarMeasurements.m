function [radarMeasurementsObjects] = getRadarMeasurements(path,objectIDs)

radarMeasurementsObjects = struct('noradID',{},'measurements',{});
for i=1:length(objectIDs)
    objectID = objectIDs(i);
    [measurements] = getRadarMeasurementsForObject(path,objectID);
    
    radarMeasurementsObjects(i).noradID = objectID;
    radarMeasurementsObjects(i).measurements = measurements.measurements;
end
