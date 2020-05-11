function [ephemerisObjects] = getEphemeris(path,objectIDs)

ephemerisObjects = struct('noradID',{},'ephemeris',{});
for i=1:length(objectIDs)
    objectID = objectIDs(i);
    [ephemeris] = getEphemerisForObject(path,objectID);
    
    ephemerisObjects(i).noradID = objectID;
    ephemerisObjects(i).ephemeris = ephemeris.ephemeris;
end
