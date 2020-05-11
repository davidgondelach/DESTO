function [ephemeris] = getEphemerisForObject(path,objectID)
% Read Leolabs ephemeris data

global GM  % Earth gravitational parameter according to accurate gravity model [m^3 s^-2]

objectIDstr = num2str(objectID);
try % Try to read a ready Struct from the hard drive to avoid time-consmuing parsing.
    matFile = fullfile(path,strcat(objectIDstr,'.mat'));
    ephemeris = load( matFile ); % We already have a Struct saved there.
    
catch % This file doesn't exist yet, parse it the hard way.
    
    % Read ephemeris file
    fname = fullfile(path,['prop_ephem_' objectIDstr '.json']);
    ephemerisDataTxt = fileread(fname);
    endIndicesTxt = strfind(ephemerisDataTxt,'}{');
    startIndicesTxt = [1, endIndicesTxt+1];
    endIndicesTxt = [endIndicesTxt, length(ephemerisDataTxt)];
    
    
    startIndices = 1;
    ephemeris = struct('timestamp',{},'position',{},'velocity',{},'covariance',{});
    for i=1:length(startIndicesTxt)
        ephemerisData = jsondecode(ephemerisDataTxt(startIndicesTxt(i):endIndicesTxt(i)));
        ephemeris = [ephemeris; ephemerisData.propagation];
        startIndices(i+1) = startIndices(i) + length(ephemerisData.propagation);
    end
    
    j=1;
    for i=1:length(ephemeris)
        et = cspice_str2et( ephemeris(i).timestamp(1:end-1) );
        ephemeris(i).timestampET = et;
        if i>=startIndices(j+1)
            j = j+1;
        end
        ephemeris(i).block = j;
        ephemeris(i).coe = pv2po(ephemeris(i).position,ephemeris(i).velocity,GM);
        ephemeris(i).mee = pv2ep(ephemeris(i).position,ephemeris(i).velocity,GM);
    end
    
    save(matFile,'ephemeris'); % Save this file for later.
end

end

