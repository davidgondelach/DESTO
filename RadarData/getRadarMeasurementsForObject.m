function [measurements] = getRadarMeasurementsForObject(path,objectID)
% Read Leolabs radar measurement data

objectIDstr = num2str(objectID);
try % Try to read a ready Struct from the hard drive to avoid time-consuming parsing.
    matFile = fullfile(path,strcat(objectIDstr,'_measurements.mat'));
    measurements = load( matFile ); % We already have a Struct saved there.
    
catch % This file doesn't exist yet, parse it.
    
    % Read measurements
    fname = fullfile(path,['measurements-' objectIDstr '.json']);
    measurementData = jsondecode(fileread(fname));
    measurements = measurementData.measurements;
    for i=1:length(measurements)
        et = cspice_str2et( measurements(i).measuredAt(1:end-1) );
        measurements(i).measuredAtET = et;
    end
    
    save(matFile,'measurements'); % Save this file for later.
end

end

