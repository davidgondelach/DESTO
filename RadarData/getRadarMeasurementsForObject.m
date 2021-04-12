function [measurements] = getRadarMeasurementsForObject(path,objectID)
%getRadarMeasurementsForObject - Read Leolabs radar measurement data from 
% json file
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
        % Compute ephemeris time
        et = cspice_str2et( measurements(i).measuredAt(1:end-1) );
        measurements(i).measuredAtET = et;
    end
    
    save(matFile,'measurements'); % Save this file for later.
end

end

%------------- END OF CODE --------------