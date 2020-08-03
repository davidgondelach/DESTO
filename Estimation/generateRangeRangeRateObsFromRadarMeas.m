function [meas,RM] = generateRangeRangeRateObsFromRadarMeas(radarMeasurementsObjects,radarStations,et0,etf)
% Compute measurements and measurement covariance

% Stations' ECEF postion
stationsPosECEF = getStationPositionECEF(radarStations);

meas = [];
RM = [];
for i=1:length(radarMeasurementsObjects)
    nofMeas = size(meas,2);
    
    epochs = [radarMeasurementsObjects(i).measurements(:).measuredAtET];
    nofObs = length(epochs);
    
    % Epochs
    meas(1,nofMeas+1:nofMeas+nofObs) = epochs;
    % Object number
    meas(2,nofMeas+1:nofMeas+nofObs) = i;
    
    % Tracking station
    stations = {radarMeasurementsObjects(i).measurements(:).instrument};
    for j=1:nofObs
        
        % Range in km
        meas(3,nofMeas+j) = [radarMeasurementsObjects(i).measurements(j).corrected.range] / 1e3;
        % Range rate in km/s
        meas(4,nofMeas+j) = [radarMeasurementsObjects(i).measurements(j).corrected.doppler] / 1e3;
        
        % Get station's position and velocity
        stationIndex = find(strcmp({radarStations.id},stations{j}));
        stationPos = stationsPosECEF(:,stationIndex);
        stationVel = zeros(3,1);
        stationECEF = [stationPos; stationVel];
        
        % Station's position and velocity in J2000 frame
        ecefToJ2000 = cspice_sxform('ITRF93', 'J2000', epochs(j) ); % ECEF to J2000 transformation matrix
        stationJ2000 = ecefToJ2000 * stationECEF; % State in ECEF
        meas(5:10,nofMeas+j) = stationJ2000;

        % Measurement noise
        % Range error in km
        RM(1,1,nofMeas+j) = [radarMeasurementsObjects(i).measurements(j).corrected.rangeError] / 1e3;
        % Range rate error in km/s
        RM(2,2,nofMeas+j) = [radarMeasurementsObjects(i).measurements(j).corrected.dopplerError] / 1e3;
    end

end

% Remove measurements outside time window of interest
withinWindow = find(meas(1,:)>=et0 & meas(1,:)<=etf);
meas = meas(:,withinWindow);
RM = RM(:,:,withinWindow);

% Sort measurements on epoch
[~, order] = sort(meas(1,:));
meas = meas(:,order);
RM = RM(:,:,order);

RM(1,1,:) = RM(1,1,:)*1; % Multiply 1-sigma range error by factor 1
RM(2,2,:) = RM(2,2,:)*2; % Multiply 1-sigma range-rate error by factor 2

% Measurement standard deviation to covaraince
RM = RM.^2;

end

