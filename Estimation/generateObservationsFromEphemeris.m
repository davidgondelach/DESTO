function [meas,RM] = generateObservationsFromEphemeris(ephemerisObjects,et0,etf)
% Compute measurements and measurement covariance

nofMinBetweenObs = 10;
% n = 6; % Measurement state length for single object

meas = [];
RM = [];
for i=1:length(ephemerisObjects)
    nofMeas = size(meas,2);
    
    epochs = [ephemerisObjects(i).ephemeris(1:nofMinBetweenObs:end).timestampET];
    nofObs = length(epochs);
    
    % Epochs
    meas(1,nofMeas+1:nofMeas+nofObs) = epochs;
    % Object number
    meas(2,nofMeas+1:nofMeas+nofObs) = i;
    % Position and velocity
    meas(3:5,nofMeas+1:nofMeas+nofObs) = [ephemerisObjects(i).ephemeris(1:nofMinBetweenObs:end).position];
    meas(6:8,nofMeas+1:nofMeas+nofObs) = [ephemerisObjects(i).ephemeris(1:nofMinBetweenObs:end).velocity];
    
    % Covariance
    RM(1:6,1:6,nofMeas+1:nofMeas+nofObs) = reshape([ephemerisObjects(i).ephemeris(1:nofMinBetweenObs:end).covariance],6,6,nofObs);
end

% Remove measurements outside time window of interest
withinWindow = find(meas(1,:)>=et0 & meas(1,:)<=etf);
meas = meas(:,withinWindow);
RM = RM(:,:,withinWindow);

% Sort measurements on epoch
[~, order] = sort(meas(1,:));
meas = meas(:,order);
RM = RM(:,:,order);

% Scale from meters to kilometers
meas(3:end,:) = meas(3:end,:) / 1e3;
RM = RM / 1e6;

end

