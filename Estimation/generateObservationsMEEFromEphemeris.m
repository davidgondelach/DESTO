function [meas,RM] = generateObservationsMEEFromEphemeris(ephemerisObjects,et0,etf,GM)
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
    meas(3:8,nofMeas+1:nofMeas+nofObs) = reshape([ephemerisObjects(i).ephemeris(1:nofMinBetweenObs:end).mee],6,nofObs);
    
    % Position and velocity in km
    posVel(1:3,1:nofObs) = [ephemerisObjects(i).ephemeris(1:nofMinBetweenObs:end).position] / 1e3;
    posVel(4:6,1:nofObs) = [ephemerisObjects(i).ephemeris(1:nofMinBetweenObs:end).velocity] / 1e3;
    % Covariance Position and velocity in km 
    posVelCov = reshape([ephemerisObjects(i).ephemeris(1:nofMinBetweenObs:end).covariance],6,6,nofObs) / 1e6;

    % Measurement noise
    RM(1:6,1:6,nofMeas+1:nofMeas+nofObs) = zeros(6,6,nofObs);
    for j=1:nofObs
        cov = posVelCov(:,:,j);
        try
            % Convert posvel covariance to mee covariance
            [~,RM(1:6,1:6,nofMeas+j)] = cartCov2meeCov(posVel(:,j),cov,GM);
        catch
            % If posvel covariance matrix is not Symmetric Positive Definite
            % then remove pos-vel cross-correlation terms
            cov(1:3,4:6) = 0;
            cov(4:6,1:3) = 0;
            [~,RM(1:6,1:6,nofMeas+j)] = cartCov2meeCov(posVel(:,j),cov,GM);
        end
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

% Scale semi-latus rectum from meters to kilometers
meas(3,:) = meas(3,:) / 1e3;
% RM(1,:,:) = RM(1,:,:) / 1e3;
% RM(:,1,:) = RM(:,1,:) / 1e3;

end

