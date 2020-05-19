function [rangeRangeRate] = extractSingleRangeRangeRate(Xp,meas,svs,GM_kms)
% - Return range and range rate of single object
% Xp = predicted states in MEE, BC and ROM
% meas = measurement: 1) object number, 2) range, 3) range rate, 4-9) Observer location

% Get object number from measurement
objectNumber = meas(1);
% Get object state
objectStateMEE = Xp(svs*(objectNumber-1)+1:svs*(objectNumber-1)+6,:);

% Convert to ECI
objectStateECI = ep2cart(objectStateMEE, GM_kms);

% Get observer position and velocity in ECI
rSiteECI = meas(end-5:end-3);
vSiteECI = meas(end-2:end);

% Compute range and range rate
rangeRangeRate = computeRangeAndRangeRate(objectStateECI(1:3,:),objectStateECI(4:6,:),rSiteECI,vSiteECI);

end