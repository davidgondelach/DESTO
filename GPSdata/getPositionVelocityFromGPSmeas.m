function [posVelAtET0] = getPositionVelocityFromGPSmeas(gpsDataPath,object,BC,et0)

% Planet satellites with GPS data
planetSkysats = struct;
planetSkysats(1).planetID = 's1';     planetSkysats(1).noradID = 39418; % SKYSAT 1 [MISLABEL]
planetSkysats(2).planetID = 's2';     planetSkysats(2).noradID = 40072; % SKYSAT 2 [MISLABEL]
planetSkysats(3).planetID = 's3';     planetSkysats(3).noradID = 41601; % SKYSAT C1
planetSkysats(4).planetID = 's4';     planetSkysats(4).noradID = 41773; % SKYSAT C2
planetSkysats(5).planetID = 's103';   planetSkysats(5).noradID = 41774; % SKYSAT C3
planetSkysats(6).planetID = 's104';   planetSkysats(6).noradID = 41771; % SKYSAT C4
planetSkysats(7).planetID = 's105';   planetSkysats(7).noradID = 41772; % SKYSAT C5
planetSkysats(8).planetID = 's106';   planetSkysats(8).noradID = 42992; % SKYSAT C6
planetSkysats(9).planetID = 's107';   planetSkysats(9).noradID = 42991; % SKYSAT C7
planetSkysats(10).planetID = 's108';  planetSkysats(10).noradID = 42990; % SKYSAT C8
planetSkysats(11).planetID = 's109';  planetSkysats(11).noradID = 42989; % SKYSAT C9
planetSkysats(12).planetID = 's110';  planetSkysats(12).noradID = 42988; % SKYSAT C10
planetSkysats(13).planetID = 's111';  planetSkysats(13).noradID = 42987; % SKYSAT C11
planetSkysats(14).planetID = 's112';  planetSkysats(14).noradID = 43797; % SKYSAT C12
planetSkysats(15).planetID = 's113';  planetSkysats(15).noradID = 43802; % SKYSAT C13

%%
index = find([planetSkysats.noradID]==object);

% Get GPS data
[gpsData] = getGPSdataFromJsonFileET(gpsDataPath,planetSkysats(index).planetID,planetSkysats(index).noradID,et0-86400,et0+86400);

[~,closestGPSmeas] = min( abs([gpsData(:).tET] - et0) );
gps_state = gpsData(closestGPSmeas).xx_J2000;
et_state = gpsData(closestGPSmeas).tET;

deltaET = et0 - et_state;
if deltaET == 0
    posVelAtET0 = gps_state;
    return
end

stateBC = [gps_state; BC];

settingsJB2000.drag=3;
settingsJB2000.thirdbody = 1;
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[~,stateBCAtET0] = ode113(@(t,x) computeAcceleration(t,x,et_state,settingsJB2000,[]),[0,deltaET],stateBC,opts);

posVelAtET0 = stateBCAtET0(end,1:6)';

end

