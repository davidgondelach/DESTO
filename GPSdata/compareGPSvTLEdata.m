% % Download GPS files
% yearStr = '2020';
% monthStr = '06';
% for i=1:18
%   websave(['gps_data_', yearStr, monthStr, num2str(i,'%02d'),'.zip'],['http://ephemerides.planet-labs.com/gps_data_', yearStr, monthStr, num2str(i,'%02d'),'.zip']);
% end

clearvars;

spicePath = fullfile('/Users','davidgondelach','Documents','mice');
addpath( 'AstroFunctions' );
addpath( 'Analysis' );
addpath( 'GPSdata' );
addpath( 'GPSdata/ECEF2ECI' );
addpath( 'TLEdata' );
addpath( fullfile(spicePath,'src','mice') );
addpath( fullfile(spicePath,'lib') );

% Constants
GM = 398600.4415;
% Setup the SGP4 propagator.
global tumin mu radiusearthkm xke j2 j3 j4 j3oj2 opsmode whichconst
opsmode = 'i';
whichconst = 72;
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc( whichconst );
xpdotp   =  1440.0 / (2.0*pi); % Conversion factor between SGP4 and TLE mean motion units [rev/day]/[rad/min]
% Clear cspice memory
cspice_kclear;
% Load SPK, PCK, LSK kernels
kernelpath  = fullfile('Data','kernel.txt');
cspice_furnsh( kernelpath );
% Load Earth orientation parameters (needed for TEME to ECI transformation)
global EOPMat
EOPpath = fullfile(pwd,'Data','EOP-All.txt');
[ EOPMat ] = inputEOP_Celestrak( EOPpath );

%%
% objectPlanetID = '0E0E'; objectNORADID = 41609;
% objectPlanetID = '0F4A'; objectNORADID = 42861;
% objectPlanetID = '0F06'; objectNORADID = 42995;
planetSkysats = struct;
planetSkysats(1).objectPlanetID = 's1';    planetSkysats(1).objectNORADID = 39418; % SKYSAT 1 [MISLABEL]
planetSkysats(2).objectPlanetID = 's2';    planetSkysats(2).objectNORADID = 40072; % SKYSAT 2 [MISLABEL]
planetSkysats(3).objectPlanetID = 's3';    planetSkysats(3).objectNORADID = 41601; % SKYSAT C1
planetSkysats(4).objectPlanetID = 's4';    planetSkysats(4).objectNORADID = 41773; % SKYSAT C2
planetSkysats(5).objectPlanetID = 's103';  planetSkysats(5).objectNORADID = 41774; % SKYSAT C3
planetSkysats(6).objectPlanetID = 's104';  planetSkysats(6).objectNORADID = 41771; % SKYSAT C4
planetSkysats(7).objectPlanetID = 's105';  planetSkysats(7).objectNORADID = 41772; % SKYSAT C5
planetSkysats(8).objectPlanetID = 's106';  planetSkysats(8).objectNORADID = 42992; % SKYSAT C6
planetSkysats(9).objectPlanetID = 's107';  planetSkysats(9).objectNORADID = 42991; % SKYSAT C7
planetSkysats(10).objectPlanetID = 's108';  planetSkysats(10).objectNORADID = 42990; % SKYSAT C8
planetSkysats(11).objectPlanetID = 's109';  planetSkysats(11).objectNORADID = 42989; % SKYSAT C9
planetSkysats(12).objectPlanetID = 's110';  planetSkysats(12).objectNORADID = 42988; % SKYSAT C10
planetSkysats(13).objectPlanetID = 's111';  planetSkysats(13).objectNORADID = 42987; % SKYSAT C11
planetSkysats(14).objectPlanetID = 's112';  planetSkysats(14).objectNORADID = 43797; % SKYSAT C12
planetSkysats(15).objectPlanetID = 's113';  planetSkysats(15).objectNORADID = 43802; % SKYSAT C13

objNr = 9;
% for objNr = 9
%% Get TLE data
filename = fullfile('TLEdata',[num2str(planetSkysats(objNr).objectNORADID),'.tle']);
[objectTLEs] = getTLEs(filename);

%     jdatesGPS_full = [];
%     rr_Diff_RTN_full = [];
%     goodOnes_full = [];
%     UsedTLEs_full = [];
%     vv_Diff_RTN_full = [];

%% Get GPS data
gpsDataPath = '/Users/davidgondelach/Documents/PostDoc/GPSdata';
%     for day = 1:30
%         tic;
%         disp(day);
[gpsData,filepath] = getGPSdataFromJsonFile(gpsDataPath,planetSkysats(objNr).objectPlanetID,planetSkysats(objNr).objectNORADID,2020,5,1,2020,5,30);

% Remove half of measurements closer than 5 seconds from each other
smallerThan5sec = find(diff([gpsData(:).tET]) <= 5);
gpsData(smallerThan5sec(1:2:end)) = [];
% Remove half of measurements closer than 10 seconds from each other
smallerThan10sec = find(diff([gpsData(:).tET]) <= 10);
gpsData(smallerThan10sec(1:2:end)) = [];
% Remove half of measurements closer than 20 seconds from each other
smallerThan20sec = find(diff([gpsData(:).tET]) <= 20);
gpsData(smallerThan20sec(1:2:end)) = [];
% Remove all measurements closer than 20 seconds from each other
% smallerThan20sec = find(diff([gpsData(:).tET]) <= 20);
% gpsData(smallerThan20sec(1:end)) = [];
figure;
histogram(diff([gpsData(:).tET]));
hold on; histogram(diff([gpsData(:).tET]));

%%
diffEpochMinutes = zeros(1,length(gpsData));
diffR = zeros(1,length(gpsData));
jdatesGPS = zeros(1,length(gpsData));
rr_Diff_RTN = zeros(length(gpsData),3);
vv_Diff_RTN = zeros(length(gpsData),3);
for i=1:length(gpsData)
    jdatestr    = cspice_et2utc( gpsData(i).tET, 'J', 12 );
    jdate       = str2double(jdatestr(4:end)); % Cut leading 'JD ' off from string
    
    % Find nearest-newer TLE and get state in J2000
    satrecIndex = find([objectTLEs.satrecs.jdsatepoch] >= jdate, 1, 'first');
    diffEpochMinutes(i) = (jdate - objectTLEs.satrecs(satrecIndex).jdsatepoch) * 24*60;
    [~, rr_TLE_TEME ,vv_TLE_TEME] = sgp4( objectTLEs.satrecs(satrecIndex), diffEpochMinutes(i) );
    [rr_TLE_J2000, vv_TLE_J2000] = convertTEMEtoJ2000(rr_TLE_TEME', vv_TLE_TEME', jdate);
    %     coeTLE(i,:) = pv2po(rr_TLE_J2000, vv_TLE_J2000, GM);
    
    rr_Diff = rr_TLE_J2000 - gpsData(i).xx_J2000(1:3);
    vv_Diff = vv_TLE_J2000 - gpsData(i).xx_J2000(4:6);
    
    diffR(i) = norm(rr_Diff);
    jdatesGPS(i) = jdate;
    
    [cart2rtnMatrix] = computeCart2RTNMatrix(gpsData(i).xx_J2000(1:3), gpsData(i).xx_J2000(4:6));
    rr_Diff_RTN(i,:) = cart2rtnMatrix*rr_Diff;
    vv_Diff_RTN(i,:) = cart2rtnMatrix*vv_Diff;
end
[UsedTLEs] = find([objectTLEs.satrecs.jdsatepoch] >= min(jdatesGPS) & [objectTLEs.satrecs.jdsatepoch] <= max(jdatesGPS));

%         goodOnes = diffR < 3;
%         gpsData = gpsData(goodOnes);
%         matFile = strcat(filepath(1:end-5),'.mat');
%         save(matFile,'gpsData');

%         jdatesGPS_full = [jdatesGPS_full,jdatesGPS];
%         rr_Diff_RTN_full = [rr_Diff_RTN_full;rr_Diff_RTN];
%         vv_Diff_RTN_full = [vv_Diff_RTN_full;vv_Diff_RTN];
%         goodOnes_full = [goodOnes_full,goodOnes];
%         UsedTLEs_full = [UsedTLEs_full,UsedTLEs];
%         toc
%     end

figure;
subplot(2,1,1);
plot((jdatesGPS-jdatesGPS(1))*24,rr_Diff_RTN(:,1),(jdatesGPS-jdatesGPS(1))*24,rr_Diff_RTN(:,2),(jdatesGPS-jdatesGPS(1))*24,rr_Diff_RTN(:,3)); hold on;
plot((jdatesGPS(~goodOnes)-jdatesGPS(1))*24,rr_Diff_RTN(~goodOnes,1),'o',(jdatesGPS(~goodOnes)-jdatesGPS(1))*24,rr_Diff_RTN(~goodOnes,2),'o',(jdatesGPS(~goodOnes)-jdatesGPS(1))*24,rr_Diff_RTN(~goodOnes,3),'o'); hold on;
if ~isempty(UsedTLEs)
    plot(([objectTLEs.satrecs(UsedTLEs).jdsatepoch]-jdatesGPS(1))*24,0,'b.');
end
xlabel('Hours'); ylabel('Position error [km]');
legend('Radial','Transverse','Normal'); legend('hide');
title('Unfiltered');
subplot(2,1,2);
plot((jdatesGPS(goodOnes)-jdatesGPS(1))*24,rr_Diff_RTN(goodOnes,1),(jdatesGPS(goodOnes)-jdatesGPS(1))*24,rr_Diff_RTN(goodOnes,2),(jdatesGPS(goodOnes)-jdatesGPS(1))*24,rr_Diff_RTN(goodOnes,3)); hold on;
if ~isempty(UsedTLEs)
    plot(([objectTLEs.satrecs(UsedTLEs).jdsatepoch]-jdatesGPS(1))*24,0,'b.');
end
xlabel('Hours'); ylabel('Position error [km]');
legend('Radial','Transverse','Normal'); legend('hide');
title('Filtered');
%     savefig(['GPSvTLE_pos_', num2str(planetSkysats(objNr).objectNORADID), '_', planetSkysats(objNr).objectPlanetID, '_2020May', num2str(day), '_filtered.fig']);

figure;
subplot(2,1,1);
plot((jdatesGPS-jdatesGPS(1))*24,vv_Diff_RTN(:,1),(jdatesGPS-jdatesGPS(1))*24,vv_Diff_RTN(:,2),(jdatesGPS-jdatesGPS(1))*24,vv_Diff_RTN(:,3)); hold on;
plot((jdatesGPS(~goodOnes)-jdatesGPS(1))*24,vv_Diff_RTN(~goodOnes,1),'o',(jdatesGPS(~goodOnes)-jdatesGPS(1))*24,vv_Diff_RTN(~goodOnes,2),'o',(jdatesGPS(~goodOnes)-jdatesGPS(1))*24,vv_Diff_RTN(~goodOnes,3),'o'); hold on;
if ~isempty(UsedTLEs)
    plot(([objectTLEs.satrecs(UsedTLEs).jdsatepoch]-jdatesGPS(1))*24,0,'b.');
end
xlabel('Hours'); ylabel('Velocity error [km/s]');
legend('Radial','Transverse','Normal'); legend('hide');
title('Unfiltered');
subplot(2,1,2);
plot((jdatesGPS(goodOnes)-jdatesGPS(1))*24,vv_Diff_RTN(goodOnes,1),(jdatesGPS(goodOnes)-jdatesGPS(1))*24,vv_Diff_RTN(goodOnes,2),(jdatesGPS(goodOnes)-jdatesGPS(1))*24,vv_Diff_RTN(goodOnes,3)); hold on;
if ~isempty(UsedTLEs)
    plot(([objectTLEs.satrecs(UsedTLEs).jdsatepoch]-jdatesGPS(1))*24,0,'b.');
end
xlabel('Hours'); ylabel('Position error [km]');
legend('Radial','Transverse','Normal'); legend('hide');
title('Filtered');
%     savefig(['GPSvTLE_vel_', num2str(planetSkysats(objNr).objectNORADID), '_', planetSkysats(objNr).objectPlanetID, '_2020May', num2str(day), '_filtered.fig']);


% end
return

%%

figure;
plot(diffEpochMinutes,diffR,'.');
xlabel('Minutes since TLE epoch'); ylabel('Position error [km]');

figure;
plot(jdatesGPS,diffR);
xlabel('Julian date'); ylabel('Position error [km]');

figure;
plot(jdatesGPS,rr_Diff_RTN(:,1),jdatesGPS,rr_Diff_RTN(:,2),jdatesGPS,rr_Diff_RTN(:,3)); hold on;
plot([objectTLEs.satrecs(UsedTLEs).jdsatepoch],0,'b.');
xlabel('Julian date'); ylabel('Position error [km]');
legend('Radial','Transverse','Normal');

meanRadialErr = mean(rr_Diff_RTN(:,1))
stdRadialErr = std(rr_Diff_RTN(:,1))
meanTransvErr = mean(rr_Diff_RTN(:,2))
stdTransvErr = std(rr_Diff_RTN(:,2))
meanNormalErr = mean(rr_Diff_RTN(:,3))
stdNormalErr = std(rr_Diff_RTN(:,3))

figure;
% plot(jdatesGPS,coeTLE(:,1));
plot([objectTLEs.satrecs(UsedTLEs).jdsatepoch],[objectTLEs.satrecs(UsedTLEs).a]*radiusearthkm,'.');
xlabel('Julian date'); ylabel('a_{TLE} [km]');