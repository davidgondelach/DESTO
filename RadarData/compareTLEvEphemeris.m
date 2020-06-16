clearvars;
clearvars -global;

addpath( 'AstroFunctions' );
addpath( 'Estimation' );
addpath( 'TLEdata' );
addpath( 'RadarData' );
spicePath = fullfile('/Users/davidgondelach/Documents','mice');
addpath( fullfile(spicePath,'src','mice') );
addpath( fullfile(spicePath,'lib') );


%% SETTINGS
% Specify the date, reduced-order model, reduction order and objects to be
% used to estimation here.

% Estimation window
% Continuous data from:2020-01-03T06:03:06.870935 , till:2020-01-28T06:10:43.642242
yr      = 2020; % Year
mth     = 1;    % Month
dy      = 3;    % Day
hr      = 6;
mn      = 0;
sc     = 0;
nofDays = 25;   % Number of days

% Use high fidelity dynamical model
highFidelity = true;

% Reduced-order model
ROMmodel = 'JB2008_1999_2010';  % Name of reduced-order density model: JB2008_1999_2010, NRLMSISE_1997_2008 or TIEGCM_1997_2008
r  = 10;                        % Reduced order

% NORAD catalog IDs of objects used for estimation
% Default: 17 objects: [63;165;614;2153;2622;4221;6073;7337;8744;12138;12388;14483;20774;23278;27391;27392;26405]
% selectedObjects = [63;165;614;2153;2622;4221;6073;7337;8744;12138;12388;14483;20774;23278;27391;27392;26405]; % TLE
% selectedObjects = [614;2153;2622;4221;12138;750;2016;2389;6073;7337;8744;12388;14483;20774;23278]; % Radar
selectedObjects = [614;2153;2622;4221;12138]; % Radar

%% Load 
% Load Earth orientation parameters (needed for TEME to ECI transformation)
global EOPMat
EOPpath = fullfile('Data','EOP-All.txt');
[ EOPMat ] = inputEOP_Celestrak( EOPpath );

% Setup the SGP4 propagator.
loadSGP4();

global mu  % Earth gravitational parameter according to SGP4 model [km^3 s^-2]
global GM  % Earth gravitational parameter according to accurate gravity model [m^3 s^-2]
GM_kms = GM*1e-9; % Earth gravitational parameter according to accurate gravity model [km^3 s^-2]

%% Datetime

% Julian date
jd0 = juliandate(datetime(yr,mth,dy,hr,mn,sc));
jdf = juliandate(datetime(yr,mth,dy+nofDays,hr,mn,sc));

% Initial Ephemeris Time (ET): ET is the number of seconds past the
% epoch of the J2000 reference frame in the time system known as
% Barycentric Dynamical Time (TDB).
et0  = cspice_str2et(strcat([num2str(jed2date(jd0),'%d %d %d %d %d %.10f') 'UTC']));
etf  = cspice_str2et(strcat([num2str(jed2date(jdf),'%d %d %d %d %d %.10f') 'UTC']));


%% SET PATHS
measurementsPath = '/Users/davidgondelach/Documents/RadarData/LeoLabsData/';
ephemerisPath = '/Users/davidgondelach/Documents/RadarData/LeoLabsEphemeris/';

% %% Get Radar measurement data
% [radarMeasurementsObjects] = getRadarMeasurements(measurementsPath,selectedObjects);
% [radarStations] = getRadarStations(measurementsPath);
% [measRangeRangeRate,RMrangeRangeRate] = generateRangeRangeRateObsFromRadarMeas(radarMeasurementsObjects,radarStations,et0,etf);

%% Get Radar ephemeris data
[ephemerisObjects] = getEphemeris(ephemerisPath,selectedObjects);

% Position and velocity measurements
[measEphem,Rmeas] = generateObservationsFromEphemeris(ephemerisObjects,et0,etf);

%% Get TLEs
[yrf, mthf, dyf, ~, ~, ~] = datevec(jdf+30-1721058.5);  % End date of TLE collection window 
getTLEsFromSingleFile = true; % If true: all TLEs are loaded from file named "estimationObjects.tle" else TLEs are loaded from individual files named "[NORADID].tle"
[objects] = getTLEsForEstimation(yr, mth, 1, yrf, mthf, dyf, selectedObjects, getTLEsFromSingleFile);

%% Compare 
TLEdata = struct('position',{},'velocity',{},'diffPos',{},'diffPosRTN',{});
% TLEstates = [];

nop = length(selectedObjects);
for i=1:nop
    noradID = selectedObjects(i);
    index = find([objects.noradID]==noradID);
    
    meas = measEphem(:,measEphem(2,:)==i);
    nofMeas = size(meas,2);
    rj2000 = zeros(3,nofMeas);
    vj2000 = zeros(3,nofMeas);
    diffR = zeros(1,nofMeas);
    rr_Diff_RTN = zeros(3,nofMeas);
    for j=1:size(meas,2)
        % Observation epoch
        et = meas(1,j);
        jdatestr    = cspice_et2utc( et, 'J', 12 );
        obsJdate    = str2double(jdatestr(4:end)); % Cut trailing 'JD ' off from string
        % Find nearest newer TLE
        satrecIndex = find([objects(index).satrecs.jdsatepoch]>=obsJdate,1,'first');
        diffObsTLEEpochMinutes = (obsJdate - objects(index).satrecs(satrecIndex).jdsatepoch) * 24*60;
        % Compute SGP4 state at epoch
        [~, rtemeObs ,vtemeObs] = sgp4( objects(index).satrecs(satrecIndex), diffObsTLEEpochMinutes );
        % Convert to J2000
        [rj2000(:,j), vj2000(:,j)] = convertTEMEtoJ2000(rtemeObs', vtemeObs', obsJdate);
        
        rr_Diff = rj2000(:,j) - meas(3:5,j);
    
        diffR(j) = norm(rr_Diff);

        [cart2rtnMatrix] = computeCart2RTNMatrix(meas(3:5,j), meas(6:8,j));
        rr_Diff_RTN(:,j) = cart2rtnMatrix*rr_Diff;
    end
    TLEdata(i).position = rj2000;
    TLEdata(i).velocity = vj2000;
    TLEdata(i).diffPos = diffR;
    TLEdata(i).diffPosRTN = rr_Diff_RTN;
end

%%
for i=1:nop
    hours = ( measEphem(1,measEphem(2,:)==i) - measEphem(1,find(measEphem(2,:)==i,1)) ) / 3600;
    
    figure; hold on;
    plot(hours,TLEdata(i).diffPosRTN(1,:),'.');
    plot(hours,TLEdata(i).diffPosRTN(2,:),'.');
    plot(hours,TLEdata(i).diffPosRTN(3,:),'.');
    legend('Radial','Transverse','Normal');
    
end
