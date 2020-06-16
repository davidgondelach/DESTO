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

dy      = 3;    % Day
hr      = 6;
nofDays = 31;   % Number of days

% Use high fidelity dynamical model
highFidelity = true;

% Reduced-order model
ROMmodel = 'JB2008_1999_2010';  % Name of reduced-order density model: JB2008_1999_2010, NRLMSISE_1997_2008 or TIEGCM_1997_2008
r  = 10;                        % Reduced order

% NORAD catalog IDs of objects used for estimation
% Default: 17 objects: [63;165;614;2153;2622;4221;6073;7337;8744;12138;12388;14483;20774;23278;27391;27392;26405]
% selectedObjects = [63;165;614;2153;2622;4221;6073;7337;8744;12138;12388;14483;20774;23278;27391;27392;26405]; % TLE
selectedObjects = [614;2153;2622;4221;12138;750;2016;2389;6073;7337;8744;12388;14483;20774;23278]; % Radar

selectedObjects = sortrows(selectedObjects);

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


%% Get Radar measurement data
[radarMeasurementsObjects] = getRadarMeasurements(measurementsPath,selectedObjects);
[radarStations] = getRadarStations(measurementsPath);
[measRangeRangeRate,RMrangeRangeRate] = generateRangeRangeRateObsFromRadarMeas(radarMeasurementsObjects,radarStations,et0,etf);

%% Get TLEs
[yrf, mthf, dyf, ~, ~, ~] = datevec(jdf+30-1721058.5);  % End date of TLE collection window 
getTLEsFromSingleFile = true; % If true: all TLEs are loaded from file named "estimationObjects.tle" else TLEs are loaded from individual files named "[NORADID].tle"
[objects] = getTLEsForEstimation(yr, mth, 1, yrf, mthf, dyf, selectedObjects, getTLEsFromSingleFile);

%% Compare 
TLEdata = struct('position',{},'velocity',{},'range',{},'rangerate',{});
% TLEstates = [];

nop = length(selectedObjects);
for i=1:nop
    meas = measRangeRangeRate(:,measRangeRangeRate(2,:)==i);
    nofMeas = size(meas,2);
    rj2000 = zeros(3,nofMeas);
    vj2000 = zeros(3,nofMeas);
    rangeAndRangeRate = zeros(2,nofMeas);
    for j=1:size(meas,2)
        % Observation epoch
        et = meas(1,j);
        jdatestr    = cspice_et2utc( et, 'J', 12 );
        obsJdate    = str2double(jdatestr(4:end)); % Cut trailing 'JD ' off from string
        % Find nearest newer TLE
        satrecIndex = find([objects(i).satrecs.jdsatepoch]>=obsJdate,1,'first');
        diffObsTLEEpochMinutes = (obsJdate - objects(i).satrecs(satrecIndex).jdsatepoch) * 24*60;
        % Compute SGP4 state at epoch
        [~, rtemeObs ,vtemeObs] = sgp4( objects(i).satrecs(satrecIndex), diffObsTLEEpochMinutes );
        % Convert to J2000
        [rj2000(:,j), vj2000(:,j)] = convertTEMEtoJ2000(rtemeObs', vtemeObs', obsJdate);
        
        rSiteECI = meas(5:7,j);
        vSiteECI = meas(8:10,j);
        rangeAndRangeRate(:,j) = computeRangeAndRangeRate(rj2000(:,j), vj2000(:,j),rSiteECI,vSiteECI);
    end
    TLEdata(i).position = rj2000;
    TLEdata(i).velocity = vj2000;
    TLEdata(i).range = rangeAndRangeRate(1,:);
    TLEdata(i).rangerate = rangeAndRangeRate(2,:);
end

%%
objectIDlabels = cell(1, nop);
for i=1:nop
    objectIDlabels(i) = {num2str(objects(i).noradID)};
end

for i=1:nop
    hours = ( measRangeRangeRate(1,measRangeRangeRate(2,:)==i) - measRangeRangeRate(1,find(measRangeRangeRate(2,:)==i,1)) ) / 3600;
    measRange = measRangeRangeRate(3,measRangeRangeRate(2,:)==i);
    measRangeRate = measRangeRangeRate(4,measRangeRangeRate(2,:)==i);
    errRange = RMrangeRangeRate(1,measRangeRangeRate(2,:)==i).^(0.5);
    errRangeRate = RMrangeRangeRate(2,measRangeRangeRate(2,:)==i).^(0.5);
    
    figure;
    subplot(2,2,1); hold on;
    plot(TLEdata(i).range,'.');
    plot(measRange,'.');
    title(num2str(radarMeasurementsObjects(i).noradID));
    subplot(2,2,3); hold on;
    plot(TLEdata(i).range - measRange,'.');
    plot(errRange,'-','Color',[0.8 0.8 0.8]);
    plot(-errRange,'-','Color',[0.8 0.8 0.8]);
    
    subplot(2,2,2); hold on;
    plot(hours,TLEdata(i).rangerate,'.');
    plot(hours,measRangeRate,'.');
    xlabel('Hours');
    subplot(2,2,4); hold on;
    plot(hours,TLEdata(i).rangerate - measRangeRate,'.');
    plot(hours,errRangeRate,'-','Color',[0.8 0.8 0.8]);
    plot(hours,-errRangeRate,'-','Color',[0.8 0.8 0.8]);
    xlabel('Hours');
end

%%
figure;
for i=1:nop
    hours = ( measRangeRangeRate(1,measRangeRangeRate(2,:)==i) - et0 ) / 3600;
    measRange = measRangeRangeRate(3,measRangeRangeRate(2,:)==i);
    
    plot(hours,TLEdata(i).range - measRange,'.');
    hold on;
    xlabel('Hours');
end
legend(objectIDlabels)
