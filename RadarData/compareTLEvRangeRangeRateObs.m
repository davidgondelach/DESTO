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
sc      = 0;
nofDays = 25;   % Number of days

dy      = 1;    % Day
hr      = 0;
nofDays = 31;   % Number of days

% Use high fidelity dynamical model
highFidelity = true;

% Reduced-order model
ROMmodel = 'JB2008_1999_2010';  % Name of reduced-order density model: JB2008_1999_2010, NRLMSISE_1997_2008 or TIEGCM_1997_2008
r  = 10;                        % Reduced order

% NORAD catalog IDs of objects used for estimation
% Default: 17 objects: [63;165;614;2153;2622;4221;6073;7337;8744;12138;12388;14483;20774;23278;27391;27392;26405]
% selectedObjects = [63;165;614;2153;2622;4221;6073;7337;8744;12138;12388;14483;20774;23278;27391;27392;26405]; % TLE
% selectedObjects = [614;2153;2622;4221;12138;750;2016;2389;6073;7337;8744;12388;14483;20774;23278]; % Radar
% selectedObjects = [614;2153;2622;4221;12138]; % Radar
% selectedObjects = [614;2153;2389;4221;7337;8744;12138;12388;14483;20774;23278]; % Radar: 11 objects
selectedObjects = [614;2153;2389;2622;4221;7337;8744;12138;12388;14483;20774;23278]; % Radar: 12 objects
% selectedObjects = [7337;8744;12388;14483;23278]; % Radar: 11 objects
% selectedObjects = [12388]; % Radar
% selectedObjects = [23278]; % Radar
% selectedObjects = [12138;6073;8744;12388;14483;20774;23278]; % Radar


selectedObjects = sortrows(selectedObjects);

%% Load 
% Load SPICE kernels and ephemerides
kernelpath  = fullfile('Data','kernel.txt');
loadSPICE(kernelpath);

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

% % Filter out first three measurements of passage
% nop = length(selectedObjects);
% for i=1:nop
%     %     instrum1 = {radarMeasurementsObjects(i).measurements(1:end-1).instrument};
%     %     instrum2 = {radarMeasurementsObjects(i).measurements(2:end).instrument};
%     %     sameStation = strcmp(instrum1,instrum2);
%     length(radarMeasurementsObjects(i).measurements)
%     for j = 1:3
%         isOutlier = [true diff( [radarMeasurementsObjects(i).measurements(:).measuredAtET] ) > 500];
%         isOutlier = isOutlier | [diff( [radarMeasurementsObjects(i).measurements(:).measuredAtET] ) > 500 true];
%         
%         radarMeasurementsObjects(i).measurements(isOutlier) = [];
%     end
%     length(radarMeasurementsObjects(i).measurements)
% end


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
    % NORAD ID
    ID = selectedObjects(i);
    % Orbital data (duplicate to ensure data is available)
    index = find([objects.noradID]==ID);
    
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
        satrecIndex = find([objects(index).satrecs.jdsatepoch]>=obsJdate,1,'first');
        diffObsTLEEpochMinutes = (obsJdate - objects(index).satrecs(satrecIndex).jdsatepoch) * 24*60;
        % Compute SGP4 state at epoch
        [~, rtemeObs ,vtemeObs] = sgp4( objects(index).satrecs(satrecIndex), diffObsTLEEpochMinutes );
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
    % NORAD ID
    ID = selectedObjects(i);
    % Orbital data (duplicate to ensure data is available)
    index = find([objects.noradID]==ID);
    objectIDlabels(i) = {num2str(objects(index).noradID)};
end

for i=1:nop
    hours = ( measRangeRangeRate(1,measRangeRangeRate(2,:)==i) - et0 ) / 3600;
%     hours = measRangeRangeRate(1,measRangeRangeRate(2,:)==i) ;
    measRange = measRangeRangeRate(3,measRangeRangeRate(2,:)==i);
    measRangeRate = measRangeRangeRate(4,measRangeRangeRate(2,:)==i);
    errRange = squeeze( RMrangeRangeRate(1,1,measRangeRangeRate(2,:)==i).^(0.5) );
    errRangeRate = squeeze( RMrangeRangeRate(2,2,measRangeRangeRate(2,:)==i).^(0.5) );
    
    diffRange = measRange - TLEdata(i).range;
%     diffRangeOutlier = isoutlier(diffRange,'movmedian',7);
%     diffRangeOutlier = abs(diffRange) > 0.3 | diffRangeOutlier;
    diffRangeOutlier = abs(diffRange) > 0.5;
%     length(diffRange)
%     sum(diffRangeOutlier)
    sum(diffRangeOutlier)/length(diffRange)*100
    
%     diffRangeOutlier = [radarMeasurementsObjects(i).measurements(:).snr] < 3;
%     diffRangeOutlier = [true diff( measRangeRangeRate(1,measRangeRangeRate(2,:)==i) ) > 1e9];
%     diffRangeOutlier = diffRangeOutlier | [diff( measRangeRangeRate(1,measRangeRangeRate(2,:)==i) ) > 500 true];
%     sum(diffRangeOutlier)/length(diffRange)*100
    
    
    diffRangeRate = measRangeRate - TLEdata(i).rangerate;
%     diffRangeRateOutlier = isoutlier(diffRangeRate,'movmedian',7);
%     diffRangeRateOutlier = abs(diffRangeRate) > 0.02 | diffRangeRateOutlier;
    diffRangeRateOutlier = abs(diffRangeRate) > 0.025;
%     length(diffRangeRate)
%     sum(diffRangeRateOutlier)
    sum(diffRangeRateOutlier)/length(diffRangeRate)*100
    
%     diffRangeRateOutlier = [radarMeasurementsObjects(i).measurements(:).snr] < 3;
%     sum(diffRangeRateOutlier)/length(diffRangeRate)*100
%     diffRangeRateOutlier = [true diff( measRangeRangeRate(1,measRangeRangeRate(2,:)==i) ) > 1e9];


    elevation = zeros(length(radarMeasurementsObjects(i).measurements),1);
    snr = zeros(length(radarMeasurementsObjects(i).measurements),1);
    for k=1:length(radarMeasurementsObjects(i).measurements)
        elevation(k) = [radarMeasurementsObjects(i).measurements(k).corrected.elevation];
        snr(k) = radarMeasurementsObjects(i).measurements(k).snr;
    end
    
    figure;
%     subplot(2,2,1); hold on;
%     plot(hours,TLEdata(i).range,'.');
%     plot(hours,measRange,'.');
%     title(num2str(radarMeasurementsObjects(i).noradID));
    subplot(2,1,1); hold on;
%     plot(hours,diffRange,'.');
    scatter(hours,diffRange,10,elevation);
    plot(hours(diffRangeOutlier),diffRange(diffRangeOutlier),'r.');
    plot(hours(diffRangeRateOutlier),diffRange(diffRangeRateOutlier),'c.');
    plot(hours,3*errRange,'-','Color',[0.8 0.8 0.8]);
    plot(hours,-3*errRange,'-','Color',[0.8 0.8 0.8]);
    title(objectIDlabels(i));
    xlabel('Hours');
    ylabel('Range diff [km]');
    
%     subplot(2,1,2); hold on;
%     plot(hours,TLEdata(i).rangerate,'.');
%     plot(hours,measRangeRate,'.');
%     xlabel('Hours');
    subplot(2,1,2); hold on;
%     plot(hours,diffRangeRate,'.');
    scatter(hours,diffRangeRate,10,elevation);
    plot(hours(diffRangeOutlier),diffRangeRate(diffRangeOutlier),'c.');
    plot(hours(diffRangeRateOutlier),diffRangeRate(diffRangeRateOutlier),'r.');
    plot(hours,3*errRangeRate,'-','Color',[0.8 0.8 0.8]);
    plot(hours,-3*errRangeRate,'-','Color',[0.8 0.8 0.8]);
    title('meas - TLE');
    xlabel('Hours');
    ylabel('Range rate diff [km/s]');
%     savefig(['RadarvTLE_range_rangeRate_', objectIDlabels{i}, '_2020Jan.fig']);
    
    etMeas = measRangeRangeRate(1,measRangeRangeRate(2,:)==i);
    etOutliers = etMeas( abs(TLEdata(i).range - measRange) > 0.7 | abs(TLEdata(i).rangerate - measRangeRate) > 0.2 );
    outliers = ismember([radarMeasurementsObjects(i).measurements(:).measuredAtET],  etOutliers);
    outliersMeasID = [radarMeasurementsObjects(i).measurements(outliers).id];
    
    print(num2str(radarMeasurementsObjects(i).noradID));
    elevation = zeros(length(radarMeasurementsObjects(i).measurements),1);
    for k=1:length(radarMeasurementsObjects(i).measurements)
        elevation(k) = [radarMeasurementsObjects(i).measurements(k).corrected.elevation];
        radarMeasurementsObjects(i).measurements(k).rangeErr = diffRange(k);
        radarMeasurementsObjects(i).measurements(k).rangeRateErr = diffRangeRate(k);
    end
%     figure;
% %     plot(elevation, abs(diffRange), '.');
%     scatter(elevation, abs(diffRange),25,hours);
end

%%
rangeErrors = [];
rangeRateErrors = [];
elevations = [];
instruments = [];
snrs = [];
for i=1:nop
    rangeErrors = [rangeErrors, [radarMeasurementsObjects(i).measurements(:).rangeErr]];
    rangeRateErrors = [rangeRateErrors, [radarMeasurementsObjects(i).measurements(:).rangeRateErr]];
    
    elevation = zeros(length(radarMeasurementsObjects(i).measurements),1);
    instrum = zeros(length(radarMeasurementsObjects(i).measurements),1);
    snr = zeros(length(radarMeasurementsObjects(i).measurements),1);
    for k=1:length(radarMeasurementsObjects(i).measurements)
        elevation(k) = [radarMeasurementsObjects(i).measurements(k).corrected.elevation];
        snr(k) = radarMeasurementsObjects(i).measurements(k).snr;
        instrument = radarMeasurementsObjects(i).measurements(k).instrument;
        switch instrument
            case 'pfisr'
                instrum(k) = 1;
            case 'msr'
                instrum(k) = 2;
            case 'ksr1'
                instrum(k) = 3;
            case 'ksr2'
                instrum(k) = 4;
        end
    end
    elevations = [elevations;elevation];
    instruments = [instruments;instrum];
    snrs = [snrs;snr];
end

%%
figure; histogram((measRangeRangeRate(1,:)-et0)/86400,[0:1/2:31]);
xlabel('Days since 1 Jan 2020 0:00'); ylabel('Number of measurements');
figure;
for i=1:4
    subplot(2,2,i); hold on;
    histogram((measRangeRangeRate(1,(instruments==i))-et0)/86400,[0:1/2:31]);
    xlabel('Days since 1 Jan 2020 0:00'); ylabel('Number of measurements');
end
    

%%
figure;
for i=1:4
    subplot(1,2,1); hold on;
% scatter(elevations(instruments==i), (rangeErrors(instruments==i)),25,measRange(instruments==i));
scatter(elevations(instruments==i), (rangeErrors(instruments==i)),25);
switch i
    case 1
        title('pfisr');
    case 2
        title('msr');
    case 3
        title('ksr1');
    case 4
        title('ksr2');
end
    subplot(1,2,2); hold on;
% scatter(elevations(instruments==i), (rangeRateErrors(instruments==i)),25,measRange(instruments==i));
scatter(elevations(instruments==i), (rangeRateErrors(instruments==i)),25);
switch i
    case 1
        title('pfisr');
    case 2
        title('msr');
    case 3
        title('ksr1');
    case 4
        title('ksr2');
end
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

figure
for i=1:nop
    errRange = squeeze( RMrangeRangeRate(1,1,measRangeRangeRate(2,:)==i).^(0.5) );
    errRangeRate = squeeze( RMrangeRangeRate(2,2,measRangeRangeRate(2,:)==i).^(0.5) );
    plot(errRange); hold on;
end

%%
for i=1:nop
    seconds = ( measRangeRangeRate(1,measRangeRangeRate(2,:)==i) - et0 );
    
    figure;
    histogram(log10(diff(seconds)));
    xticks(log10([0.001 0.01 0.1 1 10 60 300 900 3600 7200 14400 28800 86400]));
    xticklabels({'0.001', '0.01', '0.1', '1', '10', '1min', '5min', '15min', '1hr', '2hr', '4hr', '8hr', '1day'});
    hold on;
    xlabel('Hours');
    
    
end

%%
numberOfMeasPerPass = [];
for i=1:nop
    seconds = ( measRangeRangeRate(1,measRangeRangeRate(2,:)==i) - et0 );
    
    numberOfPasses(i,1) = sum(diff(seconds)>5*60);
    
    count = 0;
    currentJ = 0;
    for j=1:length(radarMeasurementsObjects(i).measurements)-1
        if ~strncmp(radarMeasurementsObjects(i).measurements(j).instrument,radarMeasurementsObjects(i).measurements(j+1).instrument,3)
            count = count+1;
            numberOfMeasPerPass(end+1) = j-currentJ;
            currentJ = j;
        end
    end
    numberOfMeasPerPass(end+1) = length(radarMeasurementsObjects(i).measurements) - currentJ;
    numberOfPasses(i,2) = count;
    
    
end
figure;
histogram(numberOfMeasPerPass);
mean(numberOfMeasPerPass)

%%
for i=1:nop    
    figure;
    histogram(log10([radarMeasurementsObjects(i).measurements(:).snr]));
end