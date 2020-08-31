clearvars -except localResultsDirPath
set(0, 'DefaultLineLineWidth', 1);
set(0,'defaultAxesFontSize',14);


% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_5obj_grav20x20';
% load(fullfile(localResultsDirPath,'workspace_runDensEstRadar_200508014710.mat'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_4obj_grav20x20';
% load(fullfile(localResultsDirPath,'workspace_runDensEstRadar_200511110150.mat'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_5obj_grav48x48_SRPSM';
% load(fullfile(localResultsDirPath,'ukf_rom_radar_JB2008_1999_2010_workspace_20200103_25d_5obj_200513145306.mat'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_5obj_grav20x20';
% load(fullfile(localResultsDirPath,'ukf_rom_radar_JB2008_1999_2010_workspace_20200103_25d_5obj_200514221532.mat'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_5obj_grav48x48_SRPSM';
% load(fullfile(localResultsDirPath,'ukf_rom_radar_JB2008_1999_2010_workspace_20200103_25d_5obj_200515074634.mat'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_13obj_grav48x48_SRPSM';
% load(fullfile(localResultsDirPath,'ukf_rom_radar_JB2008_1999_2010_workspace_20200103_25d_13obj_200615112634'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_11obj_grav48x48_SRPSM_fastProp_without7337_23278';
% load(fullfile(localResultsDirPath,'ukf_rom_radar_JB2008_1999_2010_workspace_20200103_25d_11obj_200621020933'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_11obj_grav48x48_SRPSM_fastProp';
% load(fullfile(localResultsDirPath,'ukf_rom_radar_JB2008_1999_2010_workspace_20200103_25d_11obj_200621014141'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_11obj_grav48x48_SRPSM_noOutliers/';
% load(fullfile(localResultsDirPath,'ukf_rom_radar_JB2008_1999_2010_workspace_20200103_25d_11obj_200623091517'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_13obj_grav48x48_SRPSM_noOutliers/';
% load(fullfile(localResultsDirPath,'ukf_rom_radar_JB2008_1999_2010_workspace_20200103_25d_13obj_200623210618'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_9obj_grav48x48_SRPSM_fastProp/';
% load(fullfile(localResultsDirPath,'ukf_rom_radar_JB2008_1999_2010_workspace_20200103_25d_9obj_200620105552'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/JB2008_1999_2010 - 2002 - CHAMP GOCE/HighFidelityDynamics_1percBCstd/20200103_25d_5objects/';
% load(fullfile(localResultsDirPath,'ukf_rom_tle_JB2008_1999_2010_workspace_20200103_25d_5obj_200625124925.mat'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/JB2008_1999_2010 - 2002 - CHAMP GOCE/HighFidelityDynamics_1percBCstd/20200103_25d_11objects/';
% load(fullfile(localResultsDirPath,'ukf_rom_tle_JB2008_1999_2010_workspace_20200103_25d_11obj_200626180841.mat'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_11obj_grav48x48_SRPSM_noOutliers_higherRrateCov';
% load(fullfile(localResultsDirPath,'ukf_rom_radar_JB2008_1999_2010_workspace_20200103_25d_11obj_200626151346.mat'));
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/JB2008_1999_2010 - 2002 - CHAMP GOCE/HighFidelityDynamics_1percBCstd/20200103_25d_15objects/';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/JB2008_1999_2010 - 2002 - CHAMP GOCE/HighFidelityDynamics_1percBCstd/20200103_25d_13objects/';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/JB2008_1999_2010 - 2002 - CHAMP GOCE/HighFidelityDynamics_1percBCstd/20200103_25d_12objects';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/JB2008_1999_2010 - 2002 - CHAMP GOCE/HighFidelityDynamics_1percBCstd/20200103_10d_19obj_200811132033';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/JB2008_1999_2010 - 2002 - CHAMP GOCE/HighFidelityDynamics_1percBCstd/20200103_10d_21obj_200811123552';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/JB2008_1999_2010 - 2002 - CHAMP GOCE/HighFidelityDynamics_1percBCstd/20200103_10d_20obj_200812030557';
localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/JB2008_1999_2010 - 2002 - CHAMP GOCE/HighFidelityDynamics_1percBCstd/20200103_10d_20obj_200812030557';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_5obj_grav48x48_SRPSM_noOutliers_higherProcessNoise/';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_5obj_grav48x48_SRPSM_noOutliers_higherProcessNoise2/';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_5obj_grav48x48_SRPSM_fastProp';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_13obj_grav48x48_SRPSM_noOutliers_higherProcessNoise';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_13obj_grav48x48_SRPSM_noOutliers_higherProcessNoise2';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_5diffObj_grav48x48_SRPSM_noOutliers_higherProcessNoise2';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/New Folder With Items';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/New Folder With Items 3';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/New Folder With Items 3';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_5obj_grav48x48_SRPSM_filteredAllOutliers_newProcessNoise';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_7obj_grav48x48_SRPSM_filteredAllOutliers_newProcessNoise_1rangeErr_2rangeRateErr';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_5obj_grav48x48_SRPSM_filteredAllOutliers_newProcessNoise_1rangeErr_2rangeRateErr';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_11obj_grav48x48_SRPSM_filteredAllOutliers_newProcessNoise_1rangeErr_2rangeRateErr';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/leolabs_rangeRate_6obj_grav48x48_SRPSM_filteredAllOutliers_newProcessNoise_1rangeErr_2rangeRateErr';
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/GPSmeas/LatestRun/';
% matfile = dir('/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/GPSmeas/LatestRun/ukf_rom_gps_NRLMSISE_1997_2008_workspace_20200501_10d_5obj_200718030237.mat');
% matfile = dir('/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/GPSmeas/LatestRun/ukf_rom_gps_NRLMSISE_1997_2008_workspace_20200501_10d_5obj_200718032244.mat');
% matfile = dir('/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/GPSmeas/LatestRun/20200501_10d_5obj_200718033540 GPS5m BC2%/ukf_rom_gps_NRLMSISE_1997_2008_workspace_20200501_10d_5obj_200718033540.mat');
% matfile = dir('/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/GPSmeas/LatestRun/ukf_rom_gps_NRLMSISE_1997_2008_workspace_20200511_10d_5obj_200718023911.mat');
% matfile = dir('/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/GPSmeas/LatestRun/ukf_rom_gps_NRLMSISE_1997_2008_workspace_20200521_10d_5obj_200718024833.mat');
% matfile = dir('/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/GPSmeas/LatestRun/ukf_rom_gps_JB2008_1999_2010_workspace_20200511_10d_5obj_200718052227.mat');
% matfile = dir('/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/GPSmeas/LatestRun/20200521_10d_5obj_200718031730/ukf_rom_gps_NRLMSISE_1997_2008_workspace_20200521_10d_5obj_200718031730.mat');
% localResultsDirPath = '/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/GPSmeas/20200501_10d_5obj_200718070647 5m BC10%';


% localResultsDirPath = ['/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/JB2008_1999_2010/' ...
%     'leolabs_rangeRate_21obj_grav48x48_SRPSM_movingWindowFilter_newProcessNoise_1rErr_2rrErr'];
%     'leolabs_rangeRate_20obj_grav48x48_SRPSM_movingWindowFilter_newProcessNoise_1rErr_2rrErr'];
%     'leolabs_rangeRate_21obj_grav48x48_SRPSM_movingWindowFilter_newProcessNoise_1rErr_2rrErr_full'];
%     'leolabs_rangeRate_19obj_grav48x48_SRPSM_movingWindowFilter_newProcessNoise_1rErr_2rrErr'];
% 'leolabs_rangeRate_12obj_grav48x48_SRPSM_movingWindowFilter_newProcessNoise_1rErr_2rrErr_higherRomCov'];
%     'leolabs_rangeRate_12obj_grav48x48_SRPSM_movingWindowFilter_newProcessNoise_1rErr_2rrErr'];
% localResultsDirPath = ['/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/RadarObs/' ...
%     'leolabs_rangeRate_12obj_originalRRunc'];
%     'leolabs_rangeRate_12obj_noRRfilter'];
%     'leolabs_rangeRate_12obj_lowerProcessNoise2'];
%     'leolabs_rangeRate_12obj_lowerProcessNoise'];
%     'leolabs_rangeRate_12obj_originalRRunc'];

% localResultsDirPath = ['/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/GPSmeas/' ...
%     'MS 20200417_7d_10obj_200814110127 newBC r20'];
%     'JB 20200417_7d_10obj_200814102740 newBC r20'];
%     'JB 20200501_10d_10obj_200723130615 newBC copy'];
%         'JB 20200417_9d_10obj_200806191915 newBC'];
%     'JB 20200417_7d_10obj_200805060725 newBC'];
%     'MS 20200417_7d_10obj_200805063551 newBC'];
%     'MS 20200501_30d_5obj_200801083746 oldBC'];
%     'JB 20200501_30d_5obj_200731200414 oldBC'];
%     'MS 20200417_10d_5obj_200730200620 newBC'];
%     'JB 20200501_30d_5obj_200723090440 newBC'];
%     'JB 20200501_6d_5obj_200729212329 oldBC 2%'];
%     'TLE JB 20200501_20d_10obj_200721032448 newBC'];
%     'MS 20200501_10d_10obj_200723124907 newBC'];
%     'JB 20200501_10d_10obj_200723130615 newBC'];
%     'MS 20200501_30d_5obj_200726164049 newBC'];
%     'JB 20200501_10d_5obj_200725090243 BC20%'];
%     'MS 20200501_10d_5obj_200725213121 BC20%'];
%     'TLE JB 20200416_11d_5obj_200725085922 newBC'];
%     'JB 20200417_10d_5obj_200725233051 newBC 2'];
%     'JB 20200417_10d_5obj_200725223433 newBC 1'];
%     'JB_20200511_10d_5obj_200718052227 5m BC10%'];
%     'JB 20200416_11d_5obj_200722051654 newBC'];
%     'MS 20200501_10d_10obj_200723124907 newBC'];
%     'JB 20200501_30d_5obj_200723090440 newBC'];
%     'TLE JB 20200501_20d_5obj_200720183745 newBC'];
%     'TLE JB 20200501_20d_15obj_200721200454 newBC'];
%     'MS 20200521_10d_5obj_200720141733 newBC'];
%     'JB 20200521_10d_5obj_200720135513 newBC'];
%     'MS 20200511_10d_5obj_200720134754 newBC'];
%     'JB 20200511_10d_5obj_200720125420 newBC'];
%     'MS 20200501_10d_5obj_200720121759 newBC'];
% 'JB 20200501_10d_5obj_200720134352 newBC'];
%     'MS_20200501_10d_5obj_200718033540 GPS5m BC2%'];
%     'MS_20200521_10d_5obj_200718024833 GPS5m'];
% 'MS_20200521_10d_5obj_200718070144 5m BC2%'];
% 'JB_20200501_10d_5obj_200718070647 5m BC10%'];

% 'MS_20200501_10d_5obj_200718122323 5m BC10%'];
% 'MS_20200511_10d_5obj_200718135226 5m BC10%'];
% 'MS_20200511_10d_5obj_200718071120 2m BC10%'];
% 'MS_20200511_10d_5obj_200718071840 5m BC2%'];

% matfiles = dir(fullfile(localResultsDirPath,'ukf_*.mat'));
% load(fullfile(matfiles.folder,matfiles.name));

matfiles = dir(fullfile(localResultsDirPath,'ukf_*.mat'));
load(fullfile(matfiles(1).folder,matfiles(1).name));

addpath('AstroFunctions');
addpath('/Users/davidgondelach/Dropbox (MIT)/Research-ROM/code/sqrt-ukf-rom-gps-2019/nrlmsise_matlab');
addpath('JB2008');
addpath('SpaceWeather');
spicePath = fullfile('/Users','davidgondelach','Documents','mice');
addpath( fullfile(spicePath,'src','mice') );
addpath( fullfile(spicePath,'lib') );

%%
computeJBMS = true;
SAVE_PLOTS = true;
DO_PREDICTION = false;

%% LOAD KERNELS, GRAVITY MODEL, EARTH ORIENTATION PARAMETERS AND SGP4
% Load SPICE
kernelpath  = fullfile('Data','kernel.txt');
loadSPICE(kernelpath);

%% Load GRACE/CHAMP data
champPath = fullfile('/Users/davidgondelach/Documents/PostDoc','CHAMPGRACEdata','champ','helium_heelis-gpi');
gracePath = fullfile('/Users/davidgondelach/Documents/PostDoc','CHAMPGRACEdata','grace','helium_heelis-gpi');
champPath_Mehta = fullfile('/Users/davidgondelach/Documents/PostDoc','CHAMPGRACEdata','champ_Mehta');
graceAPath_Mehta = fullfile('/Users/davidgondelach/Documents/PostDoc','CHAMPGRACEdata','graceA_Mehta');

champFilename = 'Density_3deg_';
graceAFilename = 'Density_graceA_3deg_';
graceBFilename = 'Density_graceB_3deg_';
champFilename_Mehta = 'CHAMP_Density_';
graceAFilename_Mehta = 'graceA_Density_';

nameCHAMP = 'CHAMP';
nameGRACEA = 'GRACE-A';
nameGRACEB = 'GRACE-B';

%% Load SWARM density data
swarmPath = '/Users/davidgondelach/Documents/PostDoc/SWARMdata/';

nameSWARM_A = 'A';
nameSWARM_B = 'B';
nameSWARM_C = 'C';
nameSWARM_A_full = 'SWARM A';
nameSWARM_B_full = 'SWARM B';
nameSWARM_C_full = 'SWARM C';
% SWARMdensityData = cdfread('/Users/davidgondelach/Downloads/SW_OPER_DNSAPOD_2__20200101T000000_20200101T235930_0201/SW_OPER_DNSAPOD_2__20200101T000000_20200101T235930_0201.cdf')


%% Load JB2008 space weather
[eopdata,SOLdata,DTCdata] = loadJB2008SWdata();

% Load NRLMSISE space weather data
SWpath = fullfile('Data','SW-All.txt');
[ SWmatDaily, ~ ] = inputSWnrlmsise( SWpath );
[ SWmatDailyTIEGCM, ~ ] = inputSWtiegcm( SWpath );

%% Date
% yearStart = yr;
% monthStart = mth;
% dayStart = dy;
% yearEnd = yr;
% monthEnd = mth;
% dayEnd = dayStart+nofDays-1;

% nofDays = 10;
% yearStart = 2002;
% monthStart = 8;
% dayStart = 11;
% yearEnd = 2002;
% monthEnd = 8;
% dayEnd = dayStart+nofDays;

% nofDays = 5;
% yearStart = 2004;
% monthStart = 11;
% dayStart = 7;
% yearEnd = 2004;
% monthEnd = 11;
% dayEnd = dayStart+nofDays;

% filenameBase = [localResultsDirPath ,'/', ROMmodel, '_', sprintf('%04d',yearStart), sprintf('%02d',monthStart), sprintf('%02d',dayStart), '_', num2str(nofDays), 'd_', nowTimeStr, '_'];

%% Load satellite data
% [champData] = getChampGraceData(champPath,champFilename,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);
% [graceAData] = getChampGraceData(gracePath,graceAFilename,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);
% [graceBData] = getChampGraceData(gracePath,graceBFilename,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);

% Mehta data
% [champData] = getChampGraceData_Mehta(champPath_Mehta,champFilename_Mehta,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);
% [graceAData] = getChampGraceData_Mehta(graceAPath_Mehta,graceAFilename_Mehta,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);


%% Estimated ROM state and covariance
% notSameEpoch = [diff(time)~=0 true];
% et_est = et0 + time(notSameEpoch);
% romState_est = X_est(end-r+1:end,notSameEpoch);
% romCov_est = Pv(end-r+1:end,notSameEpoch);
% BC_est = X_est(svs:svs:end-r,notSameEpoch);

%% Date
yearStart = yr;
monthStart = mth;
dayStart = dy;

%% Estimated ROM state and covariance
et_est = et0 + reshape(time,[1,length(time)]);
romState_est = X_est(end-r+1:end,:);
romCov_est = Pv(end-r+1:end,:);
BC_est = X_est(svs:svs:end-r,:);

for i=2:length(matfiles)
    load(fullfile(matfiles(i).folder,matfiles(i).name));
    
    %% Estimated ROM state and covariance
    et_est_next = et0 + reshape(time,[1,length(time)]);
    
%     [~,index] = min(abs(et_est_next - et_est(end)));
    [index] = find(et_est_next > et_est(end),1)
    
    romState_est_next = X_est(end-r+1:end,index+1:end);
    romCov_est_next = Pv(end-r+1:end,index+1:end);
    BC_est_next = X_est(svs:svs:end-r,index+1:end);
    
    et_est = [et_est, et_est_next(index+1:end)];
    romState_est = [romState_est, romState_est_next];
    romCov_est = [romCov_est, romCov_est_next];
    BC_est = [BC_est, BC_est_next];
end

notSameEpoch = [diff(et_est)~=0 , true];
et_est = et_est(notSameEpoch);
romState_est = romState_est(:,notSameEpoch);
romCov_est = romCov_est(:,notSameEpoch);
BC_est = BC_est(:,notSameEpoch);

%% Load SWARM density data
% % See: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/timout_c.html
% pic = 'DD-Mon-YYYY HR:MN:SC.##### ::RND'; % 'DD-Mon-YYYY HR:MN:SC.##### ::RND'
% 
% % Start date
% date0str = cspice_timout( et_est(1), pic );
% dateStart = datetime(date0str);
% yearStart = dateStart.Year;
% monthStart = dateStart.Month;
% dayStart = dateStart.Day;
% 
% % End date
% datefstr = cspice_timout( et_est(end), pic );
% dateEnd = datetime(datefstr);
% yearEnd = dateEnd.Year;
% monthEnd = dateEnd.Month;
% dayEnd = dateEnd.Day;

%% Date
yearEnd = yr;
monthEnd = mth;
dayEnd = dy+nofDays-1;

% SWARM data
[swarmAData] = getSwarmData(swarmPath,nameSWARM_A,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);
[swarmBData] = getSwarmData(swarmPath,nameSWARM_B,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);
% [swarmCData] = getSwarmData(swarmPath,nameSWARM_C,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);

%%
filenameBase = [localResultsDirPath ,'/', ROMmodel, '_', sprintf('%04d',yearStart), sprintf('%02d',monthStart), sprintf('%02d',dayStart), '_', num2str(nofDays), 'd_' ];

% %%
% nominalBC = 0.01376;
% yearStartDateNum = datenum(yearStart,1,1);
% plotTime = swarmAData.dates - yearStartDateNum + 1;
% [xxx,yyy] = plotDailyAveragedDensityError((et_est-et_est(1))/86400,nominalBC*ones(1,size(BC_est,2)),mean(BC_est)/1000);

%%
BCestPlot = figure;
for i=1:size(BC_est,1)
plot((et_est-et_est(1))/86400,BC_est(i,:)/median(BC_est(i,:))); hold on;
% plot((et_est-et_est(1))/86400,BC_est(i,:)/1000); hold on;
% plot((et_est-et_est(1))/86400+122,(BC_est(i,:)/1000 - nominalBC)./nominalBC*100); hold on;
end
legend(objectIDlabels,'Location','northeast');
xlabel('Day of year'); ylabel('BCest / medianBCest');
if SAVE_PLOTS
    savefig(BCestPlot,[filenameBase 'BCest.fig']);
end

%% Swarm A
goodDensities = swarmAData.densities < 1e-5;
swarmAData.longitudes = swarmAData.longitudes(goodDensities);
swarmAData.latitudes = swarmAData.latitudes(goodDensities);
swarmAData.altitudes = swarmAData.altitudes(goodDensities);
swarmAData.densities = swarmAData.densities(goodDensities);
swarmAData.localtimes = swarmAData.localtimes(goodDensities);
swarmAData.dates = swarmAData.dates(goodDensities);
[time_swarmA,rho_swarmA_real,rho_swarmA_rom,rho_swarmA_jb2,rho_swarmA_msise,SWdata_swarmA,rho_swarmA_rom_std] = getDensitiesRealRomJBMSISE(swarmAData,et_est,romState_est,r,F_U,M_U,SOLdata,DTCdata,eopdata,SWmatDaily,SWmatDailyTIEGCM,romCov_est,computeJBMS);
[orbitAvgErrChamp,dailyAvgErrChamp] = plotDensityComparison(time_swarmA,rho_swarmA_real,rho_swarmA_rom,rho_swarmA_jb2,rho_swarmA_msise,rho_swarmA_rom_std,swarmAData.latitudes,nameSWARM_A_full,SAVE_PLOTS,filenameBase);

meanStdRmsErrData(1,1) = mean(orbitAvgErrChamp(3,:)); % MS
meanStdRmsErrData(2,1) = mean(orbitAvgErrChamp(2,:)); % JB
meanStdRmsErrData(3,1) = mean(orbitAvgErrChamp(1,:)); % ROM
meanStdRmsErrData(1,2) = std(orbitAvgErrChamp(3,:)); % MS
meanStdRmsErrData(2,2) = std(orbitAvgErrChamp(2,:)); % JB
meanStdRmsErrData(3,2) = std(orbitAvgErrChamp(1,:)); % ROM
meanStdRmsErrData(1,3) = rms(orbitAvgErrChamp(3,:)); % MS
meanStdRmsErrData(2,3) = rms(orbitAvgErrChamp(2,:)); % JB
meanStdRmsErrData(3,3) = rms(orbitAvgErrChamp(1,:)); % ROM
meanStdRmsErrData(1,4) = mean(dailyAvgErrChamp(3,:)); % MS
meanStdRmsErrData(2,4) = mean(dailyAvgErrChamp(2,:)); % JB
meanStdRmsErrData(3,4) = mean(dailyAvgErrChamp(1,:)); % ROM
meanStdRmsErrData(1,5) = std(dailyAvgErrChamp(3,:)); % MS
meanStdRmsErrData(2,5) = std(dailyAvgErrChamp(2,:)); % JB
meanStdRmsErrData(3,5) = std(dailyAvgErrChamp(1,:)); % ROM
meanStdRmsErrData(1,6) = rms(dailyAvgErrChamp(3,:)); % MS
meanStdRmsErrData(2,6) = rms(dailyAvgErrChamp(2,:)); % JB
meanStdRmsErrData(3,6) = rms(dailyAvgErrChamp(1,:)); % ROM


rmsErrChamp_rom = rms((rho_swarmA_real-rho_swarmA_rom)./rho_swarmA_real*100);
rmsErrChamp_jb2 = rms((rho_swarmA_real-rho_swarmA_jb2)./rho_swarmA_real*100);
rmsErrChamp_msise = rms((rho_swarmA_real-rho_swarmA_msise)./rho_swarmA_real*100);

plotTime_champ = time_swarmA;
F10DSTplot = figure;
yyaxis left; plot(plotTime_champ,SWdata_swarmA(:,3)); ylabel('F10.7');
yyaxis right; plot(plotTime_champ,SWdata_swarmA(:,4)); ylabel('Kp');
xlabel('Day of year');  xlim([floor(plotTime_champ(1)) ceil(plotTime_champ(end))]);  xticks([floor(plotTime_champ(1)):5:ceil(plotTime_champ(end))]);
savefig(F10DSTplot,[filenameBase 'F10_Kp.fig']);


%% Swarm B
goodDensities_B = swarmBData.densities < 1e-5;
swarmBData.longitudes = swarmBData.longitudes(goodDensities_B);
swarmBData.latitudes = swarmBData.latitudes(goodDensities_B);
swarmBData.altitudes = swarmBData.altitudes(goodDensities_B);
swarmBData.densities = swarmBData.densities(goodDensities_B);
swarmBData.localtimes = swarmBData.localtimes(goodDensities_B);
swarmBData.dates = swarmBData.dates(goodDensities_B);
[time_swarmB,rho_swarmB_real,rho_swarmB_rom,rho_swarmB_jb2,rho_swarmB_msise,SWdata_swarmB,rho_swarmB_rom_std] = getDensitiesRealRomJBMSISE(swarmBData,et_est,romState_est,r,F_U,M_U,SOLdata,DTCdata,eopdata,SWmatDaily,SWmatDailyTIEGCM,romCov_est,computeJBMS);
[orbitAvgErrSwarmB,dailyAvgErrSwarmB] = plotDensityComparison(time_swarmB,rho_swarmB_real,rho_swarmB_rom,rho_swarmB_jb2,rho_swarmB_msise,rho_swarmB_rom_std,swarmBData.latitudes,nameSWARM_B_full,SAVE_PLOTS,filenameBase);

meanStdRmsErrData(4,1) = mean(orbitAvgErrSwarmB(3,:)); % MS
meanStdRmsErrData(5,1) = mean(orbitAvgErrSwarmB(2,:)); % JB
meanStdRmsErrData(6,1) = mean(orbitAvgErrSwarmB(1,:)); % ROM
meanStdRmsErrData(4,2) = std(orbitAvgErrSwarmB(3,:)); % MS
meanStdRmsErrData(5,2) = std(orbitAvgErrSwarmB(2,:)); % JB
meanStdRmsErrData(6,2) = std(orbitAvgErrSwarmB(1,:)); % ROM
meanStdRmsErrData(4,3) = rms(orbitAvgErrSwarmB(3,:)); % MS
meanStdRmsErrData(5,3) = rms(orbitAvgErrSwarmB(2,:)); % JB
meanStdRmsErrData(6,3) = rms(orbitAvgErrSwarmB(1,:)); % ROM
meanStdRmsErrData(4,4) = mean(dailyAvgErrSwarmB(3,:)); % MS
meanStdRmsErrData(5,4) = mean(dailyAvgErrSwarmB(2,:)); % JB
meanStdRmsErrData(6,4) = mean(dailyAvgErrSwarmB(1,:)); % ROM
meanStdRmsErrData(4,5) = std(dailyAvgErrSwarmB(3,:)); % MS
meanStdRmsErrData(5,5) = std(dailyAvgErrSwarmB(2,:)); % JB
meanStdRmsErrData(6,5) = std(dailyAvgErrSwarmB(1,:)); % ROM
meanStdRmsErrData(4,6) = rms(dailyAvgErrSwarmB(3,:)); % MS
meanStdRmsErrData(5,6) = rms(dailyAvgErrSwarmB(2,:)); % JB
meanStdRmsErrData(6,6) = rms(dailyAvgErrSwarmB(1,:)); % ROM

save([filenameBase 'meanStdRmsErrData.mat'],'meanStdRmsErrData');
save([filenameBase 'densities.mat'],'time_swarmA','rho_swarmA_real','rho_swarmA_rom','rho_swarmA_jb2','rho_swarmA_msise','SWdata_swarmA','rho_swarmA_rom_std', ...
                                    'time_swarmB','rho_swarmB_real','rho_swarmB_rom','rho_swarmB_jb2','rho_swarmB_msise','SWdata_swarmB','rho_swarmB_rom_std');

return;

%% Swarm C
goodDensities_C = swarmCData.densities < 1e-5;
swarmCData.longitudes = swarmCData.longitudes(goodDensities_C);
swarmCData.latitudes = swarmCData.latitudes(goodDensities_C);
swarmCData.altitudes = swarmCData.altitudes(goodDensities_C);
swarmCData.densities = swarmCData.densities(goodDensities_C);
swarmCData.localtimes = swarmCData.localtimes(goodDensities_C);
swarmCData.dates = swarmCData.dates(goodDensities_C);
[time_swarmC,rho_swarmC_real,rho_swarmC_rom,rho_swarmC_jb2,rho_swarmC_msise,SWdata,rho_swarmC_rom_std] = getDensitiesRealRomJBMSISE(swarmCData,etukf,romState_est,r,F_U,M_U,SOLdata,DTCdata,eopdata,SWmatDaily,SWmatDailyTIEGCM,romCov_est);
[orbitAvgErrSwarmC,dailyAvgErrSwarmC] = plotDensityComparison(time_swarmC,rho_swarmC_real,rho_swarmC_rom,rho_swarmC_jb2,rho_swarmC_msise,rho_swarmC_rom_std,swarmCData.latitudes,nameSWARM_C_full,SAVE_PLOTS,filenameBase);

meanStdRmsErrData(7,1) = mean(orbitAvgErrSwarmC(3,:)); % MS
meanStdRmsErrData(8,1) = mean(orbitAvgErrSwarmC(2,:)); % JB
meanStdRmsErrData(9,1) = mean(orbitAvgErrSwarmC(1,:)); % ROM
meanStdRmsErrData(7,2) = std(orbitAvgErrSwarmC(3,:)); % MS
meanStdRmsErrData(8,2) = std(orbitAvgErrSwarmC(2,:)); % JB
meanStdRmsErrData(9,2) = std(orbitAvgErrSwarmC(1,:)); % ROM
meanStdRmsErrData(7,3) = rms(orbitAvgErrSwarmC(3,:)); % MS
meanStdRmsErrData(8,3) = rms(orbitAvgErrSwarmC(2,:)); % JB
meanStdRmsErrData(9,3) = rms(orbitAvgErrSwarmC(1,:)); % ROM
meanStdRmsErrData(7,4) = mean(dailyAvgErrSwarmC(3,:)); % MS
meanStdRmsErrData(8,4) = mean(dailyAvgErrSwarmC(2,:)); % JB
meanStdRmsErrData(9,4) = mean(dailyAvgErrSwarmC(1,:)); % ROM
meanStdRmsErrData(7,5) = std(dailyAvgErrSwarmC(3,:)); % MS
meanStdRmsErrData(8,5) = std(dailyAvgErrSwarmC(2,:)); % JB
meanStdRmsErrData(9,5) = std(dailyAvgErrSwarmC(1,:)); % ROM
meanStdRmsErrData(7,6) = rms(dailyAvgErrSwarmC(3,:)); % MS
meanStdRmsErrData(8,6) = rms(dailyAvgErrSwarmC(2,:)); % JB
meanStdRmsErrData(9,6) = rms(dailyAvgErrSwarmC(1,:)); % ROM

%% GRACE-A
[time_graceA,rho_graceA_real,rho_graceA_rom,rho_graceA_jb2,rho_graceA_msise,~,rho_graceA_rom_std] = getDensitiesRealRomJBMSISE(graceAData,etukf,romState_est,r,F_U,M_U,SOLdata,DTCdata,eopdata,SWmatDaily,SWmatDailyTIEGCM,romCov_est);
[orbitAvgErrGraceA,dailyAvgErrGraceA] = plotDensityComparison(time_graceA,rho_graceA_real,rho_graceA_rom,rho_graceA_jb2,rho_graceA_msise,rho_graceA_rom_std,graceAData.latitudes,nameGRACEA,SAVE_PLOTS,filenameBase);

meanStdRmsErrData(4,1) = mean(orbitAvgErrGraceA(3,:)); % MS
meanStdRmsErrData(5,1) = mean(orbitAvgErrGraceA(2,:)); % JB
meanStdRmsErrData(6,1) = mean(orbitAvgErrGraceA(1,:)); % ROM
meanStdRmsErrData(4,2) = std(orbitAvgErrGraceA(3,:)); % MS
meanStdRmsErrData(5,2) = std(orbitAvgErrGraceA(2,:)); % JB
meanStdRmsErrData(6,2) = std(orbitAvgErrGraceA(1,:)); % ROM
meanStdRmsErrData(4,3) = mean(dailyAvgErrGraceA(3,:)); % MS
meanStdRmsErrData(5,3) = mean(dailyAvgErrGraceA(2,:)); % JB
meanStdRmsErrData(6,3) = mean(dailyAvgErrGraceA(1,:)); % ROM
meanStdRmsErrData(4,4) = std(dailyAvgErrGraceA(3,:)); % MS
meanStdRmsErrData(5,4) = std(dailyAvgErrGraceA(2,:)); % JB
meanStdRmsErrData(6,4) = std(dailyAvgErrGraceA(1,:)); % ROM

rmsErrGraceA_rom = rms((rho_graceA_real-rho_graceA_rom)./rho_graceA_real*100);
rmsErrGraceA_jb2 = rms((rho_graceA_real-rho_graceA_jb2)./rho_graceA_real*100);
rmsErrGraceA_msise = rms((rho_graceA_real-rho_graceA_msise)./rho_graceA_real*100);

%% GRACE-B
[time_graceB,rho_graceB_real,rho_graceB_rom,rho_graceB_jb2,rho_graceB_msise,~,rho_graceB_rom_std] = getDensitiesRealRomJBMSISE(graceBData,etukf,romState_est,r,F_U,M_U,SOLdata,DTCdata,eopdata,SWmatDaily,SWmatDailyTIEGCM,romCov_est);
[orbitAvgErrGraceB,dailyAvgErrGraceB] = plotDensityComparison(time_graceB,rho_graceB_real,rho_graceB_rom,rho_graceB_jb2,rho_graceB_msise,rho_graceB_rom_std,graceBData.latitudes,nameGRACEB,SAVE_PLOTS,filenameBase);

meanStdRmsErrData(7,1) = mean(orbitAvgErrGraceB(3,:)); % MS
meanStdRmsErrData(8,1) = mean(orbitAvgErrGraceB(2,:)); % JB
meanStdRmsErrData(9,1) = mean(orbitAvgErrGraceB(1,:)); % ROM
meanStdRmsErrData(7,2) = std(orbitAvgErrGraceB(3,:)); % MS
meanStdRmsErrData(8,2) = std(orbitAvgErrGraceB(2,:)); % JB
meanStdRmsErrData(9,2) = std(orbitAvgErrGraceB(1,:)); % ROM
meanStdRmsErrData(7,3) = mean(dailyAvgErrGraceB(3,:)); % MS
meanStdRmsErrData(8,3) = mean(dailyAvgErrGraceB(2,:)); % JB
meanStdRmsErrData(9,3) = mean(dailyAvgErrGraceB(1,:)); % ROM
meanStdRmsErrData(7,4) = std(dailyAvgErrGraceB(3,:)); % MS
meanStdRmsErrData(8,4) = std(dailyAvgErrGraceB(2,:)); % JB
meanStdRmsErrData(9,4) = std(dailyAvgErrGraceB(1,:)); % ROM

rmsErrGraceB_rom = rms((rho_graceB_real-rho_graceB_rom)./rho_graceB_real*100);
rmsErrGraceB_jb2 = rms((rho_graceB_real-rho_graceB_jb2)./rho_graceB_real*100);
rmsErrGraceB_msise = rms((rho_graceB_real-rho_graceB_msise)./rho_graceB_real*100);

%%
rmsErrData = [  rms(orbitAvgErrChamp'),rms(dailyAvgErrChamp'),rmsErrChamp_rom,rmsErrChamp_jb2,rmsErrChamp_msise;
                rms(orbitAvgErrGraceA'),rms(dailyAvgErrGraceA'),rmsErrGraceA_rom,rmsErrGraceA_jb2,rmsErrGraceA_msise;
                rms(orbitAvgErrGraceB'),rms(dailyAvgErrGraceB'),rmsErrGraceB_rom,rmsErrGraceB_jb2,rmsErrGraceB_msise];

save([filenameBase 'densityErrors.mat'],'orbitAvgErrChamp','dailyAvgErrChamp','orbitAvgErrGraceA','dailyAvgErrGraceA','orbitAvgErrGraceB','dailyAvgErrGraceB');

% cmap=colormap('lines');
% figure(51);
% plotOrbitAverage2(champData.latitudes,time_champ,rho_champ_real,rho_champ_rom,cmap(2,:),'-');
% figure(52);
% plotOrbitAverage2(champData.latitudes,time_champ,rho_champ_real,rho_champ_jb2,cmap(6,:),'-');
% figure(53);
% plotOrbitAverage2(champData.latitudes,time_champ,rho_champ_real,rho_champ_msise,cmap(4,:),'-');
% figure(51);
% plotOrbitAverage2(graceBData.latitudes,time_graceB,rho_graceB_real,rho_graceB_rom,cmap(2,:),'--');
% figure(52);
% plotOrbitAverage2(graceBData.latitudes,time_graceB,rho_graceB_real,rho_graceB_jb2,cmap(6,:),'--');
% figure(53);
% plotOrbitAverage2(graceBData.latitudes,time_graceB,rho_graceB_real,rho_graceB_msise,cmap(4,:),'--');

% plotOrbitAverage(champData.latitudes,time_champ,rho_champ_real,rho_champ_rom);
% plotOrbitAverage(champData.latitudes,time_champ,rho_champ_real,rho_champ_jb2);
% plotOrbitAverage(champData.latitudes,time_champ,rho_champ_real,rho_champ_msise);
% plotOrbitAverage(graceBData.latitudes,time_graceB,rho_graceB_real,rho_graceB_rom,rho_graceB_jb2,rho_graceB_msise);
% xlabel('Day of year'); ylabel('Orbit-averaged \rho [kg/m^3]');  xlim([floor(plotTime(1)) ceil(plotTime(end))]); xticks([floor(plotTime(1)):1:ceil(plotTime(end))]);
% legend(satName,'ROM','JB2008','NRLMSISE');

for k = 1:nop
xx_pv_est = zeros(6,size(meeMeas,2));
xx_pv_meas = zeros(6,size(meeMeas,2));
for j=1:size(meeMeas,2)
[pos,vel] = ep2pv(X_est((k-1)*svs+1:(k-1)*svs+6,j),GM_kms);
xx_pv_est(1:3,j) = pos;
xx_pv_est(4:6,j) = vel;
[pos,vel] = ep2pv(meeMeas((k-1)*6+1:(k-1)*6+6,j),GM_kms);
xx_pv_meas(1:3,j) = pos;
xx_pv_meas(4:6,j) = vel;
end
posErrors(k,:) = sqrt(sum( (xx_pv_est(1:3,:)-xx_pv_meas(1:3,:)) .^2,1));
end
avgPosErrPlot = figure;
plot(movmean(time,24)/86400,movmean(mean(posErrors(:,:)),24));
xlabel('Time [days]'); ylabel('Average position error');
savefig(avgPosErrPlot,[filenameBase 'avgPosErr.fig']);


if ~DO_PREDICTION
    return;
end

%% Load GRACE/CHAMP data in prediction period
yearStartPred = yearEnd;
monthStartPred = monthEnd;
dayStartPred = dayEnd+1;
yearEndPred = yearEnd;
monthEndPred = monthEnd;
dayEndPred = dayStartPred+10;

filenameBasePred = [filenameBase 'Pred_'];

[champDataPred] = getChampGraceData(champPath,champFilename,yearStartPred,monthStartPred,dayStartPred,yearEndPred,monthEndPred,dayEndPred);
[graceADataPred] = getChampGraceData(gracePath,graceAFilename,yearStartPred,monthStartPred,dayStartPred,yearEndPred,monthEndPred,dayEndPred);
[graceBDataPred] = getChampGraceData(gracePath,graceBFilename,yearStartPred,monthStartPred,dayStartPred,yearEndPred,monthEndPred,dayEndPred);



% 
% %% Predict ROM density
% startIndex = length(etukf);
% et_est_f = etukf(startIndex);
% et_pred = et_est_f + nofDays*86400;
% switch ROMmodel
%     case 'JB2008_1999_2010'
%         [romStateTime_pred] = predictROMdensity(romState_est(:,startIndex),et_est_f,et_pred,AC,BC,ROMmodel,eopdata,SOLdata,DTCdata,[]);
%     case 'TIEGCM_1997_2008_new'
%         [romStateTime_pred] = predictROMdensity(romState_est(:,startIndex),et_est_f,et_pred,AC,BC,ROMmodel,TIEGCM_SWdata,[],[],[]);
%     otherwise
%         [romStateTime_pred] = predictROMdensity(romState_est(:,startIndex),et_est_f,et_pred,AC,BC,ROMmodel,SWmatDaily,SWmatMonthlyPred,SWmatDailyTIEGCM,SWmatMonthlyPredTIEGCM);
% end
% 
% %% CHAMP prediction
% [time_champ_pred,rho_champ_real_pred,rho_champ_rom_pred,rho_champ_jb2_pred,rho_champ_msise_pred,SWdata_pred] = getDensitiesRealRomJBMSISE(champDataPred,romStateTime_pred(:,end),romStateTime_pred(:,1:r)',r,F_U,M_U,SOLdata,DTCdata,eopdata,SWmatDaily,SWmatDailyTIEGCM);
% plotDensityComparison(time_champ_pred,rho_champ_real_pred,rho_champ_rom_pred,rho_champ_jb2_pred,rho_champ_msise_pred,[],champDataPred.latitudes,nameCHAMP,SAVE_PLOTS,filenameBasePred);
% 
% % plotTime_champ_pred = (time_champ_pred-time_champ_pred(1))/86400;
% plotTime_champ_pred = time_champ_pred;
% F10DSTplot_pred = figure;
% yyaxis left; plot(plotTime_champ_pred,SWdata_pred(:,1)); ylabel('F10');
% % yyaxis right; plot(plotTime_champ_pred,SWdata_pred(:,2)); ylabel('DSTDTC');
% yyaxis right; plot(plotTime_champ_pred,SWdata_pred(:,4)); ylabel('Kp');
% xlabel('Day of year');  xlim([floor(plotTime_champ_pred(1)) ceil(plotTime_champ_pred(end))]);  xticks([floor(plotTime_champ_pred(1)):1:ceil(plotTime_champ_pred(end))]);
% savefig(F10DSTplot_pred,[filenameBasePred 'F10_Kp.fig']);
% 
% %% GRACE-A prediction
% [time_graceA_pred,rho_graceA_real_pred,rho_graceA_rom_pred,rho_graceA_jb2_pred,rho_graceA_msise_pred] = getDensitiesRealRomJBMSISE(graceADataPred,romStateTime_pred(:,end),romStateTime_pred(:,1:r)',r,F_U,M_U,SOLdata,DTCdata,eopdata,SWmatDaily,SWmatDailyTIEGCM);
% plotDensityComparison(time_graceA_pred,rho_graceA_real_pred,rho_graceA_rom_pred,rho_graceA_jb2_pred,rho_graceA_msise_pred,[],graceADataPred.latitudes,nameGRACEA,SAVE_PLOTS,filenameBasePred);
% 
% %% GRACE-B prediction
% [time_graceB_pred,rho_graceB_real_pred,rho_graceB_rom_pred,rho_graceB_jb2_pred,rho_graceB_msise_pred] = getDensitiesRealRomJBMSISE(graceBDataPred,romStateTime_pred(:,end),romStateTime_pred(:,1:r)',r,F_U,M_U,SOLdata,DTCdata,eopdata,SWmatDaily,SWmatDailyTIEGCM);
% plotDensityComparison(time_graceB_pred,rho_graceB_real_pred,rho_graceB_rom_pred,rho_graceB_jb2_pred,rho_graceB_msise_pred,[],graceBDataPred.latitudes,nameGRACEB,SAVE_PLOTS,filenameBasePred);


% rho_champ_rom_pred_err = (rho_champ_rom_pred-rho_champ_real_pred);
% rho_champ_jb2_pred_err = (rho_champ_jb2_pred-rho_champ_real_pred);

% rho_champ_rom_pred_absErr = abs(rho_champ_rom_pred-rho_champ_real_pred);
% rho_champ_jb2_pred_absErr = abs(rho_champ_jb2_pred-rho_champ_real_pred);

% plotTime_champ_pred = (time_champ_pred-time_champ_pred(1))/3600;
% figure;
% plot(plotTime_champ_pred,rho_champ_real_pred); hold on;
% plot(plotTime_champ_pred,rho_champ_rom_pred); hold on;
% plot(plotTime_champ_pred,rho_champ_jb2_pred); hold on;
% xlabel('Time [hours]'); ylabel('\rho [kg/m^3]'); xlim([0 240]); xticks([0:24:240]);
% legend('CHAMP','ROM','JB2008');

% figure;
% plot(plotTime_champ_pred,abs(rho_champ_rom_pred-rho_champ_real_pred)); hold on;
% plot(plotTime_champ_pred,abs(rho_champ_jb2_pred-rho_champ_real_pred)); hold on;
% xlabel('Time [hours]'); ylabel('Error \rho [kg/m^3]'); xlim([0 240]); xticks([0:24:240]);
% legend('ROM','JB2008');

% figure;
% plot(movmean(plotTime_champ_pred,100),movmean(abs(rho_champ_rom_pred-rho_champ_real_pred),100)); hold on;
% plot(movmean(plotTime_champ_pred,100),movmean(abs(rho_champ_jb2_pred-rho_champ_real_pred),100)); hold on;
% xlabel('Time [hours]'); ylabel('Moving average error \rho [kg/m^3]'); xlim([0 240]); xticks([0:24:240]);
% legend('ROM','JB2008');

% plotOrbitAverage(champDataPred.latitudes,plotTime_champ_pred,rho_champ_rom_pred_absErr,rho_champ_jb2_pred_absErr);
% xlabel('Time [hours]'); ylabel('Average error \rho per orbit [kg/m^3]'); xlim([0 240]); xticks([0:24:240]);
% legend('ROM','JB2008');

% plotOrbitAverage(champDataPred.latitudes,plotTime_champ_pred,rho_champ_rom_pred_err,rho_champ_jb2_pred_err);
% xlabel('Time [hours]'); ylabel('Average error \rho per orbit [kg/m^3]'); xlim([0 240]); xticks([0:24:240]);
% legend('ROM','JB2008');
