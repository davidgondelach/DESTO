clearvars
set(0, 'DefaultLineLineWidth', 1);
set(0,'defaultAxesFontSize',14);

% localResultsDirPath = '/Users/davidgondelach/Documents/GitRepos/DensityEstimation/Tests/workspace_UKF_191003093954.mat';
localResultsDirPath = '/Users/davidgondelach/Documents/GitRepos/DensityEstimation/Tests/ukf_workspace_git_200127.mat';

% matfile = dir(fullfile(localResultsDirPath,'ukf_*.mat'));
load(fullfile(localResultsDirPath));

addpath('AstroFunctions');
addpath('nrlmsise_matlab');
addpath('JacchiaBowmanAtmosphericDensityModel');
spicePath = fullfile('/Users','davidgondelach','Documents','mice');
addpath( fullfile(spicePath,'src','mice') );
addpath( fullfile(spicePath,'lib') );

%%
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

%% Load JB2008 space weather
[eopdata,SOLdata,DTCdata] = loadJB2008SWdata();

% Load NRLMSISE space weather data
SWpath = fullfile('Data','SW-All.txt');
[ SWmatDaily, ~ ] = inputSWnrlmsise( SWpath );
[ SWmatDailyTIEGCM, ~ ] = inputSWtiegcm( SWpath );

%% Date
yearStart = yr;
monthStart = mth;
dayStart = dy;
yearEnd = yr;
monthEnd = mth;
dayEnd = dayStart+nofDays-1;

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

filenameBase = [localResultsDirPath ,'/', ROMmodel, '_', sprintf('%04d',yearStart), sprintf('%02d',monthStart), sprintf('%02d',dayStart), '_', num2str(nofDays), 'd_' ];
% filenameBase = [localResultsDirPath ,'/', ROMmodel, '_', sprintf('%04d',yearStart), sprintf('%02d',monthStart), sprintf('%02d',dayStart), '_', num2str(nofDays), 'd_', nowTimeStr, '_'];

%% Load satellite data
[champData] = getChampGraceData(champPath,champFilename,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);
[graceAData] = getChampGraceData(gracePath,graceAFilename,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);
[graceBData] = getChampGraceData(gracePath,graceBFilename,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);

% Mehta data
% [champData] = getChampGraceData_Mehta(champPath_Mehta,champFilename_Mehta,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);
% [graceAData] = getChampGraceData_Mehta(graceAPath_Mehta,graceAFilename_Mehta,yearStart,monthStart,dayStart,yearEnd,monthEnd,dayEnd);

%% Estimated ROM state and covariance
 etukf = et0 + time;
romState_est = X_est(end-r+1:end,:);
romCov_est = Pv(end-r+1:end,:);
BC_est = X_est(svs:svs:end-r,:);

BCestPlot = figure;
for i=1:size(BC_est,1)
plot((etukf-etukf(1))/86400,BC_est(i,:)/median(BC_est(i,:))); hold on;
end
legend(objectIDlabels,'Location','northeast');
xlabel('Day of year'); ylabel('BCest / medianBCest');
savefig(BCestPlot,[filenameBase 'BCest.fig']);

%% CHAMP
% [time_champ,rho_champ_real,rho_champ_rom,rho_champ_jb2,SWdata,rho_champ_rom_std] = getDensitiesRealRomJB(champData,etukf,romState_est,r,F_U,M_U,SOLdata,DTCdata,eopdata,romCov_est);
[time_champ,rho_champ_real,rho_champ_rom,rho_champ_jb2,rho_champ_msise,SWdata,rho_champ_rom_std] = getDensitiesRealRomJBMSISE(champData,etukf,romState_est,r,F_U,M_U,SOLdata,DTCdata,eopdata,SWmatDaily,SWmatDailyTIEGCM,romCov_est);
[orbitAvgErrChamp,dailyAvgErrChamp] = plotDensityComparison(time_champ,rho_champ_real,rho_champ_rom,rho_champ_jb2,rho_champ_msise,rho_champ_rom_std,champData.latitudes,nameCHAMP,SAVE_PLOTS,filenameBase);

meanStdErrData(1,1) = mean(orbitAvgErrChamp(3,:)); % MS
meanStdErrData(2,1) = mean(orbitAvgErrChamp(2,:)); % JB
meanStdErrData(3,1) = mean(orbitAvgErrChamp(1,:)); % ROM
meanStdErrData(1,2) = std(orbitAvgErrChamp(3,:)); % MS
meanStdErrData(2,2) = std(orbitAvgErrChamp(2,:)); % JB
meanStdErrData(3,2) = std(orbitAvgErrChamp(1,:)); % ROM
meanStdErrData(1,3) = mean(dailyAvgErrChamp(3,:)); % MS
meanStdErrData(2,3) = mean(dailyAvgErrChamp(2,:)); % JB
meanStdErrData(3,3) = mean(dailyAvgErrChamp(1,:)); % ROM
meanStdErrData(1,4) = std(dailyAvgErrChamp(3,:)); % MS
meanStdErrData(2,4) = std(dailyAvgErrChamp(2,:)); % JB
meanStdErrData(3,4) = std(dailyAvgErrChamp(1,:)); % ROM

rmsErrChamp_rom = rms((rho_champ_real-rho_champ_rom)./rho_champ_real*100);
rmsErrChamp_jb2 = rms((rho_champ_real-rho_champ_jb2)./rho_champ_real*100);
rmsErrChamp_msise = rms((rho_champ_real-rho_champ_msise)./rho_champ_real*100);

% plotTime_champ = (time_champ-time_champ(1))/86400;
plotTime_champ = time_champ;
F10DSTplot = figure;
yyaxis left; plot(plotTime_champ,SWdata(:,3)); ylabel('F10.7');
yyaxis right; plot(plotTime_champ,SWdata(:,4)); ylabel('Kp');
xlabel('Day of year');  xlim([floor(plotTime_champ(1)) ceil(plotTime_champ(end))]);  xticks([floor(plotTime_champ(1)):5:ceil(plotTime_champ(end))]);
savefig(F10DSTplot,[filenameBase 'F10_Kp.fig']);

%% GRACE-A
[time_graceA,rho_graceA_real,rho_graceA_rom,rho_graceA_jb2,rho_graceA_msise,~,rho_graceA_rom_std] = getDensitiesRealRomJBMSISE(graceAData,etukf,romState_est,r,F_U,M_U,SOLdata,DTCdata,eopdata,SWmatDaily,SWmatDailyTIEGCM,romCov_est);
[orbitAvgErrGraceA,dailyAvgErrGraceA] = plotDensityComparison(time_graceA,rho_graceA_real,rho_graceA_rom,rho_graceA_jb2,rho_graceA_msise,rho_graceA_rom_std,graceAData.latitudes,nameGRACEA,SAVE_PLOTS,filenameBase);

meanStdErrData(4,1) = mean(orbitAvgErrGraceA(3,:)); % MS
meanStdErrData(5,1) = mean(orbitAvgErrGraceA(2,:)); % JB
meanStdErrData(6,1) = mean(orbitAvgErrGraceA(1,:)); % ROM
meanStdErrData(4,2) = std(orbitAvgErrGraceA(3,:)); % MS
meanStdErrData(5,2) = std(orbitAvgErrGraceA(2,:)); % JB
meanStdErrData(6,2) = std(orbitAvgErrGraceA(1,:)); % ROM
meanStdErrData(4,3) = mean(dailyAvgErrGraceA(3,:)); % MS
meanStdErrData(5,3) = mean(dailyAvgErrGraceA(2,:)); % JB
meanStdErrData(6,3) = mean(dailyAvgErrGraceA(1,:)); % ROM
meanStdErrData(4,4) = std(dailyAvgErrGraceA(3,:)); % MS
meanStdErrData(5,4) = std(dailyAvgErrGraceA(2,:)); % JB
meanStdErrData(6,4) = std(dailyAvgErrGraceA(1,:)); % ROM

rmsErrGraceA_rom = rms((rho_graceA_real-rho_graceA_rom)./rho_graceA_real*100);
rmsErrGraceA_jb2 = rms((rho_graceA_real-rho_graceA_jb2)./rho_graceA_real*100);
rmsErrGraceA_msise = rms((rho_graceA_real-rho_graceA_msise)./rho_graceA_real*100);

%% GRACE-B
[time_graceB,rho_graceB_real,rho_graceB_rom,rho_graceB_jb2,rho_graceB_msise,~,rho_graceB_rom_std] = getDensitiesRealRomJBMSISE(graceBData,etukf,romState_est,r,F_U,M_U,SOLdata,DTCdata,eopdata,SWmatDaily,SWmatDailyTIEGCM,romCov_est);
[orbitAvgErrGraceB,dailyAvgErrGraceB] = plotDensityComparison(time_graceB,rho_graceB_real,rho_graceB_rom,rho_graceB_jb2,rho_graceB_msise,rho_graceB_rom_std,graceBData.latitudes,nameGRACEB,SAVE_PLOTS,filenameBase);

meanStdErrData(7,1) = mean(orbitAvgErrGraceB(3,:)); % MS
meanStdErrData(8,1) = mean(orbitAvgErrGraceB(2,:)); % JB
meanStdErrData(9,1) = mean(orbitAvgErrGraceB(1,:)); % ROM
meanStdErrData(7,2) = std(orbitAvgErrGraceB(3,:)); % MS
meanStdErrData(8,2) = std(orbitAvgErrGraceB(2,:)); % JB
meanStdErrData(9,2) = std(orbitAvgErrGraceB(1,:)); % ROM
meanStdErrData(7,3) = mean(dailyAvgErrGraceB(3,:)); % MS
meanStdErrData(8,3) = mean(dailyAvgErrGraceB(2,:)); % JB
meanStdErrData(9,3) = mean(dailyAvgErrGraceB(1,:)); % ROM
meanStdErrData(7,4) = std(dailyAvgErrGraceB(3,:)); % MS
meanStdErrData(8,4) = std(dailyAvgErrGraceB(2,:)); % JB
meanStdErrData(9,4) = std(dailyAvgErrGraceB(1,:)); % ROM

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
% xlabel('Day of year'); ylabel('Orbit-averaged \rho [kg/m^2]');  xlim([floor(plotTime(1)) ceil(plotTime(end))]); xticks([floor(plotTime(1)):1:ceil(plotTime(end))]);
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
