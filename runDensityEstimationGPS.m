function runDensityEstimationGPS(yr,mth,dy,hr,mn,sc,nofDays,ROMmodel,r,selectedObjects,varargin)
%runDensityEstimationTLE - Estimate thermospheric density using TLE data.
% 
%  runDensityEstimationTLE(yr,mth,dy,nofDays,ROMmodel,r,selectedObjects) -
%     Estimate thermospheric density using TLE data: 
%     1) load TLE data, 2) load ballistic coefficient data, 3) generate 
%     observations from TLE data, 4) generate reduced-order density model,
%     5) initialize reduced-order density state, 6) set initial orbital
%     state guesses, 7) set measurement and process noise, 8) set initial  
%     state covariance, 9) run unscented Kalman filter to simultaneously
%     estimate the orbital states, ballistic coefficients and thermospheric
%     density from TLE-derived orbit observations.
% 
%  runDensityEstimationTLE(yr,mth,dy,nofDays,ROMmodel,r,selectedObjects,plotFigures)
%     Estimate thermospheric density using TLE data and plot results.
% 
%     yr                - start date: year
%     mth               - start date: month
%     dy                - start date: day
%     nofDays           - number of days of estimation window
%     ROMmodel          - name of reduced-order density model
%     r                 - dimension of reduced-order density model
%     selectedObjects   - NORAD IDs of objects used for estimation
%     plotFigures       - boolean: if true then plot results
%
% 
%     Copyright (C) 2020 by David Gondelach
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
%  Author: David Gondelach
%  Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
%  email: davidgondelach@gmail.com
%  Jan 2020; Last revision: 31-Jan-2020
%
%  Reference:
%  D.J. Gondelach and R. Linares, "Real-Time Thermospheric Density
%  Estimation Via Two-Line-Element Data Assimilation", Space Weather, 2020
%  https://doi.org/10.1029/2019SW002356 or https://arxiv.org/abs/1910.00695
% 

%------------- BEGIN CODE --------------

if nargin > 10
    plotFigures = varargin{1};
else
    plotFigures = false;
end

if nargin > 11
    highFidelity = varargin{2};
else
    highFidelity = false;
end

try
    
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
    
    % Time Interval for measurements
    dt = 3600;
    tf = (jdf-jd0)*24*60*60;
%     time = [0:dt:tf]'; m=length(time);
    
    %% Load space weather data
    SWpath = fullfile('Data','SW-All.txt');
    [ SWmatDaily, SWmatMonthlyPred ] = inputSWnrlmsise( SWpath );
    [ SWmatDailyTIEGCM, SWmatMonthlyPredTIEGCM ] = inputSWtiegcm( SWpath );
    [eopdata,SOLdata,DTCdata] = loadJB2008SWdata();
    
    
    %% Load BC estimates
    % Ballistic coefficient data: NORAD ID and BC
    BCfilePath = fullfile('Data','BCdata.txt');
    [BCdata] = loadBCdata( BCfilePath );
    
    rng(1);
    % Extract orbital data and BC data for selected objects
    for i=1:length(selectedObjects)
        % NORAD ID
        ID = selectedObjects(i);
%         % Orbital data (duplicate to ensure data is available)
%         index = find([objects.noradID]==ID);
%         if isempty(index)
%             error('No TLEs found for object %.0f.', ID);
%         end
%         % Get object's TLE data
%         newObjects(i) = objects(index);
        % Ballistic coefficient
%         BCestimates(i) = BCdata(BCdata(:,1)==ID,2);
        % Skysat BC estimate:
        BCestimates(i) = 0.01376 * (1 + randn(1)*0.2);
    end
%     objects = newObjects;
    
    % Number of objects
    nop = length(selectedObjects);
    
    % Cell with object ID text strings
    objectIDlabels = cell(1, nop);
    for i=1:nop
        objectIDlabels(i) = {num2str(selectedObjects(i))};
    end


    useMEE = true;
    
    %% Get GPS measurements
    global gpsDataPath
    [gpsMeas] = getGPSpositionMeasurements(gpsDataPath,selectedObjects,et0,etf);
    measEpochs = gpsMeas(1,:); % Time of measurement
    Meas = gpsMeas(2:end,:); % Measurements
    
    % GPS measurement covariance
    RM = zeros(3,3,size(gpsMeas,2));
    RM(1,1,:) = 0.02^2; % GPS measurement 1-sigma = 20 meters
    RM(2,2,:) = RM(1,1,:);
    RM(3,3,:) = RM(1,1,:);
    
    
    %% Load reduced-order density models
    [AC,BC,Uh,F_U,Dens_Mean,M_U,SLTm,LATm,ALTm,maxAtmAlt,SWinputs,Qrom] = generateROMdensityModel(ROMmodel,r,jd0,jdf);
    
    
    %% Generate initial state guess
    
    % Size of state vector for each object [3xpos,3xvel,1xBC]
    svs = 7;
    
    % Initial state guess: Orbits, BCs and reduced order density
    x0g = zeros(svs*nop+r,1);
    Pv = zeros(svs*nop+r,1); % state variance
    P_GPS = diag( [0.04 0.04 0.04 2e-4 2e-4 2e-4].^2 ); % position 1-sigma = 40m, velocity 1-sigma = 0.2 m/s

    for i = 1:nop
        % Get initial position and velocity from propagated GPS state
        [posVelAtET0] = getPositionVelocityFromGPSmeas(gpsDataPath,selectedObjects(i),BCestimates(i),et0);
%         x0g(svs*(i-1)+1:svs*(i-1)+6,1) = pv2ep( posVelAtET0(1:3), posVelAtET0(4:6), GM_kms );
        
        % Convert assumed pos/vel covariance to MEE
        [meeState,meeCov] = cartCov2meeCov(posVelAtET0, P_GPS, GM_kms);
        x0g(svs*(i-1)+1:svs*(i-1)+6,1) = meeState;
        Pv(svs*(i-1)+1:svs*(i-1)+6) = diag(meeCov);
        
        % Ballistic coefficient guesses
        x0g(svs*i) = BCestimates(i) * 1000;
        % Initial variance for ballistic coefficient
        Pv(svs*(i-1)+7) = (x0g(svs*i) * 0.2)^2; % 1% BC 1-sigma error
    end

    
    % Compute the initial atmosphere state from JB2008 density model
    % Seconds of day in UTC
    UT = hr*3600+mn*60+sc;
    
    % Grid points of ROM model
    sltx = reshape(SLTm,size(SLTm,1)*size(SLTm,2)*size(SLTm,3),1);
    latx = reshape(LATm,size(SLTm,1)*size(SLTm,2)*size(SLTm,3),1);
    altx = reshape(ALTm,size(SLTm,1)*size(SLTm,2)*size(SLTm,3),1);
    
    % Density at grid points according to JB2008 density model
    Dens_JB2008 = zeros(numel(sltx),1);
    for i = 1:numel(sltx)
        % Geographical longitude
        lon = 15*(sltx(i)-UT/3600);
        % Density from JB2008 density model
        Dens_JB2008(i,1) = getDensityJB2008llajd(lon,latx(i),altx(i),jd0) * 1e-9;
    end
    
    % Initialize state for ROM using JB2008
    z0_M = Uh'*(log10(Dens_JB2008)-Dens_Mean); % JB2008 initialization
    
    % Add initial ROM state to initial state guess
    x0g(end-r+1:end,1) = z0_M;
    
    
    %% Process noise Q and initial state covariance P
    
    
    Qv = zeros(svs*nop+r,1); % process variance
    for i = 1:nop
        if useMEE
            
            
%             % Initial variance for orbital state in MEE (equal to measurement noise)
%             Pv(svs*(i-1)+1) = R_TLE(6*(i-1)+1,6*(i-1)+1);
%             Pv(svs*(i-1)+2) = R_TLE(6*(i-1)+2,6*(i-1)+2);
%             Pv(svs*(i-1)+3) = R_TLE(6*(i-1)+3,6*(i-1)+3);
%             Pv(svs*(i-1)+4) = R_TLE(6*(i-1)+4,6*(i-1)+4);
%             Pv(svs*(i-1)+5) = R_TLE(6*(i-1)+5,6*(i-1)+5);
%             Pv(svs*(i-1)+6) = R_TLE(6*(i-1)+6,6*(i-1)+6);
            
            % Process noise for orbital state in MEE
            Qv(svs*(i-1)+1) = 1.5e-8;
            Qv(svs*(i-1)+2) = 2e-14;
            Qv(svs*(i-1)+3) = 2e-14;
            Qv(svs*(i-1)+4) = 1e-14;
            Qv(svs*(i-1)+5) = 1e-14;
            Qv(svs*(i-1)+6) = 1e-12;
            
            % Process noise for orbital state in MEE for 1 hour converted to per second
            % higherProcessNoise
            Qv(svs*(i-1)+1) = 9.6e-19;
            Qv(svs*(i-1)+2) = 1.3e-24;
            Qv(svs*(i-1)+3) = 1.3e-24;
            Qv(svs*(i-1)+4) = 6.4e-25;
            Qv(svs*(i-1)+5) = 6.4e-25;
            Qv(svs*(i-1)+6) = 6.4e-23;
        else
            P_GPS = diag( [0.04 0.04 0.04 2e-4 2e-4 2e-4].^2 ); % position 1-sigma = 40m, velocity 1-sigma = 0.2 m/s
            
            Pv(svs*(i-1)+1:svs*(i-1)+6) = diag(P_GPS);
            
            % Process noise for position and velocity (~ MEE covariance
            % converted to position and velocity)
            Qv(svs*(i-1)+1) = 2e-5;
            Qv(svs*(i-1)+2) = 2e-5;
            Qv(svs*(i-1)+3) = 2e-5;
            Qv(svs*(i-1)+4) = 2e-11;
            Qv(svs*(i-1)+5) = 2e-11;
            Qv(svs*(i-1)+6) = 2e-11;
        end
        
        % Process noise for ballistic coefficient
%         Qv(svs*(i-1)+7) = 1e-16; % 1-sigma error: 1e-8 per 1 hour
        Qv(svs*(i-1)+7) = 2.7e-20; % 1-sigma error: 1e-8 per 1 hour => cov: 1e-16/3600 = 2.7e-20 per second
    end
    
    % Initial variance for reduced-order density state
    Pv(end-r+1:end) = (5e0)*ones(r,1);
    Pv(end-r+1) = 2e1; % First mode
    
    % Process noise for reduced-order density state
    % Use variance of ROM model 1-hour prediction error w.r.t. the training data
    Qv(end-r+1:end) = diag(Qrom) / 3600;
    
    % Initial state covariance and process noise matrices
    P = diag(Pv); % Convert to matrix with variances on the diagonal
    Q = diag(Qv); % Convert to matrix with variances on the diagonal
    
    
    %% Density estimation
    % State estimate
    X_est = x0g; % Initial state guess
    
    % Set state propagation and measurement functions
    if useMEE
        stateFnc = @(xx,t0,tf) propagateState_MeeBcRom(xx,t0,tf,AC,BC,SWinputs,r,nop,svs,F_U,M_U,maxAtmAlt,et0,jd0,highFidelity);
    else
        stateFnc = @(xx,t0,tf) propagateState_PosVelBcRom(xx,t0,tf,AC,BC,SWinputs,r,nop,svs,F_U,M_U,maxAtmAlt,et0,jd0);
    end
    
    % Measurement and residual function
    state2measurementFcn = @(xx,meas) extractSinglePosition(xx,meas,svs,GM_kms); % Function to get position of single object without BC or ROM state
    residualFcn = @(ym,meas) getResidualDefault(ym,meas);
    
    % Run Unscented Kalman filter estimation
    time = measEpochs - et0;
    % Always start at zero
    if time(1) ~= 0
        time = [0, time];
        % Add dummy measurement at t0
        Meas = [zeros(size(Meas,1),1), Meas];
    end
%     [X_est,Pv] = UKFsingleMeasurements(X_est,Meas,time,stateFnc,state2measurementFcn,residualFcn,P,RM,Q,useMEE);
    [X_est,Pv,X_pred] = UKFsingleMeasurements_newProcessNoise(X_est,Meas,time,stateFnc,state2measurementFcn,residualFcn,P,RM,Q,svs,r,useMEE);
    
    global resultsDirPath
    nowTimeStr = datestr(now,'yymmddHHMMSS');
    filenameBase = [resultsDirPath 'ukf_rom_gps_' ROMmodel '_'];
    testCaseName = [sprintf('%04d',yr), sprintf('%02d',mth), sprintf('%02d',dy), '_', num2str(nofDays) 'd_', num2str(nop), 'obj_'];
    
    save('workspace');
    
    save([filenameBase 'workspace_' testCaseName nowTimeStr]);
    
    %% Plot results
    if plotFigures
        
        
        % Plot estimated ROM modes and corresponding uncertainty
        ROMplot = figure;
        for i = 1:r
            subplot(ceil(r/4),4,i)
            plot(time/3600,X_est(end-r+i,:),'Linewidth',1); hold on;
            plot(time/3600,X_est(end-r+i,:) + 3*Pv(end-r+i,:).^0.5,'--k','Linewidth',1);
            plot(time/3600,X_est(end-r+i,:) - 3*Pv(end-r+i,:).^0.5,'--k','Linewidth',1);
            xlabel('Time, hrs');ylabel(sprintf('z_%.0f',i));legend('Estimate','3\sigma');
            title(sprintf('Mode %.0f',i));
            axis tight;
            set(gca,'fontsize', 14);
        end
        savefig(ROMplot,[filenameBase 'ROMmodes_' testCaseName nowTimeStr '.fig']);
        
        % Plot uncertainty in estimated ROM modes
        ROMcovplot = figure;
        for i = 1:r
            subplot(ceil(r/4),4,i)
            plot(time/3600,3*Pv(end-r+i,:).^0.5,'k');
            xlabel('Time, hrs'); ylabel(['z_{' num2str(i) '} 3\sigma']);
            title(sprintf('Mode %.0f',i));
        end
        savefig(ROMcovplot,[filenameBase 'ROMmodesCov_' testCaseName nowTimeStr '.fig']);
        
        % Plot estimated ballistic coefficients and corresponding uncertainty
        BCplot2 = figure;
        for i = 1:nop
            subplot(ceil(nop/2),2,i)
            plot(time/3600,X_est(svs*i,:)/1000,'Linewidth',1); hold on;
            plot(time/3600,X_est(svs*i,:)/1000 + 3*Pv(svs*i,:).^0.5/1000,'--k','Linewidth',1);
            plot(time/3600,X_est(svs*i,:)/1000 - 3*Pv(svs*i,:).^0.5/1000,'--k','Linewidth',1);
            xlabel('Time, hrs');ylabel('BC [m^2/kg]');legend('Estimate','3\sigma');
            title(sprintf('Orbit %.0f, BC=%.2f',selectedObjects(i),X_est(svs*i,end)));
            axis tight;
            set(gca,'fontsize', 14);
        end
        savefig(BCplot2,[filenameBase 'BC_' testCaseName nowTimeStr '.fig']);
        
        % Plot uncertainty in estimated ballistic coefficients
        BCcovplot = figure;
        for i = 1:nop
            subplot(ceil(nop/2),2,i)
            plot(time/3600,3*Pv(svs*i,:).^0.5./X_est(svs*i,:)*100,'k');
            xlabel('Time, hrs');ylabel('BC 3\sigma [%]');
            title(sprintf('Orbit %.0f, BC=%.2f',selectedObjects(i),X_est(svs*i,end)));
        end
        savefig(BCcovplot,[filenameBase 'BCcov_' testCaseName nowTimeStr '.fig']);
        
        % Plot uncertainty in estimated equinoctial orbital elements
        meeCovPlot = figure;
        for i = 1:nop
            for j=1:6
                subplot(2,3,j); hold on;
                plot(time/3600,Pv((i-1)*svs+j,:).^0.5);
            end
        end
        subplot(2,3,1); xlabel('Time [hours]'); ylabel('\sigma_p [km]'); legend(objectIDlabels);
        subplot(2,3,2); xlabel('Time [hours]'); ylabel('\sigma_f [-]');
        subplot(2,3,3); xlabel('Time [hours]'); ylabel('\sigma_g [-]');
        subplot(2,3,4); xlabel('Time [hours]'); ylabel('\sigma_h [-]');
        subplot(2,3,5); xlabel('Time [hours]'); ylabel('\sigma_k [-]');
        subplot(2,3,6); xlabel('Time [hours]'); ylabel('\sigma_L [rad]');
        savefig(meeCovPlot,[filenameBase 'MEEcov_' testCaseName nowTimeStr '.fig']);
        
%         % Plot uncertainty in measurements: equinoctial orbital elements
%         measCovPlot = figure;
%         if useRangeRangeRate
%             for i = 1:nop
%                 objIndices = find(Meas(1,:)==i)-1;
%                 for j=1:2
%                     subplot(1,2,j); hold on;
%                     plot(time(objIndices)/3600,reshape(Rmeas(j,j,objIndices),length(objIndices),1).^0.5);
%                 end
%             end
%             subplot(1,2,1); xlabel('Time [hours]'); ylabel('\sigma_{range} [km]'); legend(objectIDlabels);
%             subplot(1,2,2); xlabel('Time [hours]'); ylabel('\sigma_{rangerate} [km/s]');
%         else
%             for i = 1:nop
%                 objIndices = find(Meas(1,:)==i)-1;
%                 for j=1:6
%                     subplot(2,3,j); hold on;
%                     plot(time(objIndices)/3600,reshape(Rmeas(j,j,objIndices),length(objIndices),1).^0.5);
%                 end
%             end
%             subplot(2,3,1); xlabel('Time [hours]'); ylabel('\sigma_p [km]'); legend(objectIDlabels);
%             subplot(2,3,2); xlabel('Time [hours]'); ylabel('\sigma_f [-]');
%             subplot(2,3,3); xlabel('Time [hours]'); ylabel('\sigma_g [-]');
%             subplot(2,3,4); xlabel('Time [hours]'); ylabel('\sigma_h [-]');
%             subplot(2,3,5); xlabel('Time [hours]'); ylabel('\sigma_k [-]');
%             subplot(2,3,6); xlabel('Time [hours]'); ylabel('\sigma_L [rad]');
%         end
%         savefig(measCovPlot,[filenameBase 'measCov_' testCaseName nowTimeStr '.fig']);
%         
%         if useRangeRangeRate
%             % Plot position errors w.r.t. measurements
%             rangeResidualPlot = figure;
%             for k = 1:nop
%                 ax1(k) = subplot(ceil(nop/2),2,k);
%             end
%             rangeRateResidualPlot = figure;
%             for k = 1:nop
%                 ax2(k) = subplot(ceil(nop/2),2,k);
%             end
%             for k = 1:nop
%                 objIndices = find(Meas(1,:)==k);
% 
%                 rangeRangeRate_res = zeros(2,length(objIndices));
%                 for j=1:length(objIndices)
%                     objI = objIndices(j);
%                     
%                     rangeRangeRate_est = extractSingleRangeRangeRate(X_est(:,objI),Meas(:,objI),svs,GM_kms);
%                     rangeRangeRate_res(:,j) = getResidualRangeRangeRate(rangeRangeRate_est,Meas(:,objI));
%                 end
%                 figure(rangeResidualPlot);
%                 subplot(ax1(k))
%                 plot(time(objIndices)/3600,rangeRangeRate_res(1,:)); hold on;
%                 plot(time(objIndices)/3600,3*reshape(Rmeas(1,1,objIndices-1),length(objIndices),1).^0.5,'-','Color',[0.8 0.8 0.8]);
%                 plot(time(objIndices)/3600,-3*reshape(Rmeas(1,1,objIndices-1),length(objIndices),1).^0.5,'-','Color',[0.8 0.8 0.8]);
%                 xlabel('Time, hrs');ylabel('Range error [km]');
%                 title(sprintf('Orbit %.0f, std= %.3f',objects(k).noradID,std(rangeRangeRate_res(1,:))));
%                 
%                 figure(rangeRateResidualPlot);
%                 subplot(ax2(k))
%                 plot(time(objIndices)/3600,rangeRangeRate_res(2,:)); hold on;
%                 plot(time(objIndices)/3600,3*reshape(Rmeas(2,2,objIndices-1),length(objIndices),1).^0.5,'-','Color',[0.8 0.8 0.8]);
%                 plot(time(objIndices)/3600,-3*reshape(Rmeas(2,2,objIndices-1),length(objIndices),1).^0.5,'-','Color',[0.8 0.8 0.8]);
%                 xlabel('Time, hrs');ylabel('Range rate error [km/s]');
%                 title(sprintf('Orbit %.0f, std= %.4f',objects(k).noradID,std(rangeRangeRate_res(2,:))));
%             end
%             savefig(rangeResidualPlot,[filenameBase 'rangeResiduals_' testCaseName nowTimeStr '.fig']);
%             savefig(rangeRateResidualPlot,[filenameBase 'rangeRateResiduals_' testCaseName nowTimeStr '.fig']);
%         end
        
        if useMEE
            % Plot position errors w.r.t. measurements
            posPlot = figure;
            for k = 1:nop
                objIndices = find(Meas(1,:)==k);

                xx_pv_est = zeros(3,length(objIndices));
                xx_pv_pred = zeros(3,length(objIndices));
                xx_pv_meas = zeros(3,length(objIndices));
                for j=1:length(objIndices)
                    objI = objIndices(j);
                    [pos_est,~] = ep2pv(X_est((k-1)*svs+1:(k-1)*svs+6,objI),GM_kms);
                    xx_pv_est(1:3,j) = pos_est;
                    [pos_pred,~] = ep2pv(X_pred((k-1)*svs+1:(k-1)*svs+6,objI),GM_kms);
                    xx_pv_pred(1:3,j) = pos_pred;
                    xx_pv_meas(1:3,j) = Meas(2:end,objI);
                end
                posErrors_post = sqrt(sum( (xx_pv_est(1:3,:)-xx_pv_meas(1:3,:)) .^2,1));
                posErrors_pre = sqrt(sum( (xx_pv_pred(1:3,:)-xx_pv_meas(1:3,:)) .^2,1));
                subplot(ceil(nop/2),2,k)
                plot(time(objIndices)/3600,posErrors_post); hold on;
                plot(time(objIndices)/3600,posErrors_pre); hold on;
                xlabel('Time, hrs');ylabel('Position error [km]');
                title(sprintf('Orbit %.0f, mean= %.4f',selectedObjects(k),mean(posErrors_post)));
            end
            savefig(posPlot,[filenameBase 'posErr_' testCaseName nowTimeStr '.fig']);

%             % Plot MEE errors w.r.t. measurements
%             meePlot = figure;
%             for k = 1:nop
%                 objIndices = find(Meas(1,:)==k);
%                     xx_mee_est = X_est((k-1)*svs+1:(k-1)*svs+6,objIndices);
%                     xx_mee_meas = Meas(2:end,objIndices);
%                 for j=1:5
%                     subplot(2,3,j)
%                     plot(time(objIndices)/3600,xx_mee_est(j,:)-xx_mee_meas(j,:)); hold on;
%                 end
%                 subplot(2,3,6)
%                 plot(time(objIndices)/3600,wrapToPi(xx_mee_est(6,:)-xx_mee_meas(6,:))); hold on;
%             end
%             subplot(2,3,1); xlabel('Time [hours]'); ylabel('p error [km]'); legend(objectIDlabels);
%             subplot(2,3,2); xlabel('Time [hours]'); ylabel('f error [-]');
%             subplot(2,3,3); xlabel('Time [hours]'); ylabel('g error [-]');
%             subplot(2,3,4); xlabel('Time [hours]'); ylabel('h error [-]');
%             subplot(2,3,5); xlabel('Time [hours]'); ylabel('k error [-]');
%             subplot(2,3,6); xlabel('Time [hours]'); ylabel('L error [rad]');
%             savefig(meePlot,[filenameBase 'MEEerr_' testCaseName nowTimeStr '.fig']);
        end
        
        % Plot estimated density and uncertainty on local solar time v latitude grid
        slt_plot = 0:0.5:24;
        lat_plot = -90:4.5:90;
        [SLT,LAT]=ndgrid(slt_plot,lat_plot);
        sltx = reshape(SLT,length(slt_plot)*length(lat_plot),1);
        latx = reshape(LAT,length(slt_plot)*length(lat_plot),1);
        heights = 500:-100:300; % Plot for different altitudes
        nofHeights = length(heights);
        
        % Plot estimated density
        densPlot = figure;
        for i = 1:nofHeights
            height = heights(i);
            ALT = height*ones(size(SLT,1),size(SLT,2));
            rhoVar = zeros(size(SLT,1),size(SLT,2));
            for j = 1:r
                UhI = F_U{j}(SLT,LAT,ALT);
                rhoVar = rhoVar + UhI*X_est(end-r+j,end);
            end
            MI = M_U(SLT,LAT,ALT);
            rho_rom = 10.^(rhoVar+MI);
            % Plot
            subplot(nofHeights,1,i)
            contourf(SLT,LAT,rho_rom,100,'LineStyle','none');
            max1 = max(max(rho_rom)); min1 = min(min(rho_rom));
            h = colorbar; caxis([min1 max1]);hold on;
            yticks([-90 0 90]); ylabel('Latitude [deg]');
            title(sprintf('Altitude = %.0f km',height));
            ylabel(h,'Density [kg/m^3]','FontSize',14);
            if i < nofHeights; set(gca,'XTickLabel',[]);end
            if i == nofHeights; xlabel('Local solar time'); xticks([0 6 12 18 24]); end
            set(gca,'FontSize',14);
        end
        savefig(densPlot,[filenameBase 'densAlt_' testCaseName nowTimeStr '.fig']);
        
        % Plot uncertainty in estimated density
        densCovPlot = figure;
        set(gcf,'Color','w');
        for i = 1:nofHeights
            height = heights(i);
            H_SL = zeros(numel(sltx),r);
            for ij = 1:numel(sltx)
                for jk = 1:r
                    H_SL(ij,jk) = F_U{1,jk}(sltx(ij),latx(ij),height); % ROM to grid transformation matrix
                end
            end
            Pyy = diag(H_SL(:,:) * diag(Pv(end-r+1:end,end)) * H_SL(:,:)'); % ROM covariance mapped to local solar time and latitude
            Pyyr = reshape(Pyy,length(slt_plot),length(lat_plot)); % Convert to matrix
            Pyy1 = 100*Pyyr.^0.5*log(10); % Convert covariance of log density to standard deviation of density
            
            % Plot estimated standard deviation (1-sigma error) of density
            subplot(nofHeights,1,i)
            contourf(slt_plot,lat_plot,Pyy1',100,'LineStyle','none');
            max1 = max(max(Pyy1)); min1 = min(min(Pyy1));
            h = colorbar; caxis([min1 max1]);hold on;
            yticks([-90 0 90]); ylabel('Latitude [deg]');
            title(sprintf('Altitude = %.0f km',height));
            ylabel(h,'1\sigma error [%]','FontSize',14);
            if i < nofHeights; set(gca,'XTickLabel',[]);end
            if i == nofHeights; xlabel('Local solar time'); xticks([0 6 12 18 24]); end
            set(gca,'FontSize',14);
        end
        savefig(densCovPlot,[filenameBase 'densCovAlt_' testCaseName nowTimeStr '.fig']);
        
    end
    
catch errMsg
    % Catch error
    rethrow(errMsg);
end

%------------- END OF CODE --------------