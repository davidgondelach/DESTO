function [X_est,Pv,X_pred] = UKFsingleMeasurements_newProcessNoise(X_est,Meas,time,statePredictionFnc,state2measurementFcn,residualFcn,P,RM,Q,svs,r,useMEE)
%UKF Square-root Unscented Kalman filter
%   X_est: initial state guess
%   Meas: measurements
%   time: measurement times
%   stateFnc: function to propagate state
%   measurementFnc: function to convert from state to measurement space
%   P: state covariance matrix
%   RM: observation noise
%   Q: process noise
%
% This code is licensed under the GNU General Public License version 3.
%
% Based on code by P.M. Mehta, University of Minnesota, 2018
%
% Modified by: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020
%
% Reference: Wan, E. A., & Van Der Merwe, R. (2001). The unscented Kalman filter, In: Kalman filtering and neural networks, pp. 221–280.
%

%------------- BEGIN CODE --------------

% Unscented Filter Parameter
% Compute the Sigma Points
[Wm,Wc,L,lam] = Unscented_Transform(X_est);
SR_Wc = sqrt(Wc); SR_Wm = sqrt(Wm);

S=chol(P)';
eta = sqrt(L+lam);

X_pred = X_est;
% Add initial state variance to state variance history
Pv(:,1) = diag(P);

try
    
    m = size(Meas,2);
    for i = 1:m-1
        
        fprintf('%.0f of %.0f \n',i,m-1);
        
        sigv = real([eta*S -eta*S]);
        xx = [X_est(:,i) sigv+kron(X_est(:,i),ones(1,2*L))];
        
        % Time Update
        [Xp] = statePredictionFnc(xx,time(i),time(i+1));
        
        X_est(:,i+1) = Wm(1) * Xp(:,1) + Wm(2) * sum(Xp(:,2:end),2);
        X_pred(:,i+1) = X_est(:,i+1);
        
        % Process noise (see
        % http://www.robots.ox.ac.uk/~ian/Teaching/Estimation/LectureNotes2.pdf)
        % Q is process white noise variance
        deltat = time(i+1)-time(i); % prop time in seconds
        Q_dt = zeros(length(Q));
        % The orbital state process noise grows cubically with time: Q(deltat) = Q * 1/3 * deltat^3
        Q_dt(1:end-r,1:end-r) = Q(1:end-r,1:end-r) * 1/3*deltat^3; % Discrete time process noise covariance matrix for deltat
        % The BC process noise grows linearly with time: Q(deltat) = Q * deltat
        Q_dt(svs:svs:end-r,svs:svs:end-r) = Q(svs:svs:end-r,svs:svs:end-r) * deltat; % Discrete time process noise covariance matrix for deltat
        % The ROM density process noise grows linearly with time: Q(deltat) = Q * deltat
        Q_dt(end-r+1:end,end-r+1:end) = Q(end-r+1:end,end-r+1:end) * deltat; % Discrete time process noise covariance matrix for deltat
        SR_Q = sqrt(Q_dt); % process noise
        
        % Get Propagated Square Root
        [~,S_minus] = qr([(SR_Wc(2)*(Xp(:,2:end)-kron(X_est(:,i+1),ones(1,2*L)))) SR_Q]',0);
        S_minus = cholupdate(S_minus,Wc(1)*(Xp(:,1)-X_est(:,i+1)))';
        
        % Measurement function
        [Ym] = state2measurementFcn(Xp,Meas(:,i+1)); % [nofMeas x nofSigma]
        ym = Wm(1) * Ym(:,1) + Wm(2) * sum(Ym(:,2:end),2); % [nofMeas x 1]
        
        DY = Ym(:,1)-ym; % [nofMeas x 1]
        DY2 = Ym(:,2:end)-kron(ym,ones(1,2*L)); % [nofMeas x nofSigma-1]
        if useMEE
            DY(6:6:end) = wrapToPi(DY(6:6:end)); % Wrap difference in true longitude to [-pi,pi]
            DY2(6:6:end) = wrapToPi(DY2(6:6:end)); % Wrap difference in true longitude to [-pi,pi]
        end
        
        % Measurement Update
        SR_R = sqrt(diag(diag(RM(:,:,i)))); % measurement noise
        [~,S_y] = qr([SR_Wc(2)*(DY2) SR_R]',0); % [nofMeas x nofMeas]
        S_y = cholupdate(S_y,Wc(1)*DY)'; % [nofMeas x nofMeas]
        
        % Calculate Pxy
        Pxy0 = Wc(1)*(Xp(:,1)-X_est(:,i+1))*DY'; % [nofStates x nofMeas]
        Pyymat = DY2;
        Pmat = Xp(:,2:end)-kron(X_est(:,i+1),ones(1,2*L)); % [nofStates x nofSigma-1]
        Pxy = Pxy0+Wc(2)*(Pmat*Pyymat'); % [nofStates x nofMeas]
        
        % Measurement residual
        yres = residualFcn(ym,Meas(:,i+1)); % [nofMeas x 1]
        
        % Gain and Update
        KG = real(Pxy/S_y')/S_y; % [nofStates x nofMeas]
        X_est(:,i+1) = X_est(:,i+1) + KG * yres;
        U = KG * S_y;
        S = S_minus; % [nofStates x nofStates]
        for j = 1:length(ym)
            S = cholupdate(S',U(:,j),'-')';
        end
        
        HH(:,:,i+1) = (pinv(S*S')*KG*diag(diag(RM(:,:,i))))';
        Pv(:,i+1) = diag(S*S')';
        
    end
    
catch errMsg
    nowTimeStr = datestr(now,'yymmddHHMMSS');
    save(['workspace_UKF_' nowTimeStr]);
    rethrow(errMsg);
end

end

%------------- END OF CODE --------------
