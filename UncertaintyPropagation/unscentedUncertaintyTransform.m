function [Xout,Pout] = unscentedUncertaintyTransform(Xin,transformationFnc,Pin,Q)
%unscentedUncertaintyTransform
%   X_est: initial state guess
%   Meas: measurements
%   time: measurement times
%   stateFnc: function to propagate state
%   measurementFnc: function to convert from state to measurement space
%   P: state covariance matrix
%   RM: observation noise
%   Q: process noise
%
% Copyright (C) 2021 by David Gondelach
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Aug 2020; Last revision: 31-Aug-2020

%------------- BEGIN CODE --------------

% Unscented transformation
% Compute the Sigma Points
kappa = 0;
[Wm,Wc,L,lam] = Unscented_Transform(Xin,kappa);
SR_Wc = sqrt(Wc); %SR_Wm = sqrt(Wm);

S=chol(Pin)';
SR_Q = sqrt(Q); % process noise

% Sigma points
eta = sqrt(L+lam);
sigv=real([eta*S -eta*S]);
XsigmaIn=[Xin sigv+kron(Xin,ones(1,2*L))];

% Transformation
[XsigmaOut] = transformationFnc(XsigmaIn);

% State update
Xout = Wm(1) * XsigmaOut(:,1) + Wm(2) * sum(XsigmaOut(:,2:end),2);

% Covariance update
[~,S_minus]=qr([(SR_Wc(2)*(XsigmaOut(:,2:end)-kron(Xout,ones(1,2*L)))) SR_Q]',0);
S_minus=cholupdate(S_minus,Wc(1)*(XsigmaOut(:,1)-Xout))';
S = S_minus;
Pout = (S*S')';

end

%------------- END OF CODE --------------