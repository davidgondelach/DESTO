function [Wm,Wc,L,lambda] = Unscented_Transform(x0f,varargin)
%Unscented_Transform - Compute weights for unscented transformation
%
% Copyright (C) 2021
%
% This code is licensed under the GNU General Public License version 3.
%
% Reference: Wan, E. A., & Van Der Merwe, R. (2001). The unscented Kalman filter, In: Kalman filtering and neural networks, pp. 221-280.
%

%------------- BEGIN CODE --------------

L=length(x0f);
alpha = 1;
beta = 2;
kappa = 3 - L;
if nargin > 1
    kappa = varargin{1};
end
lambda = alpha^2*(L + kappa) - L;

W0m = lambda/(L + lambda);
W0c = lambda/(L + lambda)+(1-alpha^2+beta);
Wim = 1/(2*(L + lambda));

Wm = [W0m Wim+zeros(1,2*L)];
Wc = [W0c Wim+zeros(1,2*L)];

end

%------------- END OF CODE --------------
