function [xf_pv] = propagateState_PosVelBcRom(x0,t0,tf,AC,BC,SWinputs,r,nop,svs,F_U,M_U,maxAtmAlt,et0,jdate0)
%PROPAGATESTATE_MEEBCROM - Propagate objects and ROM density
% Convert state in modified equinoctial elements to Cartesian coordinates
% propagate states and reduced-order density and convert Cartesian states
% back to modified equinoctial elements.
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and
% Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020

%------------- BEGIN CODE --------------

opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events', @(t,x) isdecayed(t,x,nop*size(x0,2),svs));
[~,xf_out]=ode113(@(t,x) computeDerivative_PosVelBcRom(t,x,AC,BC,SWinputs,r,nop,svs,F_U,M_U,maxAtmAlt,et0,jdate0),[t0 tf],x0,opts);
xf_pv = reshape(xf_out(end,:)',nop*svs+r,[]);

end

%------------- END OF CODE --------------
