function [ f ] = computeDerivative_PosVelBcRom(t,xp,AC,BC,SWinputs,r,noo,svs,F_U,M_U,maxAtmAlt,et0,jdate0,highFidelity)
% COMPUTEDERIVATIVE_POSVELBCROM - Compute derivatives of objects
% position, velocity and BC, and reduced-order state
%
% Syntax:  [ f ] = computeDerivative_PosVelBcRom(t,xp,AC,BC,Inp,rR,noo,svs,F_U,M_U,maxAtmAlt,et0,jdate0)
%
% Inputs:
%   t           current time: seconds since et0 [s]
%   xp          state vector: position and velocity (J2000) and BC of
%               multiple objects and reduced order density state
%   AC          continuous-time state transition matrix for the reduced
%               order density state
%   BC          continuous-time input matrix for the reduced order density
%               state dynamics
%   SWinputs    Space weather inputs
%   r           number of reduced order modes [integer]
%   noo         number of objects [integer]
%   svs         state size per object [integer]
%   F_U         interpolant of gridded reduced-order modes
%   M_U         interpolant of gridded mean density
%   maxAtmAlt   maximum altitude of ROM density model
%   et0         initial ephemeris time (seconds since J2000 epoch)
%   jdate0      initial Julian date
%
% Outputs:
%    f          Time derivative of xp: dxp/dt
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and
% Astronautics
% email: davidgondelach@gmail.com
% Jan 2020; Last revision: 31-Jan-2020

%------------- BEGIN CODE --------------

% Convert state from single column to multi-column matrix
x = reshape(xp,svs*noo+r,[]);

% Date and time
et = et0 + t; % Ephemeris time
jdate = jdate0 + t / 86400; % Julian date

% Space weather inputs for current time
SWinputs = interp1(SWinputs(1,:),SWinputs(2:end,:)',jdate)';

% m = state vector length = noo*svs + r
% n = number of states (i.e. number of sigma points in UKF)
[m,n]=size(x);

% State derivative f=dx/dt
f=zeros(m,n);

% J2000 to ECEF transformation matrix
xform = cspice_sxform('J2000', 'ITRF93', et ); % J2000 to ECEF transformation matrix

% Object states in ECI
x_eci = reshape(x(1:end-r,:),svs,[]);
% Object states in ECEF
x_ecef = xform*x_eci(1:6,:); % State in ECEF
rr_ecef = x_ecef(1:3,:); % Position in ECEF
vv_ecef = x_ecef(4:6,:); % Velocity in ECEF
mag_v_ecef = sqrt( sum( vv_ecef.^2, 1 )); % Magnitude of velocity in ECEF

% Gravitational accelerations in ECEF [m/s^2]
nofStates = size(rr_ecef,2);
partSize = 600;
nofParts = floor( size(rr_ecef,2) / partSize );
aa_grav_ecef_x = zeros(nofStates,1);
aa_grav_ecef_y = zeros(nofStates,1);
aa_grav_ecef_z = zeros(nofStates,1);
for i = 1:nofParts
    [aa_grav_ecef_x((i-1)*partSize+1:i*partSize,1), aa_grav_ecef_y((i-1)*partSize+1:i*partSize,1), aa_grav_ecef_z((i-1)*partSize+1:i*partSize,1)] = ...
        computeEarthGravitationalAcceleration(rr_ecef(:,(i-1)*partSize+1:i*partSize)'*1000);
end
[aa_grav_ecef_x(nofParts*partSize+1:end,1), aa_grav_ecef_y(nofParts*partSize+1:end,1), aa_grav_ecef_z(nofParts*partSize+1:end,1)] = ...
        computeEarthGravitationalAcceleration(rr_ecef(:,nofParts*partSize+1:end)'*1000);
% [aa_grav_ecef_x, aa_grav_ecef_y, aa_grav_ecef_z] = computeEarthGravitationalAcceleration(rr_ecef'*1000);
% Gravitational accelerations in ECI [km/s^2]
aa_grav_eci = xform(1:3,1:3)' * [aa_grav_ecef_x'; aa_grav_ecef_y'; aa_grav_ecef_z'] / 1000;

% Reduced order density state
romState = x(end-r+1:end,:);

% Compute accelerations per object
rho = zeros(noo,n);
b_star = zeros(noo,n);
for i = 1:noo
    
    % Atmospheric densities [kg/m^3]
    % Position in J2000
    rr_eci = x(svs*(i-1)+1:svs*(i-1)+3,:);
    % Compute density using reduced-order density model [kg/m^3]
    rho(i,:) = getDensityROM(rr_eci,jdate,romState,r,F_U,M_U,maxAtmAlt);
    
    % Ballistic coefficients (BC) [m^2/(1000kg)]
    b_star(i,:) = x(svs*i,:);
    
end

% BC * density * velocity
BCrhoV = reshape(b_star.*rho,1,[]) .* mag_v_ecef;
% Drag accelerations in ECEF [km/s^2]
aa_drag_ecef(1,:) = - 1/2.*BCrhoV.*vv_ecef(1,:); % ECEF x-direction
aa_drag_ecef(2,:) = - 1/2.*BCrhoV.*vv_ecef(2,:); % ECEF y-direction
aa_drag_ecef(3,:) = - 1/2.*BCrhoV.*vv_ecef(3,:); % ECEF z-direction
% Drag accelerations in ECI [km/s^2]
aa_drag_eci = xform(1:3,1:3)' * aa_drag_ecef;

% Total accelerations in ECI [km/s^2]
aa_grav_drag_eci = reshape(aa_grav_eci + aa_drag_eci, 3*noo,n);

% Time derivatives of position and velocity due to velocity and gravity and drag accelerations
for i = 1:noo
    % Velocities in J2000 frame [km/s]
    f(svs*(i-1)+1:svs*(i-1)+3,:) = x(svs*(i-1)+4:svs*(i-1)+6,:); % J2000 x-direction
    
    % Total accelerations in ECI [km/s^2]
    f(svs*(i-1)+4:svs*(i-1)+6,:) = aa_grav_drag_eci(3*(i-1)+1:3*i,:);
end

% Time derivative of ballistic coefficients is zero
f(7:svs:end-r,:) = 0;


% If high fidelity, add Sun, Moon and SRP perturbations
if highFidelity
    
    %%% Compute Sun Moon %%%
    % Sun position in J2000 ref frame
    rr_Sun = cspice_spkezr('Sun',et,'J2000','NONE', 'Earth');
    rr_Sun = rr_Sun(1:3,1);
    gravconst = 6.67259e-20; % [km^3/kg/s^2]
    GM_Sun    = gravconst*1.9891e30; % [kg]
    % Moon position in J2000 ref frame
    rr_Moon = cspice_spkezr('Moon',et,'J2000','NONE', 'Earth');
    rr_Moon = rr_Moon(1:3,1);
    GM_Moon = gravconst*7.3477e22; % [kg]
    
    % Solar radiation
    AU          = 1.49597870691e8; % [km]
    L_sun       = 3.846e20; % [MW] Luminosity of the Sun -> 3.846e26 W -> W = kg m/s2 m/s -(to km)-> 1e-6 kg km/s2 km/s -> 1e-6 MW/W
    c_light     = 299792.458; % [km/s] speed of light
    P_sun       = L_sun/(4*pi*AU^2*c_light);
    C_R = 1.2;
    
    %%% SRP and lunisolar perturbations %%%
    for i = 1:noo
        aa_eci = zeros(3,n);
        for j=1:n
            rr_sat = x(svs*(i-1)+1:svs*(i-1)+3,j);
            
            % Solar radiation pressure
            b_star = x(svs*i,j);
            AoMSRP = b_star/2.2 * 1e-9; % Approximate area-to-mass ratio (assuming BC=Cd*AoM and Cd=2.2)
            % SRP acceleration
            aa_SRP = AccelSolrad(rr_sat, [0;0;0],[0;0;0], rr_Sun, [0;0;0] ,AoMSRP,C_R,P_sun,AU,'cylindrical');
            aa_eci(:,j) = aa_SRP;
            
            % Moon gravitational acceleration
            aa_Moon = AccelPointMass(rr_sat,rr_Moon,GM_Moon);
            aa_eci(:,j) = aa_eci(:,j) + aa_Moon;
            
            % Sun gravitational acceleration
            aa_Sun = AccelPointMass(rr_sat,rr_Sun,GM_Sun);
            aa_eci(:,j) = aa_eci(:,j) + aa_Sun;
        end
        % Add SRP and lunisolar accelerations
        f(svs*(i-1)+4:svs*(i-1)+6,:) = f(svs*(i-1)+4:svs*(i-1)+6,:) + aa_eci;
    end
end

% Time derivative of reduced-order density state: dz/dt = Ac*z + Bc*u
f(end-r+1:end,:) = AC * romState + BC * SWinputs;

% Convert state derivative to single column
f = reshape(f,[],1);

end

%------------- END OF CODE --------------
