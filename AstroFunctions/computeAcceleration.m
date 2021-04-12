function [ f ] = computeAcceleration(t,x,et0,settings,densityModelData)
% t         current time
% x         current state of multiple objects (incl BC)
% et0       initial ephemerides time

% Date and time
et = et0 + t;

% Dynamical model
x = reshape(x,7,[]);
[m,n] = size(x);

% Ballistic coefficient
b_star = x(7,:) * 1e-6; % [km^2/kg]

% Position and velocity in ECEF ref frame
xform = cspice_sxform('J2000', 'ITRF93', et ); % J2000 to ECEF transformation matrix
x_ecef = xform * x(1:6,:);
rr_ecef = x_ecef(1:3,:);
vv_ecef = x_ecef(4:6,:);

% If high fidelity, add Sun, Moon and SRP perturbations
if settings.thirdbody == 1
    
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
    
end

% Earth gravitational acceleration
[aa_grav_ecef_x, aa_grav_ecef_y, aa_grav_ecef_z] = computeEarthGravitationalAcceleration(rr_ecef'*1000);
% Gravitational accelerations in ECI [km/s^2]
aa_eci = xform(1:3,1:3)' * [aa_grav_ecef_x'; aa_grav_ecef_y'; aa_grav_ecef_z'] / 1000;

for i=1:n
    rr_sat = x(1:3,i);
    
    % Drag acceleration
    if settings.drag ~= 0
        % Velocity w.r.t. Earth atmosphere
        mag_v_ecef = sqrt( sum( vv_ecef(:,i).^2, 1 ));
        
        % Atmospheric density
        if settings.drag == 1
            % Use reduced-order density model
            [rho] = getDensityROM(rr_sat',et,densityModelData.romStateTime,densityModelData.r,densityModelData.F_U,densityModelData.M_U);
        elseif settings.drag == 3
            % Use JB2008 atmosphere model
            [rho] = getDensityJB2008(rr_sat',et);
        end
        
        % Drag acceleration
        accDrag_ecef = [- 1/2.*b_star(1,i).*rho.*mag_v_ecef.*vv_ecef(1,i);
                        - 1/2.*b_star(1,i).*rho.*mag_v_ecef.*vv_ecef(2,i);
                        - 1/2.*b_star(1,i).*rho.*mag_v_ecef.*vv_ecef(3,i) ];
        aa_eci(:,i) = aa_eci(:,i) + xform(1:3,1:3)' * accDrag_ecef;
    end
    
    if settings.thirdbody == 1
        
        % Solar radiation pressure
        AoMSRP = b_star(1,i)/2.2; % Approximate area-to-mass ratio (assuming BC=Cd*AoM and Cd=2.2)
        % SRP acceleration
        aa_SRP = AccelSolrad(rr_sat, [0;0;0],[0;0;0], rr_Sun, [0;0;0] ,AoMSRP,C_R,P_sun,AU,'cylindrical');
        aa_eci(:,i) = aa_eci(:,i) + aa_SRP;
        
        % Moon gravitational acceleration
        aa_Moon = AccelPointMass(rr_sat,rr_Moon,GM_Moon);
        aa_eci(:,i) = aa_eci(:,i) + aa_Moon;
        
        % Sun gravitational acceleration
        aa_Sun = AccelPointMass(rr_sat,rr_Sun,GM_Sun);
        aa_eci(:,i) = aa_eci(:,i) + aa_Sun;
    end
end

% State derivative
f        = zeros(m,n);
f(1:3,:) = x(4:6,:); % Velocity
f(4:6,:) = aa_eci; % Acceleration
f(7,:)   = 0; % Constant bstar

f = reshape(f,[],1);

end


