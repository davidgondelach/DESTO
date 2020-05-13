%--------------------------------------------------------------------------
%
% AccelSolrad: Computes the acceleration due to solar radiation pressure
%			   assuming the spacecraft surface normal to the Sun direction
%
% Inputs:
%   r           Spacecraft position vector 
%   r_Earth	    Earth position vector (Barycentric)
%   r_Moon		Moon position vector (geocentric)
%   r_Sun       Sun position vector (geocentric)
%   r_SunSSB    Sun position vector (Barycentric)
%   Area        Cross-section 
%   mass        Spacecraft mass
%   Cr          Solar radiation pressure coefficient
%   P0          Solar radiation pressure at 1 AU 
%   AU          Length of one Astronomical Unit
%   shm         Shadow model (geometrical or cylindrical)
%
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Notes:
%   r, r_sun, Area, mass, P0 and AU must be given in consistent units,
%   e.g. m, m^2, kg and N/m^2. 
%
% Last modified:   2018/01/27   M. Mahooti
% 
%--------------------------------------------------------------------------
function a = AccelSolrad(r,r_Earth,r_Moon,r_Sun,r_SunSSB,AoM,Cr,P0,AU,shm)
%        Moon wrt Earth          pccor      rpc
%        Earth wrt Sun           ccor       rc
%        Moon wrt Sun            pscor      rps   
%        Satellite wrt Earth     sbcor      rsb  
%        Satellite wrt Sun       bcor       rb 
%        Satellite wrt Moon      sbpcor     rsbp
pccor = r_Moon;
ccor = r_Earth-r_SunSSB;
pscor = r_Moon-r_Sun;
sbcor = r;
bcor = r-r_Sun;
sbpcor = r-r_Moon;
if ( strcmp(shm,'cylindrical') )
    nu = Cylindrical(r,r_Sun);
else
    error('geometrical shadow not implemented');
    [nu,~] = Shadow(pccor,ccor,pscor,sbcor,bcor,sbpcor);
end
% Acceleration
a = nu*Cr*AoM*P0*(AU*AU)*bcor/(norm(bcor)^3);

% END