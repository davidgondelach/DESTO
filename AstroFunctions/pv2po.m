% function to transform cartesian coordinates in a vector of classical
% orbital parameters
function PO = pv2po(rr, vv, mu)
% input: rr position vecotor in units consistent with mu 
%        vv velocity vector in units consistent with mu  
%        mu gravity parameter in selcted units 
% output: 
%        PO vector of classical orbital parameters angles in rad and a
%        in lenght unites as rr and mu

r = norm(rr);
v = norm(vv);
a = mu/(2*(mu/r-v^2/2));
h = cross(rr,vv);

if (h(1)^2+h(2)^2)==0
   OMEGA=0;
else  
    
sinOMEGA = h(1)/sqrt(h(1)^2+h(2)^2);
cosOMEGA = -h(2)/sqrt(h(1)^2+h(2)^2);

if cosOMEGA>=0
   if sinOMEGA>=0
      OMEGA = asin(h(1)/sqrt(h(1)^2+h(2)^2));
   else
      OMEGA = 2*pi+asin(h(1)/sqrt(h(1)^2+h(2)^2));
   end 
else
   if sinOMEGA>=0
      OMEGA = acos(-h(2)/sqrt(h(1)^2+h(2)^2));
   else
      OMEGA = 2*pi-acos(-h(2)/sqrt(h(1)^2+h(2)^2));
   end
end
end

OMEGA = real(OMEGA);

ee = 1/mu*(cross(vv,h))-rr/norm(rr);
e = norm(ee);

i = acos(h(3)/norm(h));

if e<=1e-8 && i<1e-8
   e = 0;
   omega = atan2(rr(2),rr(1));
   theta = 0;
   PO = [a, e, i, OMEGA, omega, theta];
   return
end

if e<=1e-8 && i>=1e-8
    omega=0;
    K = [cos(omega)*cos(OMEGA)-sin(omega)*sin(i)*sin(OMEGA) cos(omega)*sin(OMEGA)+sin(omega)*cos(i)*cos(OMEGA) sin(omega)*sin(i);
         -sin(omega)*cos(OMEGA)-cos(omega)*cos(i)*sin(OMEGA) -sin(omega)*sin(OMEGA)+cos(omega)*cos(i)*cos(OMEGA) cos(omega)*sin(i);
         sin(OMEGA)*sin(i) -cos(OMEGA)*sin(i) cos(i)];
    rr = K*rr;
    theta = atan2(rr(2),rr(1));
    PO = [a, e, i, OMEGA, omega, theta];
    return
end

theta = acos((rr'*ee)/(norm(rr)*norm(ee)));
if rr'*vv<0
   theta=2*pi-theta;
end
theta = real(theta);


if i<=1e-8 && e>=1e-8
   i=0;
   omega = atan2(ee(2),ee(1));
   PO = [a, e, i, OMEGA, omega, theta];
   return
end

sino = rr(3)/r/sin(i);
coso = (rr(1)*cos(OMEGA)+rr(2)*sin(OMEGA))/r;

if coso>=0
   if sino>=0
      o = asin(rr(3)/r/sin(i));
   else
      o = 2*pi+asin(rr(3)/r/sin(i));
   end
   else
   if sino>=0
      o = acos((rr(1)*cos(OMEGA)+rr(2)*sin(OMEGA))/r);
   else
      o = 2*pi-acos((rr(1)*cos(OMEGA)+rr(2)*sin(OMEGA))/r);
   end
end
o = real(o);
omega = o-theta;

if omega<0
   omega = omega+2*pi;
end
omega = real(omega);

PO = [a, e, i, OMEGA, omega, theta];