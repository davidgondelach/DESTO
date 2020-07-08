function [rho] = getDensityJB2008(pos,et)

% Date and time
jdatestr    = cspice_et2utc( et, 'J', 12 );
jdate       = str2double(jdatestr(4:end)); % Cut trailing 'JD ' off from string
[yy, mm, dd, hh, mnmn, ss] = datevec(jdate-1721058.5);

% Longitude, Latitude and Height of Satellite and Density
% TODO: USE ECEF POSITION INSTEAD
% TODO: USE GEOCENTRIC LATITUDE
[lon,lat,alt]=gc2gd(pos,yy,mm,dd,hh,mnmn,ss,0,0,0);
lon(lon>180) = lon(lon>180) - 360; % Geocentric longitude

[rho] = getDensityJB2008llajd(lon,lat,alt,jdate);

end

