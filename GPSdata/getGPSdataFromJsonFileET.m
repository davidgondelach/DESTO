function [gpsData] = getGPSdataFromJsonFileET(GPSdataPath,objectPlanetIDstr,objectNORADID,et0,etf)
%getGPSdataFromJsonFileET - Load GPS data for given object from json files 
% for given date range.
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

% See: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/timout_c.html
pic = 'DD-Mon-YYYY HR:MN:SC.##### ::RND'; % 'DD-Mon-YYYY HR:MN:SC.##### ::RND'

% Start date
date0str = cspice_timout( et0, pic );
date0 = datetime(date0str);
year0 = date0.Year;
month0 = date0.Month;
day0 = date0.Day;

% End date
datefstr = cspice_timout( etf, pic );
datef = datetime(datefstr);
yearf = datef.Year;
monthf = datef.Month;
dayf = datef.Day;

% Load GPS data
[gpsData] = getGPSdataFromJsonFile(GPSdataPath,objectPlanetIDstr,objectNORADID,year0,month0,day0,yearf,monthf,dayf);

% Remove GPS data outside observation window
firstGPSindex   = find([gpsData(:).tET]>=et0,1);
lastGPSindex    = find([gpsData(:).tET]<=etf,1,'last');
gpsData         = gpsData(firstGPSindex:lastGPSindex);
    
end

%------------- END OF CODE --------------