%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%  Thermospheric Density Estimation Via Two-Line-Element Data Assimilation  %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Copyright (C) 2020 by David Gondelach
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%  Author: David Gondelach
%  Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
%  email: davidgondelach@gmail.com
%  Jan 2020; Last revision: 31-Jan-2020
%
%  Reference:
%  D.J. Gondelach and R. Linares, "Real-Time Thermospheric Density
%  Estimation Via Two-Line-Element Data Assimilation", Space Weather, 2020
%  https://doi.org/10.1029/2019SW002356 or https://arxiv.org/abs/1910.00695
% 


%------------- BEGIN CODE --------------

clearvars;
clearvars -global;


%% SETTINGS
% Specify the date, reduced-order model, reduction order and objects to be
% used to estimation here.

% Estimation window
% Continuous data from:2020-01-03T06:03:06.870935 , till:2020-01-28T06:10:43.642242
yr      = 2020; % Year
mth     = 5;    % Month
dy      = 1;    % Day
hr      = 0;
mn      = 0;
sec     = 0;
nofDays = 1;   % Number of days

% Use high fidelity dynamical model
highFidelity = true;

% Reduced-order model
ROMmodel = 'NRLMSISE_1997_2008';  % Name of reduced-order density model: JB2008_1999_2010, NRLMSISE_1997_2008 or TIEGCM_1997_2008
r  = 10;                        % Reduced order

% NORAD catalog IDs of objects used for estimation
% Default: 17 objects: [63;165;614;2153;2622;4221;6073;7337;8744;12138;12388;14483;20774;23278;27391;27392;26405]
% selectedObjects = [63;165;614;2153;2622;4221;6073;7337;8744;12138;12388;14483;20774;23278;27391;27392;26405]; % TLE
% selectedObjects = [614;2153;2622;4221;12138]; % Radar
% selectedObjects = [614;2153;2622;4221;12138;750;2016;2389;6073;7337;8744;12388;14483;20774;23278]; % Radar
selectedObjects = [39418,41601,41771,41772,41773,41774,42987,42988,42989,42990,42991,42992,43797,43802]; % GPS: Planet Skysats
selectedObjects = [41601,41771,41772,41773,43802]; % GPS: Planet Skysats
selectedObjects = sortrows(selectedObjects);

% Display date
datetime(yr,mth,dy)

%% SET PATHS
% *** SPECIFY YOUR SPICE TOOLBOX DIRECTORY HERE! ***
% spicePath = fullfile('[SPICE TOOLKIT DIRECTORY]','mice'); 
spicePath = fullfile('/Users/davidgondelach/Documents','mice'); 
global resultsDirPath gpsDataPath
resultsDirPath = ['/Users/davidgondelach/Google Drive/PostDoc/DensityEstimation/GPSmeas/',ROMmodel,'/'];
gpsDataPath = '/Users/davidgondelach/Documents/PostDoc/GPSdata';

addpath( 'AstroFunctions' );
addpath( 'Estimation' );
addpath( 'JB2008' );
addpath( 'ROMDensityModels' );
addpath( 'SpaceWeather' );
% addpath( 'TLEdata' );
addpath( 'GPSdata' );
addpath( 'UncertaintyPropagation' );
addpath( fullfile(spicePath,'src','mice') );
addpath( fullfile(spicePath,'lib') );


%% LOAD KERNELS, GRAVITY MODEL, EARTH ORIENTATION PARAMETERS AND SGP4
% Load SPICE kernels and ephemerides
kernelpath  = fullfile('Data','kernel.txt');
loadSPICE(kernelpath);

% Load gravity model
if highFidelity
    gravmodeldegree  = 48;  % Use degree and order 48 for the spherical harmonics
else
    gravmodeldegree  = 20;  % Use degree and order 20 for the spherical harmonics
end
loadGravityModel( gravmodeldegree );

% Load Earth orientation parameters (needed for TEME to ECI transformation)
global EOPMat
EOPpath = fullfile('Data','EOP-All.txt');
[ EOPMat ] = inputEOP_Celestrak( EOPpath );

% Setup the SGP4 propagator.
loadSGP4();


%% PERFORM DENSITY ESTIMATION
plotFigures = true;

runDensityEstimationGPS(yr,mth,dy,hr,mn,sec,nofDays,ROMmodel,r,selectedObjects,plotFigures,highFidelity);


%% Clear memory
% Clear cspice memory
cspice_kclear;

%------------- END OF CODE --------------