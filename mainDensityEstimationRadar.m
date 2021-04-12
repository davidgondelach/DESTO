%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%  Thermospheric Density Estimation Via Radar Tracking Data Assimilation    %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Copyright (C) 2021 by David Gondelach
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
%  Jan 2020; Last revision: 01-Mar-2021
%
%  Reference:
%  D.J. Gondelach and R. Linares, "Real-Time Thermospheric Density 
%  Estimation Via Radar And GPS Tracking Data Assimilation", Space Weather, 2021
%  https://doi.org/10.1029/2020SW002620
% 


%------------- BEGIN CODE --------------

clearvars;
clearvars -global;


%% SETTINGS
% Specify the date, reduced-order model, reduction order and objects to be
% used to estimation here.

% Estimation window
yr      = 2020; % Year
mth     = 1;    % Month
dy      = 3;    % Day
hr      = 6;
mn      = 0;
sec     = 0;
nofDays = 25;   % Number of days

% Use high fidelity dynamical model
highFidelity = true;

% Reduced-order model
ROMmodel = 'JB2008_1999_2010';  % Name of reduced-order density model: JB2008_1999_2010, NRLMSISE_1997_2008 or TIEGCM_1997_2008
r  = 10;                        % Reduced order

% NORAD catalog IDs of objects used for estimation
selectedObjects = [22;614;932;1807;2153;2389;4221;4382;7337;8744;12138;12388;14483;20774;23278;41771;41772;41773;42989;43797]; % Radar: 20 objects
selectedObjects = sortrows(selectedObjects);

% Display date
datetime(yr,mth,dy)

%% SET PATHS
% *** SPECIFY YOUR SPICE TOOLBOX DIRECTORY HERE! ***
spicePath = fullfile('[SPICE TOOLKIT DIRECTORY]','mice'); 
global resultsDirPath measurementsPath
% *** SPECIFY OUTPUT DIRECTORY HERE! ***
resultsDirPath = fullfile('[RESULTS DIRECTORY]');
% *** SPECIFY FOLDER WITH LEOLABS MEASUREMENT DATA HERE! ***
measurementsPath = fullfile('[RADAR DATA DIRECTORY]');

addpath( 'AstroFunctions' );
addpath( 'Estimation' );
addpath( 'JB2008' );
addpath( 'ROMDensityModels' );
addpath( 'SpaceWeather' );
addpath( 'TLEdata' );
addpath( 'RadarData' );
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

runDensityEstimationRadar(yr,mth,dy,hr,mn,sec,nofDays,ROMmodel,r,selectedObjects,plotFigures,highFidelity);


%% Clear memory
% Clear cspice memory
cspice_kclear;

%------------- END OF CODE --------------
