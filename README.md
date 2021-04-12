# DESTO
The Density Estimation Toolbox (DESTO) is a complete toolbox for Matlab that enables you to estimate the global thermospheric density using GPS, radar or two-line-element (TLE) data. In addition, three different reduced-order density models can be employed for the estimation. 

Copyright © 2021 by David Gondelach and Richard Linares


### License
This code is licensed under the GNU General Public License version 3 - see the [LICENSE](LICENSE) file for details.


### Acknowledgments
Initial work on the code was performed by Dr. Piyush M. Mehta, see https://doi.org/10.2514/1.G004793.

The MATLAB code for the Jacchia-Bowman 2008 model and solar radiation pressure and third-body perturbations was developed by Meysam Mahooti (copyright 2018) and was downloaded from https://www.mathworks.com/matlabcentral/fileexchange/56163-jacchia-bowman-atmospheric-density-model (version 2.0.0.0) and https://www.mathworks.com/matlabcentral/fileexchange/55167-high-precision-orbit-propagator (version 2.1.1.1).

The MATLAB code for the SGP4 model and several time and reference frame routines was developed by David Vallado (and others) and was downloaded from https://celestrak.com/software/vallado-sw.php.

The Earth Gravitational Model 2008 (EGM2008) coefficients were obtained from the NGA's Office of Geomatics: https://earth-info.nga.mil.

The toolbox makes use of NASA's SPICE Toolkit (version N0066) for Matlab, see https://naif.jpl.nasa.gov/naif/toolkit.html.


### References
The code in this repository corresponds to two publications in the AGU journal Space Weather. The thermospheric density modeling and estimation techniques using two-line element data are described in:
```
@article{gondelach2020tle,
  author = {Gondelach, David J. and Linares, Richard},
  title = {Real-Time Thermospheric Density Estimation Via Two-Line-Element Data Assimilation},
  journal = {Space Weather},
  doi = {10.1029/2019SW002356},
  url = {https://doi.org/10.1029/2019SW002356}
}
```
see https://doi.org/10.1029/2019SW002356.

The density estimation using radar range and range-rate measurements and GPS position measurements is described in:
```
@article{gondelach2021radargps,
  author = {Gondelach, David J. and Linares, Richard},
  title = {Real‐Time Thermospheric Density Estimation Via Radar And GPS Tracking Data Assimilation},
  journal = {Space Weather},
  doi = {10.1029/2020SW002620},
  url = {https://doi.org/10.1029/2020SW002620}
}
```
see https://doi.org/10.1029/2020SW002620.

### Installation instructions
1. Download the DensityEstimation Matlab code
2. Download and install SPICE Toolkit for Matlab: https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html
3. Set the path to the SPICE Toolkit directory in mainDensityEstimation[GPS/Radar/TLE].m
4. Download SPICE kernels (i.e. ephemeris files) from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/ and put them in the folder Data. See links below.
5. Download space weather file from Celestrak and put in folder Data: https://www.celestrak.com/SpaceData/SW-All.txt
6. Download Earth orientation data file from Celestrak and put in folder Data: https://www.celestrak.com/SpaceData/EOP-All.txt
7. Download 2 space weather files needed for the JB2008 model and put in folder Data: http://sol.spacenvironment.net/jb2008/indices/SOLFSMY.TXT  and  http://sol.spacenvironment.net/jb2008/indices/DTCFILE.TXT 
8. If you want to use automatic download of TLE data, then specify your “www.space-track.org” username and password in runDensityEstimationTLE.m, line 95-96. (Alternatively, you can download the TLE data manually and put them in the folder TLEdata with naming convention: [NORADID].tle, e.g. 12388.tle )
9. For each object used for estimation, specify the ballistic coefficient (BC) in the text file: Data/BCdata.txt


### Run instructions
The code is written to use GPS data in json format, radar measurements provided by Leolabs Inc. (https://platform.leolabs.space), and TLE data as provided by www.space-track.org.
To run the code, follow these instructions:
1. Open mainDensityEstimation.m (TLE), mainDensityEstimationGPS.m (GPS), or mainDensityEstimationRadar.m (Radar).
2. Specify the path to the folder with GPS or radar measurement data
3. Specify the path to the output folder where results will be saved
4. Specify the time window for density estimation by setting the start date and number of days
5. (Optionally) select the reduced-order density model and the reduction order (default order: r=10)
6. (Optionally) select the objects to use for density estimation
7. Run mainDensityEstimation[GPS/Radar/TLE]


### Ephemeris files
Download the following ephemeris files and put them in the Data folder:
* latest_leapseconds.tls:  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/
* de430.bsp:  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
* earthstns_itrf93_201023.bsp (previously: earthstns_itrf93_050714.bsp):  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/stations/
* pck00010.tpc, earth_fixed.tf, earth_200101_990628_predict.bpc (previously: earth_070425_370426_predict.bpc), earth_720101_070426.bpc, earth_latest_high_prec.bpc:  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/


### Technical notes
The speed of the code has not been optimized. A 10-day density estimation may required several hours of computation.

MATLAB R2018b (Version 9.5) was used to develop the code.



David Gondelach, April 2021
email: davidgondelach@gmail.com
