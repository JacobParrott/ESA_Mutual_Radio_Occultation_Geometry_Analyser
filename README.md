
# MEX-TGO Mutual Radio Occultation Geometry Analyser
![Cartoon of Profile formed between MEX and TGO, made with a custom Cosmographia tool](https://github.com/JacobParrott/Mutual_Radio_Occultation_Geometry_Analyser/blob/main/images/CoverImage.png)
## A tool to analyse occultation geometry to aid in observation planning.

## Usage
Interaction with this tool happens through the 'inputFile.txt' file. It accepts multiple input times along with corresponding durations (in seconds).These times are in European format with a space followed by the duration. eg '19/07/2022 19:49:30 600'
The output is in the form of an excel table, with geometric parameters for each column.
The only lines that must be altered by the end user are in *setup.py lines 29:30* as you must point to your own local repository of SPICE kernels.


## Outputs
+ Day of Year
+ Date
+ UTC Start
+ UTC Occultation epoch
+ UTC End 
+ Spacecraft Distance at Start(km)
+ Spacecraft Distance at End (km)
+ Range of MEX +Z angles (degrees)
+ Range of TGO -Y angles (degrees)
+ Latitude of tangent point at occultation epoch (degrees North)
+ Longitude of tangent point at occultation epoch (degrees East)
+ Local Time (in 24 hour and minutes)
+ Solar Zenith Angle (degrees)
+ Expected maximum Doppler shift during the test (Hz)
+ The type of occultation (ingress/egress)
+ MEX Melacom warmup time (hardcoded to output 15 minutes)
+ Climbing rate of tangent point to M2 ionospheric peak (km/s)
+ Grazing Angle of occultation profile (degrees) 
![Example Excel Output file when the example inputFile.txt is used](/images/exampleoutput.png)

>Any of these outputs can be verified with WebGeoCalc (https://wgc.jpl.nasa.gov:8443/webgeocalc/#NewCalculation).



*This project has been written for mutual radio occultation of the Martian atmosphere. Where radio signals are  sent from Mars Express (MEx) to ExoMars' Trace Gas Orbiter (TGO). This can be treated as an example, where the satellites and even the planet are up to the end user's discretion. Fortunately, with the SPICE toolset, everything in the analysis can be changed with simple alterations of the starting variables.*





## Prerequisites

+ Some knowledge of [SPICE](https://naif.jpl.nasa.gov/naif/tutorials.html). This may take some time to understand but is very much worth it. SPICE is a remarkably powerful free-to-use tool.
+ Install the correct SPICE kernels, these are the datasets that the SPICE application reads. They include information on spacecraft, ground-station and planetary ephemerides, planetary shapes, onboard instrument details, leap-second tables, reference frame conversions and *much more*. Depending on the mission the user is interested in, the download location will vary. Space agencies will manage the SPICE kernels for their own mission. [ESA](https://www.cosmos.esa.int/web/spice) manage the kernels for [MEx](https://www.cosmos.esa.int/web/spice/spice-for-mex) and [TGO](ftp://spiftp.esac.esa.int/data/SPICE/ExoMars2016/).


### Dependencies

+ Spicypy
+ Pandas
+ Numpy
+ Math
+ scipy


## Feedback
Feel free to file an issue. Contributors and feature requests are always welcome.
