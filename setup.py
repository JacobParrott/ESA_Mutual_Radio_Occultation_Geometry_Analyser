import spiceypy as spice
import spiceypy.utils.support_types as stypes
import pandas as pd
import math

import main

def setup():

    class SpiceVariables:
        obs = '-41'  # NAIF code for MEX (-41)
        # NAIF code for TGO (-143)['EARTH'/'SUN'/ a groundstation etc]
        target = '-143'
        obsfrm = 'IAU_MARS'
        abcorr = 'NONE'
        crdsys = 'LATITUDINAL'
        coord = 'LATITUDE'
        stepsz = 1.0  # Check every 1 seconds if there is an occultation
        MAXILV = 10  # Max number of occultations that can be returned by gfoclt
        bshape = 'POINT'  # Rx shape
        fshape = 'ELLIPSOID'
        front = 'MARS'
        fframe = 'IAU_MARS'
        TFMT = 'YYYY-MM-DD HR:MN:SC'  # Format that Cosmographia understands
    sv = SpiceVariables()



    PathtoMetaKernel1 = 'C:/SPICEData/exomars2016/kernels/mk/em16_ops.tm'
    PathtoMetaKernel2 = 'C:/SPICEData/mars-express/kernels/mk/MEX_OPS.tm'
    spice.furnsh(PathtoMetaKernel1)
    spice.furnsh(PathtoMetaKernel2)


    with open('inputFile.txt') as f:
        lines = f.readlines()

    occs = pd.DataFrame(index=[lines], columns=['DoY','Date','StartUTC','OccUTC','EndUTC','StartDistance (km)',
                                                'EndDistance (km)','MEX +Z','TGO -Y','Latitude(°N)','Longitude (°E)',
                                                'LocalTime','SolarZenithAngle(°)','VelocityDoppler','MaxAltitude',
                                                'OccType','MELACOMWarmup (min)','ClimbingRate(km/s)','GrazingAngle'])

    for line in lines:    
        #convert european format to us SPICE format
        SPICEInput = line[4:7] + line[1:4] + line[7:18]
        duration = int(line[19:22])
        print('Now Processing: ',SPICEInput)

        #Input times do not state the exact moment occultation occures.  
        start = spice.str2et(SPICEInput)
        scheme = main.SchemeChecker(start,duration,sv) # check if ingress or egress
        end = start +600
        etbeg = start-300
        etend = start + 1300

        # Form SPICE interval windows that gfoclt can populate
        window = stypes.SPICEDOUBLE_CELL(2)
        spice.wninsd(etbeg, etend, window)
        occwindow = stypes.SPICEDOUBLE_CELL(sv.MAXILV)

        # find occultation windows between the dates listed above [ most comp cost in this function]
        spice.gfoclt('ANY', sv.front, sv.fshape, sv.fframe, sv.target,
                    sv.bshape, '', sv.abcorr, sv.obs, sv.stepsz, window, occwindow)
        winsiz = spice.wncard(occwindow)  # Find cardinality (number of windows)
        if winsiz==0:
            print('No Occultation Present')
            continue

        for i in range(winsiz):
            [ingress, egress] = spice.wnfetd(occwindow, i) # extract the begining and ends of the windows
        OccultationEpoch = ingress if scheme =='ingress' else egress 
        UTCActualOcc = spice.et2utc(OccultationEpoch,'C',3)
        print(f'\t Occulatation at {UTCActualOcc}')
        profile, lowesttime, highesttime, platuauTime, GrazingAngle, maxVelDop,startDistance,endDistance, maxAlt  = main.Profiler(start,duration,scheme, sv)

        AngleTGOStart ,AngleTGOEnd,AngleMEXStart,AngleMEXEnd,startDistance, endDistance= main.PointingAngles(start,duration,sv)

        usefulPercentage = (highesttime/duration)*100#perhaps useful in future

        date_time,localtime,lat,lon,veldopp,dist,sza =(['']*3 for i in range(7))#initialise output lists
        i=0
        for time in [start, OccultationEpoch,end]: 

            # try to plot the location on a map with cartopy
            lon[i], lat[i], dist[i], nearestpoint, veldopp[i] = main.GeoSpec(time, sv)  # FUNCTION NEEDS NEW NAME
            sza[i] = main.SolarZenithAngles(time, nearestpoint, sv)
            hour, minute, _, _, _ = spice.et2lst(time, 499,  lon[i] * (math.pi/180), 'PLANETOCENTRIC')  # 499 = spice code for Mars
            hour = '0' + str(hour) if hour<10 else str(hour)
            minute = '0' + str(minute) if minute<10 else str(minute)
            localtime[i] = "%s:%s" % (hour, minute)
            a = (spice.et2utc(time, 'ISOD',1))
            date_time[i] = a[9:17]
            DoY = a[5:8]
            i+=1

        #output to a pandas dataframe
        occs.loc[line,'DoY'] = DoY
        occs.loc[line,'Date'] = line[1:11]
        occs.loc[line,'StartUTC'] = date_time[0]
        occs.loc[line,'OccUTC'] = date_time[1]
        occs.loc[line,'EndUTC'] = date_time[2]
        occs.loc[line,'StartDistance (km)'] = round(startDistance)
        occs.loc[line,'EndDistance (km)'] = round(endDistance)
        occs.loc[line,'MEX +Z'] = f'{AngleMEXStart:.2f} - {AngleMEXEnd:.2f}' 
        occs.loc[line,'TGO -Y'] = f'{AngleTGOStart:.2f} - {AngleTGOEnd:.2f}' 
        occs.loc[line,'Latitude(°N)'] = lat[1]
        occs.loc[line,'Longitude (°E)'] = lon[1]
        occs.loc[line,'LocalTime'] = localtime[1]
        occs.loc[line,'SolarZenithAngle(°)'] = sza[1]
        occs.loc[line,'VelocityDoppler'] =round(maxVelDop)
        occs.loc[line,'MaxAltitude'] = round(maxAlt)
        occs.loc[line,'OccType'] = scheme
        occs.loc[line,'MELACOMWarmup (min)'] = 15
        occs.loc[line,'ClimbingRate(km/s)'] =120/highesttime
        occs.loc[line,'GrazingAngle'] = round(GrazingAngle)



    # FYI -will not save to exel file if it is open in another window
    occs.to_excel("OccultationGeometryOutput1.xlsx")



if __name__ =='__main__':
    setup()
