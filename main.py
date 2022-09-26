import numpy as np
import spiceypy as spice
from scipy import constants
import math
 


'''
Profiler() Analyses the shape of the 3D profile by producing a the cartesian coordinates of the
tangent point for every second of the occultation test.
Then it finds geometric parameters of this profile
Outputs:
    profile - (unused) 3d coordinates of the tangent point for every second of the occultation 
    lowestHeightTime -(unused) seconds in the test for when the tagent point is at its lowest point. (occs can start below the horizon during signal unlock)
    peakHeightTime - seconds in the test for when the tagent point reaches roughly the peak of the m2 ionospheric peak
    plateauTime - seconds into the test when the tangent point starts to plateau out and rises slower that 50 m/s
    grazingAngle - geometric angle for the lowest point and the peak ionospheric point. close to 90 degrees means it is steap and the residuum will be amplified 
    maxVelDop - The maximum geometric doppler expected during the occultation. This is an important perameter as a large doppler measurement will effect the windows size of
        the fft used later in the processing. 
    startDistance - MEX-TGO distance in km
    endDistance - MEX-TGO distance in km
    maxAlt - The maximum tangent point altitude above the martian spheriod (dsk model is not used). The higher the better. This is capped at the hight of the lowest spacecraft.
        (usually TGO ~400 km)
'''
def Profiler(et,duration,scheme, sv):

    marsrad = spice.bodvrd(sv.front, 'RADII', 3)# Extract the triaxial dimensions of Mars
    altlist = []
    veldopp = []

    modifier = -1 if scheme =='ingress' else 1
    et = et+duration if scheme =='ingress' else et
    elapsed = modifier * duration# because we always iterate spacially upwards, this decides if we are counting back or forwards
   
    for i in range(0,elapsed, modifier):
        # Find relative positions of TGO and MEX
        [targetpos, _] = spice.spkpos(sv.front, et + i, sv.fframe, 'LT+S', sv.target)
        [states, _] = spice.spkezr(sv.target, et + i, sv.fframe, 'NONE', sv.obs)
        [obspos, _] = spice.spkpos(sv.front, et + i, sv.fframe, 'LT+S', sv.obs)
        
        #calculate maximum expected geometric doppler
        velocity = states[3:6]
        relativespeed = np.linalg.norm(velocity)
        vd = (relativespeed/constants.c) * 437.1e9# e9 because we are converting from km to m (SPICE outputs km, but constants in m)
        veldopp = np.append(veldopp,vd)
        maxVelDop = np.max(veldopp)
        
        # Find the unit vector between the SCs
        sc2scvector = states[0:3]
        displacement = np.linalg.norm(sc2scvector)
        unitvector = np.true_divide(sc2scvector, displacement)
        
        #find the start and end distances
        if i==0:
            startDistance = displacement
        elif i==elapsed-modifier:
            endDistance = displacement

        # Find the point this sc vector is closest to the Mars
        [profilesurfacepoint, alt] = spice.npedln(
            marsrad[1][0], marsrad[1][1], marsrad[1][2], targetpos, unitvector)

        altlist = np.append(altlist, alt)

        # Find distance between the tangent point and MEX (obs)
        tangent2mexdist = np.linalg.norm(obspos - profilesurfacepoint)

        # Need to add the altitude to the surface bound point. 'profilesurfacepoint' is a vector
        # from Mars' barrycenter to the surface, therefore
        # The vector is also the direction of altitude, so product must be added to the surface point
        tangentpointunitvector = profilesurfacepoint / \
            np.linalg.norm(profilesurfacepoint)
        tangentpoint = (tangentpointunitvector * alt) + profilesurfacepoint

        # the SPICE function 'npedln' does not calculate along a vector between the two SCs, the vector continues through TGO. This is correced
        # by the following. If the tangent2mex distance > than the tgo2mex distance, then the vector has passed though the TGO. End function when this occurs
        if tangent2mexdist > (displacement-50):
            tangentpoint = targetpos
            highesttime = i
            # break

        # The peak height can also be defined by when the profile begins to planteau [when the profile raises by less than
        # 50 m/s ]
        if altlist[i-1]+0.5 > altlist[i]:
            plateauTime = i

        # Form array of profile coordinates
        if i == 0:
            x = int(tangentpoint[0])
            y = int(tangentpoint[1])
            z = int(tangentpoint[2])
        else:
            x = np.append(x, int(tangentpoint[0]))
            y = np.append(y, int(tangentpoint[1]))
            z = np.append(z, int(tangentpoint[2]))

    # Form a array of cartesian coordinates, transpose for convienience
    profile = [x.T, y.T, z.T]#profile not required for Occ Geometry Tool, but useful for future

    peakHeightTime = np.where(altlist>120)#find the peak ionospheric intercept (roughly 120km)
    peakHeightTime = peakHeightTime[0][0]

    lowestHeightTime = np.where(altlist>1)#start and end times may be below the surface (electra unlock)
    lowestHeightTime = np.min(lowestHeightTime)

    #calculate the grazing angle    
    lowestpoint = np.asarray([x[0], y[0], z[0]])
    highestpoint = np.asarray([x[peakHeightTime], y[peakHeightTime], z[peakHeightTime]])
    # find the vector from start to end of profile
    maxsc2scvector = (highestpoint - lowestpoint)
    norm = np.linalg.norm(maxsc2scvector)
    endunitvector = maxsc2scvector/norm
    norm = np.linalg.norm(lowestpoint)
    # provide unit vector for the lowestpoint, prepresenting the direction from Mars' center to the begining point
    startunitvector = lowestpoint/norm
    # producing the anlge between Mars' center ->start & start-> end
    angle = spice.vsep(startunitvector, endunitvector)
    grazingAngle = angle * (180/math.pi)

    maxAlt = np.max(altlist)

    return profile, lowestHeightTime,peakHeightTime, plateauTime, grazingAngle, maxVelDop, startDistance,endDistance, maxAlt


'''
SchemeChecker() checks if this occultation is an egress or an egress. It simply checks if the tangent point is higher or lower at the start 
when compared to the end. A typical occultation test will be shorter than a TGO orbit, so there is no risk of the tangent point coming back down after an egress.
Because the other functions work by counting up in altitude, this scheme output is important because it states whether the time should be going up or down 
outputs:
    scheme - a string which classifies whether this specific occultation is an egress or ingress
'''
def SchemeChecker(et,duration,sv):
    alt = [0,0]
    i=0
    for time in [et,et+duration]:
        marsrad = spice.bodvrd(sv.front, 'RADII', 3)
        [targetpos, _] = spice.spkpos(sv.front, time, sv.fframe, 'None', sv.target)
        [sc2scvector, _] = spice.spkpos(sv.target, time, sv.fframe, 'NONE', sv.obs)

        # Find the unit vector between the SCs
        displacement = np.linalg.norm(sc2scvector)
        unitvector = np.true_divide(sc2scvector, displacement)

        # Find the point this unit vector is closest to the Mars
        [_, alt[i]] = spice.npedln(
            marsrad[1][0], marsrad[1][1], marsrad[1][2], targetpos, unitvector)
        i+=1
    scheme ='ingress' if alt[0]>alt[1] else 'egress'

    return scheme


'''

'''
def PointingAngles(start,duration,sv):
    #start
    relativePositionofTGO = spice.spkpos(sv.target, start, 'TGO_SPACECRAFT', 'LT', sv.obs)
    startDistance = np.linalg.norm(relativePositionofTGO[0])
    unitvector = relativePositionofTGO[0]/startDistance
    angleRad  = spice.vsep([0, -1, 0],unitvector) # TGO pointing ref is -y
    TGOStartAngle = 180 -np.rad2deg(angleRad)

    relativePositionofMEX = spice.spkpos(sv.obs, start, 'MEX_SPACECRAFT', 'LT', sv.target)
    unitvector = relativePositionofMEX[0]/ np.linalg.norm(relativePositionofMEX[0])
    angleRad  = spice.vsep([0, 0, 1],unitvector)#MEX pointing ref is +Z
    MEXStartAngle = 180 -np.rad2deg(angleRad)   

    #end
    end = start+duration
    relativePositionofTGO = spice.spkpos(sv.target, end, 'TGO_SPACECRAFT', 'LT', sv.obs)
    endDistance = np.linalg.norm(relativePositionofTGO[0])
    unitvector = relativePositionofTGO[0]/endDistance 
    angleRad  = spice.vsep([0, -1, 0],unitvector)# TGO pointing ref is -y
    TGOEndAngle = 180 -np.rad2deg(angleRad)

    relativePositionofMEX = spice.spkpos(sv.obs, end, 'MEX_SPACECRAFT', 'LT', sv.target)
    unitvector = relativePositionofMEX[0]/ np.linalg.norm(relativePositionofMEX[0])
    angleRad  = spice.vsep([0, 0, 1],unitvector)#MEX pointing ref is +Z
    MEXEndAngle = 180 -np.rad2deg(angleRad)

    #print(f'TGO angle = {TGOStartAngle} -> {TGOEndAngle} \n MEX angle = {MEXStartAngle} -> {MEXEndAngle}')
    return TGOStartAngle,TGOEndAngle,MEXStartAngle,MEXEndAngle, startDistance, endDistance


'''
GeoSpec() is a simplified version profiler, which provies some geometric parameters but only for specific input times (does not iterate through the whole test)
Output:
    lon - longitude of tangent point in degrees
    lat - latitude of tangent point in degrees
    displacement - distance between the two spacecraft
    nearestpoint - 3D coordinate of tangent point in the martian intertial frame
    veldopp - Geometric doppler in Hz assuming the Tx is 437.1 MHz 
'''
def GeoSpec(et, sv):
    
    [tgopos, _] = spice.spkpos(sv.front, et, sv.fframe, 'NONE', sv.target)
    [sc2scvector, _] = spice.spkezr(sv.target, et, sv.fframe, 'NONE', sv.obs) # states relative to eachother
    
    velocity = sc2scvector[3:6]
    relativespeed = np.linalg.norm(velocity)
    veldopp = (relativespeed/constants.c) * 437.1e9# e9 because we are converting from km to m (SPICE outputs km, but constants in m)
    
    sc2scvector = sc2scvector[0:3]
    displacement = np.linalg.norm(sc2scvector)
    sc2scunitvector = np.true_divide(sc2scvector, displacement)
    
    marsrad = spice.bodvrd(sv.front, 'RADII', 3)# Extract the triaxial dimensions of Mars
    # For the ray that connects MEX and TGO, find the point on this ray that is closest to the Martian surface
    [nearestpoint, alt] = spice.npedln(marsrad[1][0], marsrad[1][1], marsrad[1][2], tgopos, sc2scunitvector)
    [radius, lon, lat] = spice.reclat(nearestpoint)
    lon = 180 - (lon * (-180 / math.pi))# Rad -> Deg , frame inversion required (hence the negative 180)
    lat = lat * (-180 / math.pi)

    return lon, lat, displacement, nearestpoint, veldopp

''' 
SolarZenithAngle() finds the **SIGNED** SZA for a given 3D point

Normally SZA is an unsigned integer, indicating the angle to the sun is. However this doesnt give
any indication on whether it is AM and PM). Therefore the SZA is calculated twice to see if the 
angle is getting bigger or smaller at a later time. If smaller, then it is the AM (as
the tangent point is moving towards the noon) This is important for martian occultation as the
ionosphere can vary wildly from sunrise-> sunset.
Outputs:
    sza = angle to the sun from a given 3D coordinate. Positive value means pm. Negative means am
'''

def SolarZenithAngles(et, nearestpoint, sv):

    subsolarpoint, _, _ = spice.subslr('INTERCEPT/ELLIPSOID', sv.front, et, sv.fframe, 'NONE', sv.target)  # Where is the sun?
    sza = spice.vsep(-subsolarpoint, nearestpoint)  # angle between sun and sc
    latersubsolarpoint, _, _ = spice.subslr('INTERCEPT/ELLIPSOID', sv.front, et + 30, sv.fframe, 'NONE', sv.target)  # where is the sun later?
    latersza = spice.vsep(-latersubsolarpoint, nearestpoint)

    if sza < latersza:  # if the sun is moving away from the tangent point (PM)
        sza = sza * (180/math.pi)
        
    else:
        sza = sza * (-180 / math.pi)  # neg SZA means mornings

    return sza
