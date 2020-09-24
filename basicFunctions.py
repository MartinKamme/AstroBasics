# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 07:46:34 2020

@author: marti
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as NP
from scipy.integrate import ode

def orbitPropEQ(t,stateIn,mu):
    # Description: Fundamental orbit equations of motion to propagated in an ODE propagator
    
    # Inputs:
    # t = time                                      [s]
    # stateIn = input state vector [COLUMN 6x1]     [km km/s] 
    # mu = gravitational parameter of central body  [km^3/s^2]
    
    rx,ry,rz,vx,vy,vz = stateIn
    r = np.array([rx,ry,rz])
    
    norm_r = np.linalg.norm(r)
    
    ax,ay,az = -r*mu/norm_r**3
    
    stateOut = [vx,vy,vz,ax,ay,az] 
    return stateOut


def coes2RV(h,inc,RAAN,ecc,per,theta,mu):
    # Description: Converts classical orbital elements into position(R) and velocity vectors(V) in ECI
    # Sources: Curtis Orbital Mechanics for Engineering Students
    
    # Inputs:
    # h = Specific angular momentum             [km^2/s]
    # inc = inclination                         [rad]
    # RAAN = Right Acsension of ascending node  [rad]
    # ecc = Eccentricity 
    # per = argument of perigee                 [rad]
    # theta = true anomaly                      [rad]
    # mu = Planet's gravitational parameter     [km^3/s^2]
    
    # Outputs:
    # R: Postion vector in ECI                  [km]
    # V: Velocity vector in ECI                 [km/s]
    
    # State Vectors in Perifocal Coordinates
    rx = ((h**2)/mu)*(1/(1+ecc*np.cos(theta)))*np.array( ((np.cos(theta)),(np.sin(theta)),(0)) )
    vx = (mu/h)*np.array( ((-np.sin(theta)),(ecc+np.cos(RAAN)),(0)) )
    
    
    # Direction Cosine Matrix
    DCM = np.array( ((np.cos(per),np.sin(per),0), (-np.sin(per),np.cos(per),0), (0,0,1)) ) \
    @ np.array( ((1,0,0), (0,np.cos(inc),np.sin(inc)), (0,-np.sin(inc),np.cos(inc))) ) \
    @ np.array( ((np.cos(RAAN),np.sin(RAAN),0), (-np.sin(RAAN),np.cos(RAAN),0), (0,0,1)) )
    
    # Transformation Matrix
    Dcm = NP.pinv(DCM)
    
    # ECI R
    R = Dcm @ (rx) #np.transpose
    
    # ECI V
    V = Dcm @ (vx)
    
    return [R,V]


def RV2coes(R,V,mu):
    # Description: Converts R ECI and V ECI into classical orbital elements
    # Sources: Curtis Orbital Mechanics for Engineering Students
    
    # Inputs:
    # R = position vector in ECI    [COLUMN]    [km]
    # V = velocity vector in ECI    [COLUMN]    [km]
    # mu = planet's gravitational parameter     [km^3/s^2]
    
    # Outputs:
    # a = semi-major axis                       [km]
    # ecc = eccentricity
    # inc = inclination                         [rad]
    # RAAN = Right Ascension of Ascending Node  [rad]
    # per = argument of perigee                 [rad]
    # TA = Mean Anomaly                         [rad]
    # argLat = Argument of Latitude             [rad]
    # lonPer = Longitude of Periapsis           [rad]
    # p = semilatus Rectum                      [km]
 #   R = np.array(list(R.flat))
 #   V = np.array(list(V.flat))
    
    eps = np.exp(-10)
    
    r = NP.norm(R)
    v = NP.norm(V)
    
    vr = (np.matmul(R.T,V))/r
    
    H = np.cross(R,V)
    h = NP.norm(H)
    
    inc = np.arccos(H[2]/h)
    
    N = np.cross(np.array( [0,0,1] ),H)
    n = NP.norm(N)
    
    if n != 0:
        RAAN = np.arccos(N[0]/n)
        if N[1] < 0:
            RAAN = 2*np.pi - RAAN
    else:
        RAAN = 0
    
    E = (1/mu)*((v**2 - mu/r)*R - r*vr*V)   # eccentricity vector
    ecc = NP.norm(E)
    
    if n != 0:
        if ecc > eps:
            per = np.arccos(np.dot(N,E)/n/ecc)
            if E[2] < 0:
                per = 2*np.pi - per
        else:
            per = 0
    else: 
        per = 0
          
    if ecc > eps:
        TA = np.arccos(np.dot(E,R)/ecc/r)
        if vr < 0 :
            TA = 2*np.pi - TA
    else:
        cp = NP.cross(N,R)
        if cp[2] >= 0:
            TA = np.arccos(np.dot(N,R)/n/r)
        else: 
            TA = 2*np.pi - np.arccos(np.dot(N,R)/n/r)
    
    a = (h**2)/mu/(1-ecc*2)
    
    coes = [h, ecc, RAAN, inc, per, TA, a]
    
    return coes