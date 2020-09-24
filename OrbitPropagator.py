# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 01:03:03 2020

@author: marti
"""

import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import ode

import planetData as pd

class OrbitPropagator:
    def __init__(self,stateIn,tspan,dt,cb=pd.earth):
        self.r0 = stateIn[:3]
        self.v0 = stateIn[3:]
        self.stateIn = stateIn
        self.tspan = tspan
        self.dt = dt
        self.cb = cb
        
        

    def propOrbit(self):
        
        # Number of steps
        self.nSteps = int(np.ceil(self.tspan/self.dt))
    
        # Setup size of initial arrays
        self.state = np.zeros((self.nSteps,6))
        self.ts = np.zeros((self.nSteps,1))
        
        #initial conditions
        self.state[0] = self.stateIn
        self.ts[0] = 0
        self.step = 1
        
        # initiate solver
        self.solver = ode(self.orbitPropEQ)
        self.solver.set_integrator('lsoda')
        self.solver.set_initial_value(self.stateIn,0)
       # self.solver.set_f_params(self.cb['mu'])
        
        #propagate orbit
        while self.solver.successful() and self.step<self.nSteps:
            self.solver.integrate(self.solver.t+self.dt)
            self.ts[self.step] = self.solver.t
            self.state[self.step] = self.solver.y
            self.step += 1
    
    
        self.stateout = self.state[:,:]
        self.rf = self.state[:3]
        self.vf = self.state[3:]
        
        #plot(self.rf)
    
    
    def orbitPropEQ(self,t,state):
        # Description: Fundamental orbit equations of motion to propagated in an ODE propagator
        
        # Inputs:
        # t = time                                      [s]
        # stateIn = input state vector [COLUMN 6x1]     [km km/s] 
        # mu = gravitational parameter of central body  [km^3/s^2]
        
        rx,ry,rz,vx,vy,vz = state
        r = np.array([rx,ry,rz])
        
        norm_r = np.linalg.norm(r)
        
        ax,ay,az = -r*self.cb['mu']/norm_r**3
        
        stateOut = [vx,vy,vz,ax,ay,az] 
        return stateOut
    
    
    def plotOrbit(self, show_plot=False, save_plot=False, title = 'Test Title'):
        fig = plt.figure(figsize=(18,6)) 
        ax = fig.add_subplot(111,projection='3d')
        
        #plot trajectory
        ax.plot(self.rs[:,0],self.rs[:,1],self.rs[:,2],'w',label='Trajectory')
        ax.plot([self.rs[0,0]],[self.rs[0,1]],[self.rs[0,2]],'wo',label='Initial Position')
        
        # plot central body
        _u,_v = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
        _x = self.cb['radius']*np.cos(_u)*np.sin(_v)
        _y = self.cb['radius']*np.sin(_u)*np.sin(_v)
        _z = self.cb['radius']*np.cos(_v)
        ax.plot_surface(_x,_y,_z,cmap='Blues')
        
        # plot the x,y,z vectors
        l = self.cb['radius']*2
        x,y,z = [[0,0,0],[0,0,0],[0,0,0]]
        u,v,w = [[l,0,0],[0,l,0],[0,0,l]]
        ax.quiver(x,y,z,u,v,w,color='k')
        
        max_val = np.max(np.abs(self.rs))
        
        ax.set_xlim([-max_val,max_val])
        ax.set_ylim([-max_val,max_val])
        ax.set_zlim([-max_val,max_val])
        
        ax.set_xlabel(['X (km)'])
        ax.set_ylabel(['Y (km)'])
        ax.set_zlabel(['Z (km)'])
        
        ax.set_title(title)
        plt.legend()
        
        if show_plot:
            plt.show()
        
        if save_plot:
            plt.savefig(title+'.png',dpi = 300)