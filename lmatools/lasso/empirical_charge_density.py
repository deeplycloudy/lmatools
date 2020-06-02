'''
Charge Density Power Law Function: 
(modified from: http://www.bdnyc.org/2012/04/fitting-a-power-law-to-data/)
----------------------------------

As referenced in Bruning and MacGorman (2013), LMA data can be used to find
two-dimensionally channel lengths from the convex hull area of the LMA source points.

This is defined as:

        L = sqrt(A_hull)

This length can be used to estimate the energy dissipated by a lightning flash using a 
capacitor model:

        W = [(rho**2 * d**3)/(2 * e_0)] * L**2

        where: rho == charge density
               d == distance between charge region centers (capacitor plates)
               e_0 == permitivity of vacuum
               L == the length derived from the convex hull area

From this estimation, however, the variable (rho) was not computed or obtained. The 
following method was emperically obtained from work done to compare the influence of 
geometry on the electrostatics for variable flash sizes (Vicente, 2016 - MS Thesis).
Charge densities obtained from least squares optimization under the constraints of Charge
Conservation and an Electric Field bounded by the Runaway Breakeven Threshold showed that (rho)
was correlated to the size of the flash (defined by it's convex hull volume).

From these obtained charge densities an expression for (rho) was emperically defined as a power law 
function of (L), expressed as:


            rho(L) = a*L^(b) + c

#####EDIT: THESE COEFFICIENTS USED KM not M, SO, CHANGED:
Where from least squares optimization, reducing the error between (rho_data - rho_guess) gave the 
coefficients:

            #a = 1.79643909279e-08
            #b = -0.00411767345222
            #c = -1.74885409779e-08 (almost neglible)

             a = 2.83385909218e-09
             b = -0.0227528656802
             c = -1.89660631464e-09

Therefore:

            #rho(L) = (1.796e-8)L^(-0.004)+(-1.7489e-8) (KM) ---> C/km^3
            rho(L) = (2.834e-09)L^(-0.023)+(-1.897e-09) (M)  ---> C/m^3

This function allows the retrieval of charge densities for an arbitrary, or known set of lighting 
channel lengths obtained from the convex hull area method, resembling a capacitor plate). The electrostatic
energy dissipated between these plates can the be found as in Bruning and MacGorman (2013).
'''
from __future__ import absolute_import
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class rho_retrieve(object):
    '''
    Optimization Method:
    -------------------
    from scipy.optimize import leastsq
   
    powerlaw = lambda x, amp, coeff, index: amp * (x**index) +coeff
    fitfunc = lambda p, x: p[0]*(x**p[1]) + p[2]
    errfunc = lambda p, x, y: (y - fitfunc(p, x))

    optimized coefficients for fit:
    index = pfinal[1]
    amp = pfinal[0]
    coeff = pfinal[2]

    '''
    
    def __init__(self, area, d, zinit,separation,constant, arbitrary_rho):
        self.area = area
        self.d = d  #distance between charge region centers
        self.l = np.sqrt(area*1e6)  #in m
        self.zinit = zinit
        self.separation=separation
        # self.rho = []
        # self.w = []#electrostatic energy for capcitor plates
        self.e_0 = 8.85418782e-12
        
        self.rho = 0
        self.w = 0
        self.constant = False
        self.arbitrary_rho = 0.4e-9

    def rho_br(self,z,separation):
        '''
        Compute critical electric field from initiation altitudes.
        From the electric field, compute critical charge density.
        Note: does not assume RBE as breakdown mechanism, simply
          allows for approximation of electric field upon
          flash initiation. 
          
          args: z = initiation altitude in km 
                     MUST CONVERT FROM INDEX TO ALTITUDE 
                     USING GRID SPACING:
                     
                     z(km) = (z*125) * 1e-3
                     
          returns: rho_critical in kV/m
        '''
        z = z*1e-3
        rho_a = 1.208 * np.exp(-(z/8.4))
        efield=(167.*rho_a)*1e3
        sig_crit=2*efield*self.e_0
        rho_crit=sig_crit/separation
        return(rho_crit)
    
    #The power law idea is no longer used: Use methods in Salinas et al [In Progress]    
    def rho_powerlaw(self, l):
        #Convergence after 801 iterations of the solution: From leastsq output
        a = 2.83385909218e-09
        b = -0.0227528656802
        c = -1.89660631464e-09
        p_law = a * l **(b) + c
        self.rho = p_law
        return p_law
        
    def energy_estimation(self, rho, d, area):
        w = (rho**2. * d**3. * area**2.)/(2. * self.e_0) 
        self.w = w
        #print(w,area,area*1e-6,d,rho)
        return w
        
    def calculate(self):
        self.rho = self.rho_br(self.zinit,self.separation)
        if self.constant == False:
            self.w = self.energy_estimation(self.rho, self.separation, self.area)
        else: 
            self.w = self.energy_estimation(self.arbitrary_rho, self.separation, self.area)
        return self.rho, self.w
        
    # def plot_rho(self):
    #     rho = self.rho_powerlaw(self.l*1e-3)
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111)
    #     sort = np.argsort(self.l)
    #     ax.plot(self.l[sort], self.rho_powerlaw(self.l[sort]), 'k-', alpha=0.5)
    #     ax.plot(self.l, rho, 'ok', alpha=0.5)
    #     plt.show()
        
    # def plot_w(self):
   #      rho = self.rho_powerlaw(self.l*1e-3)
   #      w = self.energy_estimation(rho, self.d, self.l)
   #
   #      fig = plt.figure()
   #      ax = fig.add_subplot(111)
   #      sort = np.argsort(self.l)
   #      ax.semilogy(self.l, w, 'or', alpha=0.5)
   #      plt.grid(alpha=0.5)
   #      plt.show()









