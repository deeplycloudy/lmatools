'''
Charge Density Initiation Altitude: 
----------------------------------

As referenced in Bruning and MacGorman (2013), LMA data can be used to find
two-dimensionally channel lengths and areas of the LMA source points.

The flash area can be used to estimate the flash energy purely based on the geometry of 
individual flashes.

The energy dissipated by a lightning flash is defined using a capacitor model:

        W = [(rho**2 * d**3)/(2 * e_0)] * A

        where: rho == charge density
               d == distance between charge region centers (capacitor plates)
               e_0 == permitivity of vacuum
               A == convex hull area

The variable that cannot be simply obtained from LMA data is the charge density (rho). The 
following routines are used to estimate rho by method introduced in Boccippio [2002]. From the
initiation altitude of a flash, the critical electric field (defined by the Runaway Breakeven
Threshold) and critical surface and space charge densities can be retrieved.

The energy values obtained using this model have been shown to follow closely to estimates 
of inidivudal flash energies using a cloud-resolving model COMMAS, allowing for application of
the capacitor model of Bruning and MacGorman [2013] while minimizing the use of arbitrary physical
quantities and units. All energy estimates are in Joules [J].
'''
from __future__ import absolute_import
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class rho_retrieve(object):
    '''
    Critical Electric Field Method [Boccippio 2002]:
    -------------------
    '''
    
    def __init__(self, area, d, zinit, separation, constant=False, arbitrary_rho=0.4e-9):
        self.area      = area  #convex hull area
        self.d         = d     #distance between charge region centers
        self.l         = np.sqrt(area*1e6)  #flash length
        self.zinit     = zinit #altitude of initiation
        self.separation= separation
        self.e_0       = 8.85418782e-12
        
        #Initialized rho and energy values
        self.rho      = 0
        self.w        = 0
        self.constant = constant
        self.arbitrary_rho = arbitrary_rho

    def rho_br(self, z, separation):
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
    
    def energy_estimation(self, rho, d, area):
        w = (rho**2. * d**3. * area)/(2. * self.e_0) 
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









