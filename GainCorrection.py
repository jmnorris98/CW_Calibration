# -*- coding: utf-8 -*-
"""
Created on Tue May  9 11:03:19 2023

Minute by Minute Gain Calibration

This code takes the outputs of the power meters as well as the spectrometer 
.npy files as generated by Phillip Blacks python package, LBASS.sh and
calculates the gains and gain corrected data.

See 'Implementation of a Continuous-Wave Calibration Scheme for the 
L-Band All Sky Survey' MPhys Report for the mathematical details of how this
works.

Requires files for the power meter data as well as the data produced by LBASS.sh

Import CW_ReadIn_Funcs.py to read data from the power meter and LBASS.sh files.

@author: Jordan Norris
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import dates as mpl_dates
from astropy.time import Time
from astropy.time import TimeDelta
from scipy.fft import fft, fftfreq, rfft, rfftfreq
from scipy.optimize import curve_fit
import csv
from datetime import datetime, timedelta, date, time
from matplotlib.ticker import AutoMinorLocator
import os

import CW_ReadIn_Funcs

#------------------------------------------------------------------------------
#Directories
#Example Directory: Change if needed


PM_FILE_PATH = r'D:\ExampleDirectory\PowerMeterFiles/'
RPG_FILE_PATH = r'D:\ExampleDirectory\RPGFiles/'
IMAGE_SAVE_PATH = r'D:\ExampleDirectory\Images'


#PARAMETERS -----------

"""
Enter Paramters such as the desired averaging time, the start date time of the 
observation, its duration as well as the power meter ID.
"""

AVERAGING_TIME = 60 # Seconds
START_DATE_TIME = datetime.combine(date(2023, 4, 3), time(18, 00)) #UTC TIMES
END_DATE_TIME = datetime.combine(date(2023, 4, 4), time(6, 00))
DURATION = END_DATE_TIME - START_DATE_TIME
POWER_METER_ID1 = 2213504
POWER_METER_ID2 = 2213505
OBS_ID = POWER_METER_ID1 #May be useful to assign a power meter to each horn
#and label it as East Power Meter or West Power Meter with its serial number.

#------------------
out_time, pm_cw, a1p1_rpg_cw, a2p2_rpg_cw, out_a1p1_wo_cw, out_a2p2_wo_cw = \
    CW_ReadIn_Funcs.zip_rpg_and_pm(RPG_FILE_PATH, PM_FILE_PATH, START_DATE_TIME, DURATION,'W', 
                                   OBS_ID, AVERAGING_TIME)
out_time, pm_cw, a1p2_rpg_cw, a2p1_rpg_cw, out_a1p2_wo_cw, out_a2p1_wo_cw = \
    CW_ReadIn_Funcs.zip_rpg_and_pm(RPG_FILE_PATH, PM_FILE_PATH, START_DATE_TIME, DURATION,'E', 
                                   OBS_ID, AVERAGING_TIME)


gains_a1p1 = a1p1_rpg_cw / pm_cw
gains_a2p2 = a2p2_rpg_cw / pm_cw
#prelim
gains_a1p2 = a1p2_rpg_cw / pm_cw
gains_a2p1 = a2p1_rpg_cw / pm_cw

CorrectedTotPowerWest1 = out_a1p1_wo_cw * pm_cw / a1p1_rpg_cw
CorrectedTotPowerWest2 = out_a2p2_wo_cw * pm_cw / a2p2_rpg_cw
    
fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)
fig.autofmt_xdate(rotation=20)
ax2.set_ylabel(r'$P_{PM}$ [mW]')
ax1.set_ylabel(r'$P_{RPG}$ [arb.]')
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.grid()
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.grid()
ax1.plot(out_time, a1p1_rpg_cw, label=r'W P(l, $\pi$)')
ax1.plot(out_time, a2p2_rpg_cw, label='W P(r, 0)')
ax2.plot(out_time, pm_cw)
ax1.legend(loc='upper right')
fig.tight_layout()
fig.savefig(IMAGE_SAVE_PATH+'WestHornPowersNew.png', dpi=400)
plt.show()

fig, (ax1) = plt.subplots(1,1)
fig.autofmt_xdate(rotation=20)
#ax2.set_ylabel(r'$P_{PM}$ [mW]')
ax1.set_ylabel(r'$P_{CW, E}$ [arb.]')
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.grid()
#ax2.yaxis.set_minor_locator(AutoMinorLocator())
#ax2.grid()
ax1.plot(out_time, a1p2_rpg_cw, label=r'E P(l, 0)')
ax1.plot(out_time, a2p1_rpg_cw, label=r'E P(r, $\pi$)')
#ax2.plot(out_time, pm_cw)
ax1.legend(loc='upper right')
fig.tight_layout()
fig.savefig(IMAGE_SAVE_PATH+'EastHornPowersNew.png', dpi=400)
plt.show()


plt.plot(out_time, CorrectedTotPowerWest1, label=r'W P(l, $\pi$)')
plt.plot(out_time, CorrectedTotPowerWest2, label=r'W P(r, 0)')
plt.legend()
plt.ylabel('Power [arb.]')
plt.tight_layout()
plt.xticks(rotation=20)
plt.grid()
plt.title('Corrected Signal - West')
plt.tight_layout()
plt.savefig(IMAGE_SAVE_PATH+'Corrected_SignalWest.png', dpi=400)
plt.show()

