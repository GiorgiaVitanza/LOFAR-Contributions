import sys
import pylab as lb
import numpy as np
import pyrap.tables as pt
import coordinates_mode as co
from scipy import *
import matplotlib.font_manager as font_manager

# Print closure phase vs time/elevation for three selected antennas and channel range/polarisation

pol = 0
start_chan = 0
end_chan = 14
ant1 = 0
ant2 = 1
ant3 = 2
RAdeg = 123.40025  # phase centre RA,DEC in degrees. Needed for elevation calculation
DECdeg = 48.2179139
xaxis = 'time'

# Open table
t = pt.table('/home/giorgia/LOFAR-Contributions/Maintained/sco_174_iter0.ms')

lb.figure(1)
lb.clf()


def phase(antenna1, antenna2):
    t1 = t.query('ANTENNA1= ' + str(antenna1) + ' AND ANTENNA2= ' + str(antenna2))
    print('***FLAG ANALYSIS***')
    datapoint = 0
    noflags = 0
    flag = t1.getcolslice("FLAG", [start_chan, pol], [end_chan, pol])
    phase_data = t1.getcolslice("DATA", [start_chan, pol], [end_chan, pol])
    time = t1.getcol("TIME")
    time = time / (24 * 3600)  # Convert MJD in seconds to days
    timehr = (time - time[0]) * 24  # Hours since observation start
    ant1_cell = t1.getcell("ANTENNA1", 0)
    ant2_cell = t1.getcell("ANTENNA2", 0)
    print('Baseline = ' + str(ant1_cell) + ' ' + str(ant2_cell))
    length = len(flag)
    mphase = []
    elevation = []

    for x in range(len(time)):
        total = 0 + 0j
        count = 0
        tmpflag = flag[x]
        tmpphase = phase_data[x]
        altaz = co.altaz(time[x].item(), RAdeg, DECdeg)
        alt = altaz[0]
        elevation.append(alt)

        for i in range(0, end_chan):
            if tmpflag[i]:
                noflags += 1
                datapoint += 1
            else:
                test = tmpphase[i].item()
                total += test
                count += 1
                datapoint += 1

        if count != 0:
            spectral_mean = total / count
            phase_element = np.arctan2(spectral_mean.imag, spectral_mean.real)
            mphase.append(phase_element)
        else:
            mphase.append(np.nan)

    print(len(mphase))
    print(len(time))
    print(len(elevation))
    print('Total number of datapoints = ' + str(datapoint))
    print('Total number of flags = ' + str(noflags))
    percentage = (float(noflags) / float(datapoint)) * 100
    print('Percentage baseline flagged for correlation ' + str(pol) + ' = ' + str(percentage) + '%')

    return elevation, timehr, mphase


elev, timehr, mphase01 = phase(ant1, ant2)
elev2, timehr2, mphase12 = phase(ant2, ant3)
elev3, timehr3, mphase02 = phase(ant1, ant3)

# The closure phase is reduced modulo 360deg to the range -180 to +180. 
# The closure phase is undefined (NaN) if any of the constituent visibility phases is absent or flagged
min_len = min(len(mphase01), len(mphase12), len(mphase02))
closure = (np.array(mphase01[:min_len]) +
           np.array(mphase12[:min_len]) -
           np.array(mphase02[:min_len])) % (2 * np.pi)

#closure = (array(mphase01) + array(mphase12) - array(mphase02)) % (2 * pi)
closure[np.where(closure > np.pi)] = closure[np.where(closure > np.pi)] - (2 * np.pi)
closure[np.where(closure < -np.pi)] = closure[np.where(closure < -np.pi)] + (2 * np.pi)

titlestring = "Closure Phase Analysis"
ylabelstring = "Closure Phase/Radians"

if xaxis == 'time':
    xlabelstring = "Time Since Obs Start/hours"
    min_len = min(len(timehr), len(closure))
    lb.plot(timehr[:min_len], closure[:min_len], ',')

    #lb.plot(timehr, closure, ',')
elif xaxis == 'elevation':
    xlabelstring = "Elevation Angle/degrees"
    lb.plot(elev, closure, ',')

lb.ylabel(ylabelstring)
lb.xlabel(xlabelstring)
lb.title(titlestring)

lb.show()
