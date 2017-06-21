#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Jan 31, 2012 by: rch

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as md
import datetime as dt

import pylab as p

if __name__ == '__main__':
    file_name = 'mroof2_temp_humidity.txt'

    xfmt = md.DateFormatter('%Y-%m-%d %H:%M:%S')

    record = [('n', np.int), ('d', dt.date),
          ('time', dt.datetime),
          ('t', np.float), ('h', np.float) ]

    data = np.loadtxt(file_name, skiprows = 8,
                      dtype = record,
                      usecols = (0, 1, 2, 3, 4))

    dformat = '%d.%m.%Y %H:%M:%S'

    dates = np.array([ dt.datetime.strptime(d[1] + ' ' + d[2], dformat) for d in data ])
    temp = np.array([ d[3] for d in data ], dtype = float)
    mask = temp <= 10.0
    hum = np.array([ d[4] for d in data ])

    dates = dates[mask]
    temp = temp[mask]
    hum = hum[mask]

    ax = p.gca()
    ax.xaxis.set_major_formatter(xfmt)
    ax.grid()
    p.plot(dates, temp, color = 'blue')
    p.ylabel('temperature')
    p.xlabel('date')

    p.twinx()
    p.plot(dates, hum, color = 'red')
    p.ylabel('relative humidity')

    p.show()
