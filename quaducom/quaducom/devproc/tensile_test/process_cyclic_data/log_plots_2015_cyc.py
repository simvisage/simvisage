'''
Created on Jun 15, 2015

@author: rch
'''

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

series_label = ['10e+6', '10e+5', '10e+4', '10e+3', '10e+0']
series_values = [5, 4, 3, 2, 1, 0]

colors = ['black', 'green', 'blue', 'red', 'orange', 'gray']

markers = []
for m in Line2D.markers:
    try:
        if len(m) == 1 and m != ' ':
            markers.append(m)
    except TypeError:
        pass

styles = markers + [
    r'$\lambda$',
    r'$\bowtie$',
    r'$\circlearrowleft$',
    r'$\clubsuit$',
    r'$\checkmark$',
    r'$\times$']

powers = np.array(
    [1.00E+00,    1.00E+03,   1.00E+04,  1.00E+05,   1.00E+06], dtype='float')
upper_strain = np.array([
    [0.0076368,    0.008164,    0.008664,   0.0088204, 0.0092032],
    [0.0080936,    0.0087892,    0.0091796,    0.0096092, 0],
    [0.0075392,    0.0082032,    0.0086016,    0,    0],
    [0.008078,    0.0096952,    0,    0,    0],
    [0.0075508,    0,    0,    0,    0]], dtype='float')
lower_strain = np.array([
    [0.0057852,    0.0063476,    0.0067032,    0.0070664,    0.0073084],
    [0.0061796,    0.0069608,    0.0073084,    0.0076876,    0],
    [0.0058632,    0.0061916,    0.0066444,    0,    0],
    [0.0061756,    0.0078204,    0,    0,    0],
    [0.0061252,    0,    0,    0,    0]], dtype='float')
stiffness = np.array([
    [893441.8205876,    912368.730731117,    843684.656772746,
        943713.476339795,    876164.370909858],
    [868905.923458725,    908517.180321593,
        885649.082139803,    862928.353975854, 0],
    [988215.99045346,    824320.50357924,    846733.375485388, 0, 0],
    [872150.802933138,    884990.232291445, 0, 0, 0],
    [803538.860830528, 0, 0, 0, 0]], dtype='float')
ultimate_stress = np.array([16.594, 16.344,   19.359,  18.453,  17.516],
                           dtype='float')

for i, (label, n_p, c, s) in enumerate(zip(series_label,
                                           series_values, colors, styles)):
    plt.semilogx(
        powers[:n_p], upper_strain[i, :n_p], color=c, marker=s, label=label)
    plt.semilogx(powers[:n_p], lower_strain[i, :n_p], color=c, marker=s)

plt.legend()
plt.show()
