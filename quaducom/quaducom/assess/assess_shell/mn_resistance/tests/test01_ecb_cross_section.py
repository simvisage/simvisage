'''
Created on Nov 21, 2013

@author: rch
'''

from quaducom.assess.assess_shell.mn_resistance  import \
    ECBCrossSection

import numpy as np

def test_ecb_cross_section_mn():
    '''Test the moment and normal force calculated for a cross section.
    '''
    cp = ECBCrossSection(n_layers=3,
                         thickness=0.05,
                         width=0.1
                         )

    assert np.allclose([cp.M, cp.N], [1.14513592334, -22.1303533699])
    cp.n_layers = 5
    assert np.allclose([cp.M, cp.N], [1.29225385264, -6.60917224146])

if __name__ == '__main__':
    test_ecb_cross_section_mn()
