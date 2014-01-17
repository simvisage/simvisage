'''
Created on Nov 21, 2013

@author: rch
'''

from quaducom.assess.assess_shell.mn_resistance  import \
    ECBCrossSection

import numpy as np

def test_ecb_cross_section_mn():
    '''Test the mappings including neighbors and connectivity.
    '''
    cp = ECBCrossSection(n_layers=3,
                         thickness=0.05,
                         width=0.1
                         )

    assert np.allclose([cp.M, cp.N], [1.14513592334, -22.1303533699])

if __name__ == '__main__':
    test_ecb_cross_section_mn()
