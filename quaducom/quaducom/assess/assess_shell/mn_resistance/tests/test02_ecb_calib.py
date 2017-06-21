'''
Created on Nov 21, 2013

@author: rch
'''

from quaducom.assess.assess_shell.mn_resistance  import \
    ECBLCalib

import numpy as np

def test_ecb_law_calib():
    '''Test the calibrated crack bridge law.
    '''
    calib = ECBLCalib()
    assert np.allclose(calib.calibrated_ecb_law.mfn.ydata[:10],
                       np.array([    0.        , 51.09432702, 98.10817152, 142.21106285,
                                 183.90544436, 223.50840093, 261.24757242, 297.2982659 ,
                                 331.801653  , 364.87500021        ],
                                dtype=float
                                ))

if __name__ == '__main__':
    test_ecb_law_calib()
