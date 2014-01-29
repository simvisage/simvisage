'''
Created on Nov 21, 2013

@author: rch
'''

from quaducom.assess.assess_shell.mn_resistance  import \
    ECBCrossSection, ECBReinfTexUniform, ECBCrossSectionGeo

import numpy as np

def test_ecb_cross_section_mn():
    '''Test the moment and normal force calculated for a cross section.
    '''
    cp = ECBCrossSection(components=[ECBReinfTexUniform(n_layers=3),
                                     ECBCrossSectionGeo(width=0.1, n_cj=20)],
                         eps_lo=0.014,
                         eps_up= -0.0033,
                         height=0.05,
                         )

    print 'eps_ti_arr', cp.components_with_state[0].eps_ti_arr
    print 'sig_ti_arr', cp.components_with_state[0].sig_ti_arr
    print 'f_ti_arr', cp.components_with_state[0].f_ti_arr

    print 'eps_cj_arr', cp.components_with_state[1].eps_ti_arr
    print 'sig_cj_arr', cp.components_with_state[1].sig_ti_arr
    print 'f_cj_arr', cp.components_with_state[1].f_ti_arr

    print cp.M, cp.N
    assert np.allclose([cp.M, cp.N], [1.14513592334, -22.1303533699])
    cp.n_layers = 5
    assert np.allclose([cp.M, cp.N], [1.29225385264, -6.60917224146])
    cp.eps_lo = 0.010
    assert np.allclose([cp.M, cp.N], [1.52939155655, -28.4691640432])

if __name__ == '__main__':
    test_ecb_cross_section_mn()
