'''
Created on Aug 23, 2012

@author: rch
'''

import sympy as sp

a_, b_, sig_tt_, eps_max_, eps_ = sp.symbols('a,b,sig_tt,eps_max, eps')

sig = -a_ + sp.sqrt(a_ ** 2 + b_ * eps_)  


ex = a_ * eps_ - b_

print sp.solve(ex, b_)
print sp.solve(sig - sig_tt_, b_)
print sp.simplify(sig ** 2)
