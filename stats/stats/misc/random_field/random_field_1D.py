'''
Created on 24.06.2011

@author: rrypl
'''

from etsproxy.traits.api import HasTraits, Float, Array, Int, Property, \
    cached_property, Bool, Event, Enum
from math import e
from numpy import dot, transpose, ones, array, eye, linspace, reshape
from numpy.linalg import eig
from numpy.random import shuffle
from scipy.linalg import toeplitz
from scipy.stats import norm, weibull_min
import numpy as np

class RandomField(HasTraits):
    '''Generating a random field by scaling a standardized normal distribution random
     field. The random field array is stored in the property random_field'''

    #Parameters to be set
    lacor = Float(1. , auto_set = False, enter_set = True,
                   desc = 'autocorrelation  of the field', modified = True)
    nsim = Int(1 , auto_set = False, enter_set = True,
                desc = 'No of Fields to be simulated', modified = True)
    mean = Float(0, auto_set = False, enter_set = True,
                  desc = 'mean value', modified = True)
    stdev = Float(1., auto_set = False, enter_set = True,
                   desc = 'standard deviation', modified = True)
    shape = Float(10., auto_set = False, enter_set = True,
                  desc = 'shape for 3 params weibull', modified = True)
    scale = Float(5., auto_set = False, enter_set = True,
                   desc = 'scale at total length L for 3 params weibull', modified = True)
    loc = Float(1., auto_set = False, enter_set = True,
                   desc = 'location for 3 params weibull', modified = True)

    xgrid = Array
    distribution = Enum('Gauss', 'Weibull', modified = True)

    scale_gridpoints = Property(depends_on='scale, lacor, xgrid')
    @cached_property
    def _get_scale_gridpoints(self):
        return self.scale * (self.xgrid[-1] / self.lacor) ** (1. / self.shape)

    non_negative_check = False
    reevaluate = Event
    seed = Bool(False)

    def acor(self, dx, lcorr):
        '''autocorrelation function'''
        return e ** (-(dx / lcorr) ** 2)

    eigenvalues = Property(depends_on = 'lacor')
    @cached_property
    def _get_eigenvalues(self):
        '''evaluates the eigenvalues and eigenvectors of the autocorrelation matrix'''
        #creating a symm. toeplitz matrix with (xgrid, xgrid) data points
        Rdist = toeplitz(self.xgrid, self.xgrid)
        #apply the autocorrelation func. to get the correlation matrix
        R = self.acor(Rdist , self.lacor)
        #evaluate the eigenvalues and eigenvectors of the autocorrelation matrix
        print 'evaluating eigenvalues for random field...'
        eigenvalues = eig(R)
        print 'complete'
        return eigenvalues

    random_field = Property(Array , depends_on = '+modified, reevaluate')
    @cached_property
    def _get_random_field(self):
        if self.seed == True:
            np.random.seed(1540)
        '''simulates the Gaussian random field'''
        #evaluate the eigenvalues and eigenvectors of the autocorrelation matrix
        _lambda, phi = self.eigenvalues
        #simulation points from 0 to 1 with an equidistant step for the LHS
        randsim = linspace(0, 1, len(self.xgrid) + 1) - 0.5 / (len(self.xgrid))
        randsim = randsim[1:]
        #shuffling points for the simulation
        shuffle(randsim)
        #matrix containing standard Gauss distributed random numbers
        xi = transpose(ones((self.nsim, len(self.xgrid))) * array([ norm().ppf(randsim) ]))
        #eigenvalue matrix 
        LAMBDA = eye(len(self.xgrid)) * _lambda
        #cutting out the real part
        ydata = dot(dot(phi, (LAMBDA) ** 0.5), xi).real
        if self.distribution == 'Gauss':
            # scaling the std. distribution
            scaled_ydata = ydata * self.stdev + self.mean
        elif self.distribution == 'Weibull':
            # setting Weibull params
            Pf = norm().cdf(ydata)
            scaled_ydata = weibull_min(self.shape, scale = self.scale_gridpoints, loc = self.loc).ppf(Pf)
        self.reevaluate = False
        rf = reshape(scaled_ydata, len(self.xgrid))
        if self.non_negative_check == True:
            if (rf < 0).any():
                raise ValueError, 'negative value(s) in random field'
        return rf

if __name__ == '__main__':

    from matplotlib import pyplot as p
    rf = RandomField(lacor = 6., xgrid = linspace(0, 100., 300), mean = 4., stdev = 1.5)
    x = rf.xgrid
    rf.distribution = 'Weibull'
    rf.loc = .0
    rf.shape = 10.
    rf.scale = 5.
    p.plot(x, rf.random_field, lw = 2, color = 'black', label = 'Weibull')
    rf.distribution = 'Gauss'
    p.plot(x, rf.random_field, lw = 2, label = 'Gauss')
    p.legend(loc = 'best')
    p.show()
