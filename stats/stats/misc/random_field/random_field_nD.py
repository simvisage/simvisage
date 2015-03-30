'''
Created on Jan 25, 2015
implements the 3D random field
@author: rostislavrypl
'''

from etsproxy.traits.api import HasTraits, Float, Property, \
                                cached_property, Array, Enum, \
                                Event, Bool, List
from scipy.linalg import eigh, eig
from scipy.linalg import toeplitz
import numpy as np
from math import e
from scipy.stats import norm, weibull_min


class RandomField(HasTraits):
    '''
    This class implements a 3D random field on a regular grid
    and allows for interpolation using the EOLE method
    '''
    lacor_arr = Array(Float, modified=True) #(nD,1) array of autocorrelation lengths
    nDgrid = List(Array, modified=True) # list of nD entries: each entry is an array of points in the part. dimension
    non_negative_check = False
    reevaluate = Event
    seed = Bool(False)
    distr_type = Enum('Gauss', 'Weibull', modified=True)
    stdev = Float(1.0, modified=True)
    mean = Float(0.0, modified=True)
    shape = Float(5.0, modified=True)
    scale = Float(1.0, modified=True)
    loc = Float(0.0, modified=True)
    
    def acor(self, dx, lacor):
        '''autocorrelation function'''
        C = e ** (-(dx / lacor) ** 2)
        return C

    eigenvalues = Property(depends_on='+modified')
    @cached_property
    def _get_eigenvalues(self):
        '''evaluates the eigenvalues and eigenvectors of the covariance matrix'''
        # creating distances from the first coordinate
        for i, grid_i in enumerate(self.nDgrid):
            self.nDgrid[i] -= grid_i[0]
        # creating a symm. toeplitz matrix with (xgrid, xgrid) data points
        coords_lst = [toeplitz(grid_i) for grid_i in self.nDgrid]
        # apply the autocorrelation func. on the coord matrices to obtain the covariance matrices
        C_matrices = [self.acor(coords_i, self.lacor_arr[i]) for i, coords_i in enumerate(coords_lst)]
        # evaluate the eigenvalues and eigenvectors of the autocorrelation matrices
        eigen_lst = []
        for i, C_i in enumerate(C_matrices):
            print 'evaluating eigenvalues for dimension ' + str(i+1)
            lambda_i, Phi_i = eigh(C_i)
            # truncate the eigenvalues at 99% of tr(C)
            truncation_limit = 0.99 * np.trace(C_i)
            argsort = np.argsort(lambda_i)
            cum_sum_lambda = np.cumsum(np.sort(lambda_i)[::-1])
            idx_trunc = int(np.sum(cum_sum_lambda < truncation_limit)) 
            eigen_lst.append([lambda_i[argsort[::-1]][:idx_trunc], Phi_i[:, argsort[::-1]][:,:idx_trunc]])
        print 'complete'
        
        Lambda_C = 1.0
        Phi_C = 1.0
        for lambda_i, Phi_i in eigen_lst:
            Lambda_i = np.diag(lambda_i)
            Lambda_C = np.kron(Lambda_C, Lambda_i)
            Phi_C = np.kron(Phi_C, Phi_i)
            
        Lambda_1D = Lambda_C.diagonal() # take the diagonal of the Lambda matrix
        argsort = Lambda_1D.argsort()[::-1] # indices of Lambda_1D sorted in descending order
        Lambda_C_sorted = np.diag(Lambda_1D[argsort]) # diagonal matrix with sorted entries
        Phi_C_sorted = Phi_C[:,argsort] # sort the columns of phi_C according to the eigenvalues
        return Lambda_C_sorted, Phi_C_sorted
    
    random_field = Property(Array, depends_on='+modified')
    @cached_property
    def _get_random_field(self):
        if self.seed == True:
            np.random.seed(141)
        '''simulates the Gaussian random field'''
        # evaluate the eigenvalues and eigenvectors of the autocorrelation matrix
        Lambda_C, phi_C = self.eigenvalues
        # simulation points from 0 to 1 with an equidistant step for the LHS
        npts = Lambda_C.shape[0]
        randsim = np.linspace(0, 1, npts + 1) - 0.5 / npts
        randsim = randsim[1:]
        # shuffling points for the simulation
        np.random.shuffle(randsim)
        # matrix containing standard Gauss distributed random numbers
        xi = norm().ppf(randsim)
        # eigenvalue matrix
        # cutting out the real part
        ydata = np.dot(np.dot(phi_C, (Lambda_C) ** 0.5), xi)
        if self.distr_type == 'Gauss':
            # scaling the std. distribution
            scaled_ydata = ydata * self.stdev + self.mean
        elif self.distr_type == 'Weibull':
            # setting Weibull params
            Pf = norm().cdf(ydata)
            scaled_ydata = weibull_min(self.shape, scale=self.scale, loc=self.loc).ppf(Pf)
        shape = tuple([len(grid_i) for grid_i in self.nDgrid])
        rf = np.reshape(scaled_ydata, shape)
        return rf
    
    def plot_rf(self):
        if len(self.nDgrid) == 1:
            ''' 1D plot '''
            import matplotlib
            matplotlib.use('WxAgg')
            import matplotlib.pyplot as plt
            plt.plot(self.nDgrid[0], self.random_field)
            plt.show()
        elif len(self.nDgrid) == 2:
            ''' 2D plot '''
            import os
            os.environ['ETS_TOOLKIT'] = 'qt4'
            os.environ['QT_API'] = 'pyqt'
            from mayavi import mlab
            rand_field_2D = self.random_field
            x, y = self.nDgrid
            mlab.surf(x, y, rand_field_2D)
            mlab.show()
            mlab.close()
        elif len(self.nDgrid) == 3:
            ''' 3D plot '''
            import os
            os.environ['ETS_TOOLKIT'] = 'qt4'
            os.environ['QT_API'] = 'pyqt'
            from mayavi import mlab
            rand_field_3D = self.random_field
            x, y, z = self.nDgrid
            #x, y, z = numpy.ogrid[-5:5:64j, -5:5:64j, -5:5:64j]
            mlab.contour3d(rand_field_3D, contours=7, transparent=True)
            mlab.show()

if __name__ == '__main__':
    example1D = False
    example2D = True
    example3D = False
    
    if example1D is True:
        rf = RandomField(distr_type='Gauss',
                         lacor_arr=np.array([2.0]),
                         nDgrid=[np.linspace(0.0, 30., 300)]
                         )
        rf.plot_rf()
    
    if example2D is True:
        rf = RandomField(distr_type='Gauss',
                         lacor_arr=np.array([2.0, 1.2]),
                    nDgrid=[np.linspace(0.0, 30., 100),
                            np.linspace(0.0, 30., 100)]
                    )
        rf.plot_rf()
        
    if example3D is True:
        rf = RandomField(distr_type='Gauss',
                         lacor_arr=np.array([2.0, 10.0, 3.0]),
                    nDgrid=[np.linspace(0.0, 50., 30),
                            np.linspace(0.0, 50., 25),
                            np.linspace(0.0, 70., 20)]
                    )
        rf.plot_rf()