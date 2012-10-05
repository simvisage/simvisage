from material import Material
from etsproxy.traits.api import Float


def H(x):
    return x >= 0


class ACK(HasTraits):
    '''
    ACK model with one matrix breaking strain all over the specimen

    material parameters

    tension strength matrix:          sigma_mu[MPa]
    E-Modulus matrix:                 Em [MPa]
    E-Modulus fibers:                 Ef [MPa]
    reinforcement ratio:              Vf [-]

    program parameters

    plot range:                       sigma_max [MPa]
    '''

    material = Instance(Material)
    sigma_max = Float
    sigma_mu = Float

    def eps_1(self, sigma):
        #sigma<sigma_mu
        eps_c = self.sigma_mu / self.material.Em
        return eps_c

    def eps_2(self, sigma):
        #sigma=sigma_mu
        alpha = self.Em * self.Vm / (self.Ef * self.Vf)
        eps_m_c = (1 + 0.666 * alpha) * self.sigma_mu / self.Em
        return eps_m_c

    def eps_3(self, sigma):
        #sigma>sigma_mu
        Kc = (self.Em * self.Vm + self.Ef * self.Vf)
        sigma_cu = self.eps_1(self.sigma_mu) * Kc
        eps_diff = sigma_cu / self.Ef / self.Vf - self.eps_2(sigma_cu)
        eps_c_u = sigma / self.Ef / self.Vf - eps_diff
        return eps_c_u

    def plot_diagram(self):
        Kc = (self.Em * self.Vm + self.Ef * self.Vf)
        sigma_cu = self.eps_1(self.sigma_mu) * Kc
        eps_list = [0, self.eps_1(sigma_cu)
                    , self.eps_2(sigma_cu), self.eps_3(self.sigma_max)]
        sigma_list = [0, sigma_cu, sigma_cu, self.sigma_max]

        ''''''
        test_eps = [0, self.sigma_max / self.Ef / self.Vf]
        test_sigma = [0, self.sigma_max]
        plt.plot(test_eps, test_sigma)
        ''''''

        plt.plot(eps_list, sigma_list)
        plt.ylabel('$\sigma_c$ in [MPa]', fontsize=16)
        plt.xlabel('$\epsilon_c$ in [-]', fontsize=16)
        plt.title('ACK-Model ')
        plt.show()

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    material = Material(V_f=0.99,
                        E_f=100e3,
                        E_m=30e3,)
    a = ACK(sigma_max=1000.,
            sigma_mu=5.0)

    a.plot_diagram()
