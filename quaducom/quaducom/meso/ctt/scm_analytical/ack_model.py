from material import Material
from etsproxy.traits.api import HasTraits, Float, Instance


def H(x):
    return x >= 0


class ACK(HasTraits):
    '''
    ACK model with one matrix breaking strain all over the specimen

    material parameters

    tension strength matrix:          sigma_mu[MPa]
    E-Modulus matrix:                 E_m [MPa]
    E-Modulus fibers:                 E_f [MPa]
    reinforcement ratio:              V_f [-]

    program parameters

    plot range:                       sigma_max [MPa]
    '''

    material = Instance(Material)
    sigma_max = Float
    sigma_mu = Float

    def eps_1(self, sigma):
        #sigma<sigma_mu
        eps_c = self.sigma_mu / self.material.E_m
        return eps_c

    def eps_2(self, sigma):
        #sigma=sigma_mu
        alpha = self.material.E_m * self.material.V_m / (self.material.E_f
                                                       * self.material.V_f)
        eps_m_c = (1 + 0.666 * alpha) * self.sigma_mu / self.material.E_m
        return eps_m_c

    def eps_3(self, sigma):
        #sigma>sigma_mu
        K_c = (self.material.E_m * self.material.V_m
               + self.material.E_f * self.material.V_f)
        sigma_cu = self.eps_1(self.sigma_mu) * K_c
        eps_diff = sigma_cu / (self.material.E_f
                               * self.material.V_f) - self.eps_2(sigma_cu)
        eps_c_u = sigma / self.material.E_f / self.material.V_f - eps_diff
        return eps_c_u

    def plot_diagram(self):
        K_c = (self.material.E_m * self.material.V_m
              + self.material.E_f * self.material.V_f)
        sigma_cu = self.eps_1(self.sigma_mu) * K_c
        eps_list = [0, self.eps_1(sigma_cu)
                    , self.eps_2(sigma_cu), self.eps_3(self.sigma_max)]
        sigma_list = [0, sigma_cu, sigma_cu, self.sigma_max]

        ''''''
        test_eps = [0, self.sigma_max / self.material.E_f / self.material.V_f]
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
    material = Material(V_f=0.3,
                        E_f=100e3,
                        E_m=30e3)
    a = ACK(material=material,
            sigma_max=10.,
            sigma_mu=3.0
            )

    a.plot_diagram()
