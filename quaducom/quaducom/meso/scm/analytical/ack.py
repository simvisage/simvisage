from .material import Material
from etsproxy.traits.api import HasTraits, Float, Instance, Property


def H(x):
    return x >= 0.0


class ACK(HasTraits):
    '''
    Stochastic cracking model due to Aveston, Cooper and Kelly
    assumed deterministic matrix strength
    '''

    material = Instance(Material)
    def _material_default(self):
        return Material()

    sigma_max = Float
    sigma_mu = Float

    eps1 = Property()
    def _get_eps1(self):
        #sigma<sigma_mu
        return self.sigma_mu / self.material.E_m

    eps2 = Property()
    def _get_eps2(self):
        #sigma=sigma_mu
        Ec = self.material.E_c
        return (Ec/self.material.E_m - .334 * self.material.V_m) * self.sigma_mu /  self.material.E_f / self.material.V_f

    eps3 = Property()
    def _get_eps3(self):
        #sigma>sigma_mu
        return self.sigma_max / (self.material.E_f * self.material.V_f) - .334 * self.material.V_m * self.sigma_mu /  self.material.E_f / self.material.V_f

    def plot_diagram(self):
        Ec = self.material.E_c
        sigma_cu = self.sigma_mu * Ec / self.material.E_m
        print('ecxact = ', sigma_cu - self.eps2 * self.material.E_f * self.material.V_f)
        print(.334 * self.sigma_mu * self.material.V_m)
        eps_list = [0, self.eps1, self.eps2, self.eps3]
        sigma_list = [0, sigma_cu, sigma_cu, self.sigma_max]

        reinf_eps = [0, self.sigma_max / self.material.E_f / self.material.V_f]
        reinf_sigma = [0, self.sigma_max]
        plt.plot(reinf_eps, reinf_sigma)

        plt.plot(eps_list, sigma_list)
        plt.ylabel('$\sigma_c$ in [MPa]', fontsize=16)
        plt.xlabel('$\epsilon_c$ in [-]', fontsize=16)
        plt.title('ACK-Model')
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
