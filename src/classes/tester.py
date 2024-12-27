import numpy as np
def get_flame_speed_reaction_sensitivities(self):
    r"""
    Compute the normalized sensitivities of the laminar flame speed
    :math:`S_u` with respect to the reaction rate constants :math:`k_i`:

    .. math::

        s_i = \frac{k_i}{S_u} \frac{dS_u}{dk_i}
    """

    def g(sim):
        return sim.velocity[0]

    Nvars = sum(D.n_components * D.n_points for D in self.domains)

    # Index of u[0] in the global solution vector
    i_Su = self.inlet.n_components + self.flame.component_index('velocity')

    dgdx = np.zeros(Nvars)
    dgdx[i_Su] = 1

    Su0 = g(self)

    def perturb(sim, i, dp):
        sim.gas.set_multiplier(1 + dp, i)

    return self.solve_adjoint(perturb, self.gas.n_reactions, dgdx) / Su0

def get_species_reaction_sensitivities(self, species:str, ix:int):
    r"""
    Compute the normalized sensitivities of a species
    :math:`S_u` with respect to the reaction rate constants :math:`k_i`:

    .. math::
        s_i = \frac{k_i}{S_u} \frac{dS_u}{dk_i}
    :param species:
        Species of interest for sensitivity calculation (i.e "NO")
    :param ix:
        grid point index (0 is the first, -1 is the last point in the domain)
    """
    def g(sim):
        return self.flame.X[self.gas.species_index(species), ix]

    Nvars = sum(D.n_components * D.n_points for D in self.flame.domains)
    # Nvars = sum(D.n_components * D.n_points for D in self.f.domains)

    # Index of u[0] in the global solution vector modified
    # i_Su = self.inlet.n_components + self.flame.component_index('velocity')
    # i_Su = self.f.inlet.n_components + self.f.flame.component_index(self.species) + self.f.domains[
    #     1].n_components * self.grid_point
    i_Su = self.inlet.n_components + self.flame.component_index(species) + self.flame.domains[1].n_components * ix

    dgdx = np.zeros(Nvars)
    dgdx[i_Su] = 1

    Su0 = g(self)

    def perturb(sim, i, dp):
        sim.gas.set_multiplier(1 + dp, i)

    # return self.solve_adjoint(perturb, self.gas.n_reactions, dgdx) / Su0
    return self.solve_adjoint(perturb, self.gas.n_reactions, dgdx) / Su0
