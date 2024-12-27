def get_sens_adjoint(self):
    """
    Adjoint method to get sensitivities for a specific species + plot
    @return: None
    """
    self.solver_adjoint_init()
    sens_df = self.calculate_solver_adjoint_sens()
    self.plot_sens(sens_df, type_f='adjoint_reactions')

    def solver_adjoint_init(self):
        """
        Setup for getting solver adjoint values
        @return:
        """
        self.grid_point = min(range(len(self.f.grid)), key=lambda i: abs(self.f.grid[i] - config.SENS_SPECIES_LOC))
        logger.info(f"using an adjoint method on grid point {self.grid_point}")

        # Index of self.species in the global solution vector
        if self.species == 'lbv':
            # auto target grid point 0 in the solution vector for lbv:
            i_spec = self.f.inlet.n_components + self.f.flame.component_index('velocity')
            self.spec_0 = self.f.velocity[0]
            logger.info(f"running sensitivity analysis for lbv")

        else:
            # specify species and which grid point to target in solution vector:
            i_spec = self.f.inlet.n_components + self.f.flame.component_index(self.species) + self.f.domains[1].n_components*self.grid_point
            self.spec_0 = self.f.X[self.gas.species_index(self.species), self.grid_point]
            logger.info(f"running sensitivity analysis for species {self.species}")

        # ImpingingJet flame has three domains (inlet, flame, surface), but only the flame domain stores information:
        Nvars = sum(D.n_components * D.n_points for D in self.f.domains)

        self.dgdx = np.zeros(Nvars)
        self.dgdx[i_spec] = config.SENS_PERTURBATION


    def calculate_solver_adjoint_sens(self):
        """
        Compute the normalized sensitivities of the species production, taken at the final grid point
        :math:`s_{i, spec}` with respect to the reaction rate constants :math:`k_i`:
        .. math::
            s_{i, spec} = \frac{k_i}{[X]} \frac{d[X]}{dk_i}
        @return:
        """
        def perturb(sim, i, dp):
            sim.gas.set_multiplier(1 + dp, i)

        sens_vals = self.f.solve_adjoint(perturb, self.gas.n_reactions, self.dgdx) / self.spec_0
        return pd.DataFrame(index=self.gas.reaction_equations(), columns=['base_case'], data = sens_vals)