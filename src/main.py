import cantera as ct
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('macosx')
import pandas as pd

def main():

    # Define the gas mixture
    gas = ct.Solution("mech/okafor_2017.yaml")  # Use the GRI-Mech 3.0 reaction mechanism
    gas.TP = 295, ct.one_atm  # Methane-air mixture
    air = "O2:0.21,N2:0.79"
    gas.set_equivalence_ratio(phi=1.0, fuel="NH3:0.7, H2:0.3", oxidizer=air)
    flame = ct.FreeFlame(gas=gas)
    flame.inlet.mdot = 0.35 * gas.density
    # flame.surface.T = 493.5
    # flame.set_initial_guess("equil")

    # Refine grid to improve accuracy
    flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)

    # Solve the flame
    flame.solve(loglevel=1, auto=True)  # Increase loglevel for more output

    # plot temperature:
    # Extract temperature and distance
    temperature = flame.T
    distance = flame.grid  # Grid points (distance in meters)

    # Plot temperature profile
    plt.figure(figsize=(8, 6))
    plt.plot(distance * 1e3, temperature, label="Flame Temperature", color="red")
    plt.xlabel("Distance (mm)")
    plt.ylabel("Temperature (K)")
    plt.title("Temperature Profile of a Flame")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.show()


    # Create a DataFrame to store sensitivity-analysis data
    sens = pd.DataFrame(index=gas.reaction_equations(), columns=["sensitivity"])

    # Use the adjoint method to calculate sensitivities
    sens.sensitivity = flame.get_species_reaction_sensitivities("N2O", 0)

    sens.head(10)
    sens = sens.iloc[(-sens['sensitivity'].abs()).argsort()]
    fig, ax = plt.subplots()

    sens.head(15).plot.barh(ax=ax, legend=None)
    ax.invert_yaxis()  # put the largest sensitivity on top
    ax.set_title("Sensitivities for N2O Using GRI Mech")
    ax.set_xlabel(r"Sensitivity: $\frac{\partial\:\ln X_{N2O}}{\partial\:\ln k}$")
    ax.grid(axis='x')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()