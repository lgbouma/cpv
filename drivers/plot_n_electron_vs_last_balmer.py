"""
The Inglis-Teller equation exists, and it provides one simplified constraint on
the electron density.
"""
import numpy as np
import matplotlib.pyplot as plt

def calculate_n_electron(last_visible_balmer_line: int) -> float:
    """
    Estimate the electron density (n_electron) using the Inglis-Teller equation.

    Args:
        last_visible_balmer_line (int): The principal quantum number of the last visible Balmer line (N).

    Returns:
        float: The estimated electron density in cm^-3.
    """
    # Constants
    a_0 = 5.29177e-9  # Bohr radius in cm

    # Inglis-Teller equation for electron density
    n_electron = (0.027 * a_0**-3) / (last_visible_balmer_line**(15/2))

    return n_electron

# Range of last visible Balmer lines
balmer_lines = np.arange(10, 31)
electron_densities = [calculate_n_electron(n) for n in balmer_lines]

# Plotting
fig, ax = plt.subplots()
ax.plot(balmer_lines, electron_densities, marker='o', linestyle='-')
ax.set_yscale('log')
ax.set_xlabel('Last Visible Balmer Line (N)')
ax.set_ylabel('Electron Density (cm^-3)')
ax.set_title('Electron Density vs. Last Visible Balmer Line')
plt.grid(True, which="both", ls="--")

# Save plot to hard drive
plt.savefig('../results/inglis_teller/electron_density_vs_balmer_line.png')
plt.close()

