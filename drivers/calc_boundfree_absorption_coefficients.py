"""
Calculate bound-free absorption coefficient (cm2 per neutral H atom), for
neutral hydrogen.  As a function of wavelength.

Follows Gray, Eq 8.7

[[note: doesn't actually reproduce Gray's plots.  Gray 8.8 migth be the better
one to implement]]
"""
from astropy import units as u, constants as c
import numpy as np, matplotlib.pyplot as plt

# 100 nm to 2 microns... -> UV to NIR
λ = np.linspace(0, 2e4, int(1e4))

α_0 = 1.0449e-26 # for λ in angstrom

α_bf = np.zeros(len(λ))

T = 5000*u.K

for n in range(2,5):

    # eq 8.5
    R = 1.0968*1e-3 / (u.angstrom) # Rydberg constant
    λR = (λ*u.angstrom * R).cgs.value

    this = 0
    for _n in range(n,100):

        # eq 8.3
        χ = 13.60 * (1 - 1/_n**2) * u.eV
        kT = (c.k_B * T).to(u.eV)
        expfactor = (-χ / kT).cgs.value

        g_bf = (
            1 -
            (0.3456 / ((λR)**(1/3))) * (
                (λR)/_n**2 - 0.5
            )
        )

        this += (
            α_0 * g_bf * λ**3 / (_n**3) * np.exp(expfactor)
        ) / (1e-17)

    plt.plot(λ, this, label=f"{n}")

    α_bf += this

#plt.plot(λ, α_bf, label='total')

plt.legend(loc='best')
plt.xlabel('wavelen [A]')
plt.ylabel('bf absorption coeff (cm2/H atom)')
#plt.yscale('log')
#plt.xscale('log')
plt.savefig('../results/radiative_transfer/bf_hydrogen.png', dpi=300,
            bbox_inches='tight')
