"""
Calculate bound-free absorption coefficient (cm2 per neutral H atom), for
neutral hydrogen.  As a function of wavelength.

Follows Gray, Eq 8.8

Reproduces Fig 8.2 (great!), but not quite the subsequent figures, perhaps
because of the P_e factor.

Contents:
    κ_Hbf
"""
from astropy import units as u, constants as c
import numpy as np, matplotlib.pyplot as plt
from aesthetic.plot import set_style, savefig

α_0 = 1.0449e-26 # for λ in angstrom
R = 1.0968e-3 # per angstrom; rydberg constant

def χ3(I, n0):
    return I * (
        1 - 1/((n0+3)**2)
    )

def g_bf(λ, n):

    λR = λ*R

    term = (0.3456 / (λR)**(1/3) ) * (
        λR/(n**2) - 0.5
    )

    max_λ = n**2 / R

    result = 1 - term

    result[λ > max_λ] = 0

    return result


def χ(n):
    # returns eV
    return 13.60*(1 - 1/(n**2))

def θ(T):
    # Gray eq8.7, underneath; units inverse eV
    return 5040 / T


def κ_Hbf(λ, n_0, T):
    # Gray Eq 8.8
    # Absorption coefficient in cm2 per neutral H atom, for all continua
    # starting at n_0.

    I = 13.6 # eV

    prefactor = α_0 * λ**3

    term1 = np.log10(np.e) / (2*θ(T)*I) * (
        10**(-χ3(I, n_0) * θ(T)) - 10**(-I*θ(T))
    )

    term0 = 0
    for n in range(n_0, n_0+2):
        term0 += g_bf(λ,n) / (n**3) * 10**(-θ(T)*χ(n))

    parenth = term0 + term1

    return prefactor * parenth


def main():

    # wavelength in angstroms
    λ = np.linspace(0, 2e4, int(1e4))
    T = 10000 # kelvin

    # Gray Fig 8.2, 'alpha_bound-free for hydrogen'
    plt.close("all")
    set_style("clean")
    fig, ax = plt.subplots(figsize=(4,3))
    for n in range(1,6):
        ax.plot(λ, α_0*g_bf(λ,n)*λ**3/(n**5) / (1e-17), label=f'n={n}')
    ax.legend(loc='best', fontsize='small')
    ax.set_ylabel('α_bf(H) [1e-17 cm2/H atom]')
    ax.set_xlabel('λ [angstrom]')
    ax.set_title(f"T = {T} K")

    # convert to opacity
    plt.close("all")
    set_style("clean")
    fig, ax = plt.subplots(figsize=(4,3))
    for n0 in range(1,6):
        ax.plot(λ, κ_Hbf(λ, n0, T) / (1e-17), label=f'n0={n0}')
    ax.set_yscale('log')
    ax.legend(loc='best', fontsize='small')
    ax.set_ylabel('κ_bf(H) [1e-17 cm2/H atom]')
    ax.set_xlabel('λ [angstrom]')
    ax.set_title(f"T = {T} K")

    outpath = '../results/radiative_transfer/bf_opacity_hydrogen.png'
    savefig(fig, outpath)


if __name__ == "__main__":
    main()
