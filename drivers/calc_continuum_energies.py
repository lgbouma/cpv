# hydrogen

import numpy as np
from astropy import units as u, constants as c

from cdips_followup.spectools import vac_to_air

n_1 = 3 # "Paschen" continua

n_1_dict = {
    1: 'Lyman',
    2: 'Balmer',
    3: 'Paschen/Bohr',
    4: 'Brackett',
    5: 'Pfund',
    6: 'Humphreys',
    7: 'seven'
}

for n_1, name in n_1_dict.items():
    print(42*'-')
    print(f'n_1: {n_1} ({name})')
    for n_2 in range(n_1+1,20):
        ΔE = (13.6*u.eV) * (1/n_1**2 - 1/n_2**2)
        λ = c.h * c.c / ΔE
        txt = (
            f'{str(n_1).zfill(2)}->{str(n_2).zfill(2)}: ΔE = {ΔE.to(u.eV).value:.2f}eV, '
            f'λvac = {λ.to(u.Angstrom).value:.1f} A, '
            f'λair = {vac_to_air(λ).to(u.Angstrom).value:.1f} A'
        )
        print(txt)
