import matplotlib.pyplot as plt, numpy as np
from matplotlib.gridspec import GridSpec
from aesthetic.plot import set_style, savefig

set_style('science')

factor = 0.8
fig = plt.figure(figsize=(factor*7,factor*8.7))
#fig.suptitle("Controlling spacing around and between subplots")

# left, right, top, bottom --- float, optional
# Extent of the subplots as a fraction of figure width or height. Left cannot
# be larger than right, and bottom cannot be larger than top. If not given, the
# values will be inferred from a figure or rcParams at draw time. See also
# GridSpec.get_subplot_params.

# wspace float, optional
# The amount of width reserved for space between subplots, expressed as a
# fraction of the average axis width. If not given, the values will be inferred
# from a import figure or rcParams when necessary. See also
# GridSpec.get_subplot_params.

x0 = 0.03
smx = 0.03
dx = 0.30
hspace = 0.2
left, right = x0, x0+dx
print(left, right)
gs1 = GridSpec(9, 2, left=left, right=right, wspace=0.0, hspace=hspace)
for ix in range(9):
    ax = fig.add_subplot(gs1[ix, 0])
    ax = fig.add_subplot(gs1[ix, 1])

left, right = x0+smx+dx, x0+smx+2*dx
print(left, right)
gs2 = GridSpec(9, 2, left=left, right=right, wspace=0.0, hspace=hspace)
for ix in range(9):
    ax = fig.add_subplot(gs2[ix, 0])
    ax = fig.add_subplot(gs2[ix, 1])

left, right = x0+2*smx+2*dx, x0+2*smx+3*dx
print(left, right)
gs3 = GridSpec(9, 2, left=left, right=right, wspace=0.0, hspace=hspace)
for ix in range(9):
    ax = fig.add_subplot(gs3[ix, 0])
    ax = fig.add_subplot(gs3[ix, 1])

for ix, ax in enumerate(fig.axes):

    ax.text(0.5, 0.5, f"ax{ix}", va="center", ha="center")

    if ix % 2 == 0:
        # or 53 whitespaces otherwise
        spstr = 50*' '
        ax.set_title(spstr+f'{ix} TIC123456789', fontsize=3.5, pad=0.05)

    if ix not in np.array([16,17, 34,35, 52,53]):
        ax.set_xticklabels([])
    if ix % 2 == 1:
        ax.set_yticklabels([])


fs = 'medium'
fig.text(0.5,0.05, r"Phase, φ", fontsize=fs, va='bottom', ha='center')
fig.text(-0.01,0.5, r"Flux [%]", va='center', ha='center', rotation=90, fontsize=fs)


#plt.tight_layout()

savefig(fig, 'temp.png', dpi=400, writepdf=0)
