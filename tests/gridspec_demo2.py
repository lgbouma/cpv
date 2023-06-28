import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def annotate_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)


fig = plt.figure()
fig.suptitle("Controlling spacing around and between subplots")

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

gs1 = GridSpec(3, 3, left=0.05, right=0.48, wspace=0.05)
ax1 = fig.add_subplot(gs1[:-1, :])
ax2 = fig.add_subplot(gs1[-1, :-1])
ax3 = fig.add_subplot(gs1[-1, -1])

gs2 = GridSpec(3, 3, left=0.55, right=0.98, hspace=0.05)
ax4 = fig.add_subplot(gs2[:, :-1])
ax5 = fig.add_subplot(gs2[:-1, -1])
ax6 = fig.add_subplot(gs2[-1, -1])

annotate_axes(fig)

plt.savefig('temp.png', bbox_inches='tight', dpi=400)
