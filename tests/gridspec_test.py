import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def do_stuff(cell): #just so the plots show up
    ax = plt.subplot(cell)
    ax.plot()
    #ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)


#plt.subplots_adjust(hspace=0.0)

#make outer gridspec
nrows = 2
ncols = 1

#gridspec.GridSpec(nrows, ncols, figure=None, left=None, bottom=None, right=None, top=None, wspace=None, hspace=None, width_ratios=None, height_ratios=None)
outer = gridspec.GridSpec(nrows, ncols, height_ratios = [1, 6], hspace=0)

#make nested gridspecs
gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[0])
gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = outer[1], hspace = .05)

for cell in gs1:
    do_stuff(cell)

for cell in gs2:
    do_stuff(cell)

plt.savefig('temp.png', bbox_inches='tight')


#gridspec = dict(hspace=0.0, height_ratios=[1, 1, 0.4, 3])
#fig, axs = plt.subplots(nrows=4, ncols=1, gridspec_kw=gridspec)
