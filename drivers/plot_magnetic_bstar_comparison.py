import os
import complexrotators.plotting as rp
from complexrotators.paths import RESULTSDIR

def main():

    PLOTDIR = os.path.join(RESULTSDIR, 'magnetic_bstar_comparison')
    if not os.path.exists(PLOTDIR):
        os.mkdir(PLOTDIR)

    rp.plot_magnetic_bstar_comparison(PLOTDIR, selfn='simpleB')
    rp.plot_magnetic_bstar_comparison(PLOTDIR, selfn='complexB')
    rp.plot_magnetic_bstar_comparison(PLOTDIR, selfn='complexM')
    assert 0
    rp.plot_magnetic_bstar_comparison(PLOTDIR, selfn='hd64740')
    rp.plot_magnetic_bstar_comparison(PLOTDIR, selfn='hd37776')
    assert 0
    rp.plot_magnetic_bstar_comparison(PLOTDIR, selfn='all')
    rp.plot_magnetic_bstar_comparison(PLOTDIR, selfn='both')


    #NOTE: not implemented, b/c i would have probably autolabelled it an EB
    #rp.plot_magnetic_bstar_comparison(PLOTDIR, selfn='sigmaoriE')

if __name__ == "__main__":
    main()
