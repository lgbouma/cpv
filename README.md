# complex rotators

== complex periodic variables
== Stauffer stars
== transient-dipping T Tauris
== "scallops", but at 2 minute cadence (cf. Gunther+2020)
== "transient periodic dip" stars

Young stars are a lot of fun, is all I'm saying.


# install

`conda create --name cpv python=3.12;  conda activate cpv; pip install -e .`

extra dependencies:

`pip install --no-build-isolation aesthetic==0.9`

`pip install --no-build-isolation cdips`

`pip install pygam`

astrobase:
-> `pip install astrobase; pip install numpy==2.2.6` (to fix what the original
pip install breaks)

(this is b/c it depends on pyeebls 0.2.0, which requires numpy<2, but the rest of this
environment is built on numpy 2.2.6 which is incompatible.)

`pip install mpl-scatter-density`
