# env: pymc_env (on wh3; mac install having compiler issues)

import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt

def sinusoid(t: np.ndarray,
             amplitude: float,
             freq: float,
             phase: float,
             offset: float) -> np.ndarray:
    """Compute a simple sinusoid.

    Args:
        t: Time array.
        amplitude: Amplitude of the sinusoid.
        freq: Frequency of the sinusoid (cycles per unit time).
        phase: Phase offset in radians.
        offset: Vertical offset.

    Returns:
        The sinusoidal values at each time t.
    """
    return offset + amplitude * np.sin(2.0 * np.pi * freq * t + phase)

df = pd.read_csv('TIC141146667_rel_rv.csv')

# Replace df below with your actual DataFrame
# It should have columns: 'mjd', 'mean_rel_rv', and 'std_rel_rv'
t = df['mjd'].values
rv = df['mean_rel_rv'].values
erv = df['std_rel_rv'].values

# Build the PyMC model
P = 3.930 / 24
with pm.Model() as model:
    amplitude = pm.Normal('amplitude', mu=0.0, sigma=5.0)
    freq = pm.Normal('freq',
                     mu=1.0 / P,
                     sigma=1e-5)
    phase = pm.Uniform('phase', lower=0.0, upper=2.0 * np.pi)
    offset = pm.Normal('offset', mu=np.mean(rv), sigma=5.0 * np.std(rv))

    mu = offset + amplitude * pm.math.sin(2.0 * np.pi * freq * t + phase)

    pm.Normal('obs', mu=mu, sigma=erv, observed=rv)

    # Sample from the posterior
    trace = pm.sample(draws=2000, tune=2000, target_accept=0.95, random_seed=42)

# Extract posterior samples for amplitude, freq, phase, offset
amplitude_samples = trace.posterior['amplitude'].values.flatten()
freq_samples = trace.posterior['freq'].values.flatten()
phase_samples = trace.posterior['phase'].values.flatten()
offset_samples = trace.posterior['offset'].values.flatten()

# Summarize the posterior samples using ArviZ
summary = az.summary(trace, hdi_prob=0.95)

# Display the summary
print(summary)

# 95% confidence upper limit on amplitude
amplitude_997 = np.percentile(amplitude_samples, 99.7)
summary_3sigma = az.summary(trace, hdi_prob=0.997)
print(summary_3sigma)

# Take the median of freq, phase, offset for plotting the sinusoid
freq_median = np.median(freq_samples)
phase_median = np.median(phase_samples)
offset_median = np.median(offset_samples)

# Evaluate the sinusoid at this 95% amplitude
t_model = np.linspace(t.min(), t.max(), 500)
y_model = sinusoid(t_model,
                   amplitude_997,
                   freq_median,
                   phase_median,
                   offset_median)

fig, ax = plt.subplots()
ax.errorbar(df['mjd'],
            df['mean_rel_rv'],
            yerr=df['std_rel_rv'],
            c='k',
            marker='o',
            elinewidth=1,
            capsize=1,
            lw=0,
            mew=0.5,
            markersize=0,
            label='Data')

ax.plot(t_model,
        y_model,
        'r-',
        lw=1,
        label='95% credible max amplitude sinusoid')

ax.legend(loc='best', fontsize='xx-small')
ax.set_xlabel('time')
ax.set_ylabel('rv km/s')
fig.savefig('max_amplitude_sin_pymc.png', dpi=300)
plt.show()

