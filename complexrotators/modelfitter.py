"""
Contents:

ModelFitter
generate_synthetic_data_specific_params

"""
import numpy as np
from numpy import array as nparr
import matplotlib.pyplot as plt
import pymc as pm
import pickle
import os
from os.path import join

from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import differential_evolution, dual_annealing
from scipy.ndimage import gaussian_filter1d

class ModelFitter:
    def __init__(self, xval, yval, flux_arr, err_flux_arr,
                 modelid='1_gaussians', N_samples=1000, N_cores=2, N_chains=2,
                 plotdir=None, overwrite=False, guessdict=None,
                 map_guess_method='raw', velmask=1., flux_arr_nomask=None):
        """
        Initialize the ModelFitter class with data and parameters.

        Parameters:
        - xval: array-like, phase normalized from 0 to 1 (one orbital period)
        - yval: array-like, velocity in units of equatorial velocity
        - flux_arr: 2D array, flux values with shape (len(xval), len(yval))
        - velmask: float, fit only yval < velmask
        - err_flux_arr: 2D array, uncertainties with the same shape as flux_arr
        - modelid: str, specifies the model to use ('1_gaussians', '2_gaussians', etc.)
        - N_samples: int, number of samples for MCMC
        - N_cores: int, number of CPU cores to use
        - N_chains: int, number of chains for MCMC
        - plotdir: str, directory where plots and results will be saved
        - overwrite: bool, whether to overwrite existing results
        - guessdict: dict, keys 'period','K','phi','A','sigma' values the ndarr guesses
        """

        assert velmask >= 1. # don't fit data within line core
        self.VELMASK = velmask
        self.xval = xval
        self.yval = yval
        self.flux_arr = flux_arr
        self.err_flux_arr = err_flux_arr
        self.N_samples = N_samples
        self.N_cores = N_cores
        self.N_chains = N_chains
        self.plotdir = plotdir
        self.overwrite = overwrite
        self.modelid = modelid
        self.map_guess_method = map_guess_method
        self.guessdict = guessdict

        self.flux_arr_nomask = flux_arr_nomask

        # Determine the number of Gaussians (N) from modelid
        self.N = int(self.modelid.split('_')[0])

        # Prepare the results file path
        self.result_file = join(self.plotdir, f'results_{self.modelid}.pkl')

        # Check for existing results
        if os.path.exists(self.result_file) and not overwrite:
            print(f"Results already exist at {self.result_file}. Set overwrite=True to overwrite.")
        else:
            # Create the plot directory if it doesn't exist
            os.makedirs(self.plotdir, exist_ok=True)
            # Run the model fitting
            self.run_ngaussian_model()

    def run_ngaussian_model(self):
        """
        Set up and run the model summing N gaussians using pymc3.
        """

        yval_mask = np.abs(self.yval) >= self.VELMASK  # Shape: (len(self.yval),)

        # Apply the mask to flux_arr and err_flux_arr
        flux_arr_filtered = self.flux_arr[:, yval_mask] # Shape: (len(self.xval), N_filtered_yval)
        err_flux_arr_filtered = self.err_flux_arr[:, yval_mask]  # Same shape as flux_arr_filtered

        yval_filtered = self.yval[yval_mask]  # Shape: (N_filtered_yval,)

        with pm.Model() as model:
            # Lists to hold parameters for each Gaussian component
            period_n = []
            K_n = []
            phi_n = []
            A_n = []
            sigma_n = []

            # Define priors for each Gaussian component
            for n in range(self.N):

                # Prior for the orbital period
                period = pm.Uniform(f'period_{n}', lower=0.7, upper=1.3, initval=1)
                period_n.append(period)

                # K_n: Amplitude scaling of the mean position
                K = pm.Uniform(f'K_{n}', lower=1, upper=5, initval=2)
                K_n.append(K)

                # φ_n: Phase shift
                phi = pm.Uniform(f'phi_{n}', lower=0, upper=2*np.pi, initval=0.5)
                phi_n.append(phi)

                # A_n: Fixed amplitude across all epochs
                A = pm.Uniform(
                    f'A_{n}', lower=0.1, upper=2, initval=1
                )
                A_n.append(A)

                # σ_n: Width of the Gaussian (constant over time)
                #sigma = pm.Uniform(f'sigma_{n}', lower=0.01, upper=1, initval=0.1)
                sigma = pm.TruncatedNormal(f'sigma_{n}', mu=0.1, sigma=0.01,
                                           initval=0.1, lower=0.01)
                sigma_n.append(sigma)

            # Initialize the expected flux model
            flux_model = 0

            # Compute the expected flux for each Gaussian component
            for n in range(self.N):

                omega_n = 2*np.pi / period_n[n]
                mu_n_t = K_n[n] * pm.math.sin(omega_n * self.xval + phi_n[n])  # Shape: (len(self.xval),)
                mu_n_t_reshaped = mu_n_t[:, None]  # Shape: (len(self.xval), 1)
                yval_reshaped = yval_filtered[None, :]  # Shape: (1, N_filtered_yval)

                # Gaussian flux contribution from the nth component
                # Fitted to only the yval >= 1 data...
                gaussian_n = (
                    A_n[n]
                    * pm.math.exp(
                        -0.5 * ((yval_reshaped - mu_n_t_reshaped) / sigma_n[n]) ** 2
                    ) / (sigma_n[n] * np.sqrt(2 * np.pi))
                )

                flux_model += gaussian_n

            # Define the likelihood function (Gaussian)
            flux_obs = pm.Normal('flux_obs', mu=flux_model,
                                 sigma=err_flux_arr_filtered,
                                 observed=flux_arr_filtered)

            # Find the Maximum A Posteriori estimate
            # NOTE This is purely for the purpose of determining the MAP
            # solution to know where to start sampling from.  1_gaussians can
            # use "raw"
            print("Finding MAP estimate...")

            if self.map_guess_method == 'raw':
                map_soln = pm.find_MAP(method='Nelder-Mead')

            elif self.map_guess_method == 'handtuned':

                freeparams = 'period,K,phi,A,sigma'.split(",")
                start = {}
                for n in range(self.N):
                    for param in freeparams:
                        start[f"{param}_{n}"] = self.guessdict[param][n]

                map_soln = pm.find_MAP(start=start)

            self.map_estimate = map_soln

            self.visualize_map_estimate()

            import IPython; IPython.embed()
            assert 0 #FIXME

            # Sample from the posterior distribution
            print("Sampling from posterior...")
            self.trace = pm.sample(
                self.N_samples,
                cores=self.N_cores,
                chains=self.N_chains,
                return_inferencedata=False
            )

            # Cache the results
            with open(self.result_file, 'wb') as f:
                pickle.dump({'map_estimate': self.map_estimate, 'trace': self.trace}, f)

            print(f"Results saved to {self.result_file}")

    def load_results(self):
        """
        Load the cached results from the pickle file.
        """
        if os.path.exists(self.result_file):
            with open(self.result_file, 'rb') as f:
                data = pickle.load(f)
                self.map_estimate = data['map_estimate']
                self.trace = data['trace']
            print("Results loaded successfully.")
        else:
            print("No results found. Please run the model first.")


    def visualize_map_estimate(self):
        """
        Visualize the MAP estimate by plotting the data, model, and residual.
        """

        # Check if MAP estimate is available
        if not hasattr(self, 'map_estimate'):
            print("MAP estimate not found. Please run the model first.")
            return

        # Extract MAP parameter values
        map_params = self.map_estimate

        # Generate higher-resolution grid for xval (phase/time) and yval (velocity)
        xval_model = np.linspace(np.min(self.xval), np.max(self.xval), 100)  # High-res xval
        yval_model = np.linspace(np.min(self.yval), np.max(self.yval), 200)  # High-res yval

        # Initialize the expected flux model over the high-res grid
        flux_model = np.zeros((len(xval_model), len(yval_model)))

        # Extract parameters from MAP estimate
        # Prepare lists for parameters
        period_n = []
        K_n = []
        phi_n = []
        A_n = []
        sigma_n = []

        # For each Gaussian component
        for n in range(self.N):
            period_n.append(map_params[f'period_{n}'])
            K_n.append(map_params[f'K_{n}'])
            phi_n.append(map_params[f'phi_{n}'])
            A_n.append(map_params[f'A_{n}'])
            sigma_n.append(map_params[f'sigma_{n}'])

        # Compute the expected flux for each Gaussian component
        for n in range(self.N):

            omega_n = 2*np.pi / period_n[n]

            mu_n_t = K_n[n] * np.sin(omega_n * xval_model + phi_n[n])  # Shape: (100,)
            mu_n_t_reshaped = mu_n_t[:, None]  # Shape: (100, 1)
            yval_reshaped = yval_model[None, :]  # Shape: (1, 200)

            t_tra = phi_n[n] / (2*np.pi)
            t_sec_ecl = t_tra + period_n[n]/2
            half_tdur = 0.05

            # Mask velocities behind or in front of the star
            mask0 = (np.abs(yval_reshaped) <= 1)

            # Mask for times outside the secondary eclipse
            phase_model = (
                (xval_model - t_tra) / period_n[n]
                -
                np.floor((xval_model - t_tra) / period_n[n])
            )
            mask1 = (
                (phase_model[:, None] > 0.5-half_tdur) &
                (phase_model[:, None] < 0.5+half_tdur)
            )

            # Combine masks
            mask = ~( mask0 & mask1 )

            # Gaussian flux contribution from the nth component
            gaussian_n = A_n[n] * np.exp(
                -0.5 * ((yval_reshaped - mu_n_t_reshaped) / sigma_n[n]) ** 2
            ) / (sigma_n[n] * np.sqrt(2 * np.pi))

            # Apply mask
            mask_gaussian_n = gaussian_n * mask

            # Sum the contributions
            flux_model += mask_gaussian_n

        # Interpolate the model flux onto the data grid
        from scipy.interpolate import RegularGridInterpolator

        # Create the interpolator
        interpolator = RegularGridInterpolator((xval_model, yval_model), flux_model)

        # Create meshgrid of data points
        # xi.shape == (len(self.xval), len(self.yval))
        xi, yi = np.meshgrid(self.xval, self.yval, indexing='ij')

        # Flatten the grid points
        xi_flat = xi.flatten()
        yi_flat = yi.flatten()
        points_to_interpolate = np.column_stack((xi_flat, yi_flat))

        # Evaluate the interpolator
        flux_model_on_data_grid_flat = interpolator(points_to_interpolate)
        flux_model_on_data_grid = flux_model_on_data_grid_flat.reshape(xi.shape)

        # Compute the residuals
        residuals = self.flux_arr_nomask - flux_model_on_data_grid


        #plt.close("all")
        #f_50 = np.nanpercentile(
        #    residuals[:, (np.abs(self.yval) < 1) ], 50, axis=0
        #)
        #fn = lambda x: gaussian_filter1d(x, sigma=5)
        #plt.plot(self.yval[np.abs(self.yval) < 1], fn(f_50), c='k', lw=0.5)
        #plt.show()
        #import IPython; IPython.embed()
        #assert 0


        # Plotting
        fig, axs = plt.subplots(1, 3, figsize=(18, 6))

        # Plot the data
        ax = axs[0]
        cmap = plt.get_cmap('Greys')
        c = ax.pcolor(self.xval, self.yval, self.flux_arr_nomask.T, cmap=cmap,
                      shading='auto', rasterized=True)
        cb = fig.colorbar(c, ax=ax)
        cb.set_label('Flux (Data)')
        ax.set_xlabel('Time (P)')
        ax.set_ylabel('Δv / v_eq')
        ax.set_title('Data (only Δv/veq>1 fitted)')

        # Plot the model
        ax = axs[1]
        c = ax.pcolor(xval_model, yval_model, flux_model.T, cmap=cmap,
                      shading='auto', rasterized=True)
        cb = fig.colorbar(c, ax=ax)
        cb.set_label('Flux (Model)')
        ax.set_xlabel('Time (P)')
        ax.set_title('Model*mask')

        # Plot the residuals
        ax = axs[2]
        c = ax.pcolor(self.xval, self.yval, residuals.T, cmap='RdBu_r',
                      shading='auto', rasterized=True)
        cb = fig.colorbar(c, ax=ax)
        cb.set_label('Flux (Residuals)')
        ax.set_xlabel('Time (P)')
        ax.set_title('Residuals (Data - Model)')

        plt.tight_layout()
        savpath = join(self.plotdir, f"{self.modelid}_map_viz.png")
        plt.savefig(savpath, bbox_inches='tight')
        print(f"wrote {savpath}")


def generate_synthetic_data_specific_params(K_n_true, phi_n_true, sigma_n_true, A_n_true,
                                            num_time_points,
                                            period_true=1.0, num_velocity_points=100,
                                            noise_level=0.1, random_seed=42):
    """
    Generate synthetic data using specific true parameter values.
    """
    np.random.seed(random_seed)
    N = len(K_n_true)

    # Generate xval and yval
    xval = np.linspace(0, 1.3*period_true, num_time_points)
    yval = np.linspace(-5, 5, num_velocity_points)

    # Compute omega_n
    omega_n_true = 2 * np.pi / period_true

    # Initialize flux_arr
    flux_arr = np.zeros((num_time_points, num_velocity_points))

    # Compute flux_arr
    for n in range(N):
        mu_n_t = K_n_true[n] * np.sin(omega_n_true * xval + phi_n_true[n])
        mu_n_t_reshaped = mu_n_t[:, None]
        yval_reshaped = yval[None, :]
        sigma_n = sigma_n_true[n]

        gaussian_n = A_n_true[n] * np.exp(
            -0.5 * ((yval_reshaped - mu_n_t_reshaped) / sigma_n) ** 2
        ) / (sigma_n * np.sqrt(2 * np.pi))

        flux_arr += gaussian_n

    # Add noise
    err_flux_arr = noise_level * np.ones_like(flux_arr)
    noise = np.random.normal(0, noise_level, size=flux_arr.shape)
    flux_arr_noisy = flux_arr + noise

    return xval, yval, flux_arr_noisy, err_flux_arr



