"""
Quick implementation of Euler inversion
"""

import concurrent.futures
import numpy as np
import scipy as sp
import verde as vd
import xarray as xr


class EulerDeconvolution:

    def __init__(self, structural_index):
        self.structural_index = structural_index

    def fit(self, coordinates, data):
        """
        Perform Euler Deconvolution on a single window spanning the entire data
        """
        east, north, up = coordinates
        field, deriv_east, deriv_north, deriv_up = data
        jacobian = np.empty((field.size, 4))
        jacobian[:, 0] = deriv_east
        jacobian[:, 1] = deriv_north
        jacobian[:, 2] = deriv_up
        jacobian[:, 3] = self.structural_index
        pseudo_data = (
            east * deriv_east
            + north * deriv_north
            + up * deriv_up
            + self.structural_index * field
        )
        hessian = jacobian.T @ jacobian
        cofactor = sp.linalg.inv(hessian)
        parameters = cofactor @ jacobian.T @ pseudo_data
        pseudo_residuals = pseudo_data - jacobian @ parameters
        chi_squared = np.sum(pseudo_residuals**2) / (pseudo_data.size - parameters.size)
        self.cofactor_ = cofactor
        self.covariance_ = chi_squared * cofactor
        self.location_ = parameters[:3]
        self.base_level_ = parameters[-1]
        return self

    def fit_grid(
        self,
        grid,
        data_names=("field", "deriv_east", "deriv_north", "deriv_up"),
        coordinate_names=("easting", "northing", "height"),
    ):
        """
        Perform Euler Deconvolution on data on a regular grid
        """
        shape = grid[data_names[0]].shape
        table = vd.grid_to_table(grid)
        coordinates = [table[n] for n in coordinate_names]
        data = [table[n] for n in data_names]
        self.fit(coordinates, data)
        return self


class EulerDeconvolutionWindowed:

    def __init__(self, window_size, window_step, structural_index):
        self.structural_index = structural_index
        self.window_size = window_size
        self.window_step = window_step

    def fit(self, coordinates, data):
        _, windows = vd.rolling_window(
            coordinates,
            size=self.window_size,
            spacing=self.window_step,
            adjust="region",
        )
        self.solutions_ = []
        for window in windows.ravel():
            ed = EulerDeconvolution(self.structural_index)
            window_coordinates = [c[window[0]] for c in coordinates]
            ed.fit(
                window_coordinates,
                [d[window[0]] for d in data],
            )
            small_variance = (
                np.sqrt(ed.covariance_[2, 2]) < np.abs(ed.location_[2]) / 10
            )
            inside_window = vd.inside(ed.location_, vd.get_region(window_coordinates))
            if small_variance and inside_window:
                self.solutions_.append(ed)
        self.locations_ = np.transpose([ed.location_ for ed in self.solutions_])
        self.base_levels_ = np.array([ed.base_level_ for ed in self.solutions_])
        return self

    def fit_grid(
        self,
        grid,
        data_names=("field", "deriv_east", "deriv_north", "deriv_up"),
        coordinate_names=("easting", "northing", "height"),
    ):
        """
        Perform Euler Deconvolution on data on a regular grid
        """
        shape = grid[data_names[0]].shape
        table = vd.grid_to_table(grid)
        coordinates = [table[n] for n in coordinate_names]
        data = [table[n] for n in data_names]
        self.fit(coordinates, data)
        return self


class EulerInversion:

    def __init__(
        self, structural_index, max_iterations=20, tol=0.1, euler_misfit_balance=0.1
    ):
        self.structural_index = structural_index
        self.max_iterations = max_iterations
        self.tol = tol
        self.euler_misfit_balance = euler_misfit_balance

    def fit_grid(
        self,
        grid,
        weights=(1, 0.1, 0.1, 0.05),
        data_names=("field", "deriv_east", "deriv_north", "deriv_up"),
        coordinate_names=("easting", "northing", "height"),
    ):
        """
        Perform Euler Inversion on data on a regular grid
        """
        shape = grid[data_names[0]].shape
        table = vd.grid_to_table(grid)
        coordinates = [table[n] for n in coordinate_names]
        data = [table[n] for n in data_names]
        self.fit(coordinates, data, weights)
        self.predicted_grid_ = xr.full_like(
            grid,
            {
                data_names[0]: self.predicted_field_.reshape(shape),
                data_names[1]: self.predicted_deriv_east_.reshape(shape),
                data_names[2]: self.predicted_deriv_north_.reshape(shape),
                data_names[3]: self.predicted_deriv_up_.reshape(shape),
            },
        )
        return self

    def fit(self, coordinates, data, weights=(1, 0.1, 0.1, 0.05)):
        """
        Perform Euler Inversion on a single window spanning the entire data
        """
        n_data = data[0].size
        # The data are organized into a single vector because of the maths
        data_observed = np.concatenate(data)
        data_predicted = self._initial_data(coordinates, data)
        parameters, cofactor = self._initial_parameters(coordinates, data)
        # Data weights
        data_weights = np.empty_like(data_predicted)
        data_weights[:n_data] = weights[0]
        data_weights[n_data : 2 * n_data] = weights[1]
        data_weights[2 * n_data : 3 * n_data] = weights[2]
        data_weights[3 * n_data : 4 * n_data] = weights[3]
        Wd_inv = sp.sparse.diags(1 / data_weights, format="csc")
        euler = self._eulers_equation(coordinates, data_predicted, parameters)
        residuals = data_observed - data_predicted
        # Keep track of the way these things vary with iteration
        self.euler_misfit_ = [np.linalg.norm(euler)]
        self.data_misfit_ = [np.linalg.norm(residuals * data_weights)]
        self.merit_ = [
            self.data_misfit_[-1] + self.euler_misfit_balance * self.euler_misfit_[-1]
        ]
        self.location_per_iteration_ = [parameters[:3].copy()]
        self.base_level_per_iteration_ = [parameters[3]]
        for iteration in range(self.max_iterations):
            parameter_step, data_step, cofactor_step = self._newton_step(
                coordinates,
                data_observed,
                data_predicted,
                parameters,
                euler,
                Wd_inv,
            )
            parameters += parameter_step
            data_predicted += data_step
            cofactor += cofactor_step
            euler = self._eulers_equation(coordinates, data_predicted, parameters)
            residuals = data_observed - data_predicted
            # Check stopping criteria
            new_euler_misfit = np.linalg.norm(euler)
            new_data_misfit = np.linalg.norm(residuals * data_weights)
            new_merit = new_data_misfit + self.euler_misfit_balance * new_euler_misfit
            if new_merit > self.merit_[-1]:
                data_predicted -= data_step
                parameters -= parameter_step
                self.stopping_reason_ = "Merit function increased"
                break
            # Update tracked variables
            self.euler_misfit_.append(new_euler_misfit)
            self.data_misfit_.append(new_data_misfit)
            self.merit_.append(new_merit)
            self.location_per_iteration_.append(parameters[:3].copy())
            self.base_level_per_iteration_.append(parameters[3])
            # Check convergence
            merit_change = abs((self.merit_[-2] - self.merit_[-1]) / self.merit_[-2])
            if merit_change < self.tol:
                self.stopping_reason_ = f"Merit change ({merit_change:.2e}) below tolerance ({self.tol:.2e})"
                break
        # Save output attributes
        self.iterations_ = iteration + 1
        self.predicted_field_ = data_predicted[:n_data]
        self.predicted_deriv_east_ = data_predicted[n_data : 2 * n_data]
        self.predicted_deriv_north_ = data_predicted[2 * n_data : 3 * n_data]
        self.predicted_deriv_up_ = data_predicted[3 * n_data :]
        self.location_ = parameters[:3]
        self.base_level_ = parameters[3]
        chi_squared = np.sum(residuals**2) / (residuals.size - parameters.size)
        self.covariance_ = chi_squared * cofactor
        return self

    def _newton_step(
        self, coordinates, data_observed, data_predicted, parameters, euler, Wd_inv
    ):
        """
        Calculate the step in parameters and data in the Gauss-Newton iteration
        """
        # Weights don't seem to make a difference
        east, north, up = coordinates
        field, deriv_east, deriv_north, deriv_up = np.split(data_predicted, 4)
        xo, yo, zo, base_level = parameters
        A = self._parameter_jacobian(deriv_east, deriv_north, deriv_up)
        B = self._data_jacobian(coordinates, parameters[:3])
        residuals = data_observed - data_predicted
        Q = B @ Wd_inv @ B.T
        Q_inv = sp.sparse.linalg.inv(Q)
        ATQ = A.T @ Q_inv
        BTQ = Wd_inv @ B.T @ Q_inv
        Br = B @ residuals
        cofactor_step = sp.linalg.inv(ATQ @ A)
        parameter_step = -cofactor_step @ ATQ @ (euler + Br)
        data_step = residuals - BTQ @ Br - BTQ @ (euler + A @ parameter_step)
        return parameter_step, data_step, cofactor_step

    def _initial_data(self, coordinates, data):
        # Initial estimate for the predicted data is close to the observed
        # data. It can't be exactly the observed data, zero or any other single
        # value because those lead to singular matrices
        data_predicted = 0.9 * np.concatenate(data)
        return data_predicted

    def _initial_parameters(self, coordinates, data):
        # Make an initial estimate for the parameters using the Euler
        # Deconvolution results
        euler_deconv = EulerDeconvolution(structural_index=self.structural_index)
        euler_deconv.fit(coordinates, data)
        parameters = np.empty(4)
        parameters[:3] = euler_deconv.location_
        parameters[3] = euler_deconv.base_level_
        return parameters, euler_deconv.cofactor_

    def _parameter_jacobian(
        self,
        deriv_east,
        deriv_north,
        deriv_up,
    ):
        """
        Calculate the model parameter Jacobian for Euler Inversion
        """
        jacobian = np.empty((deriv_east.size, 4))
        jacobian[:, 0] = -deriv_east
        jacobian[:, 1] = -deriv_north
        jacobian[:, 2] = -deriv_up
        jacobian[:, 3] = -self.structural_index
        return jacobian

    def _data_jacobian(self, coordinates, source_location):
        """
        Calculate the data Jacobian for Euler Inversion
        """
        east, north, up = coordinates
        east_s, north_s, up_s = source_location
        nequations = east.size
        jacobian = sp.sparse.hstack(
            [
                sp.sparse.diags(np.full(nequations, self.structural_index)),
                sp.sparse.diags(east - east_s),
                sp.sparse.diags(north - north_s),
                sp.sparse.diags(up - up_s),
            ],
            format="csc",
        )
        return jacobian

    def _eulers_equation(self, coordinates, data, parameters):
        """
        Evaluate Euler's homogeneity equation
        """
        east, north, up = coordinates
        field, deriv_east, deriv_north, deriv_up = np.split(data, 4)
        east_s, north_s, up_s = parameters[:3]
        base_level = parameters[-1]
        euler = (
            (east - east_s) * deriv_east
            + (north - north_s) * deriv_north
            + (up - up_s) * deriv_up
            + self.structural_index * (field - base_level)
        )
        return euler


class EulerInversionWindowed:

    def __init__(
        self,
        window_size,
        window_step,
        structural_index=None,
        max_variance=0.1,
        max_iterations=20,
        tol=0.1,
        euler_misfit_balance=0.1,
    ):
        self.structural_index = structural_index
        self.window_size = window_size
        self.window_step = window_step
        self.max_iterations = max_iterations
        self.tol = tol
        self.euler_misfit_balance = euler_misfit_balance
        self.max_variance = max_variance

    def fit(self, coordinates, data, weights=(1, 0.1, 0.1, 0.05)):
        _, windows = vd.rolling_window(
            coordinates,
            size=self.window_size,
            spacing=self.window_step,
            adjust="region",
        )
        if self.structural_index is None:
            structural_indices = [1, 2, 3]
        else:
            structural_indices = [self.structural_index]
        pool = concurrent.futures.ProcessPoolExecutor()
        futures = []
        for window in windows.ravel():
            window_coordinates = [c[window[0]] for c in coordinates]
            window_data = [d[window[0]] for d in data]
            future = pool.submit(
                fit_window,
                window_coordinates,
                window_data,
                weights,
                self.max_variance,
                structural_indices,
                dict(
                    tol=self.tol,
                    max_iterations=self.max_iterations,
                    euler_misfit_balance=self.euler_misfit_balance,
                ),
            )
            futures.append(future)
        self.solutions_ = []
        for future in concurrent.futures.as_completed(futures):
            model = future.result()
            if model is not None:
                self.solutions_.append(model)
        self.locations_ = np.transpose([ei.location_ for ei in self.solutions_])
        self.base_levels_ = np.array([ei.base_level_ for ei in self.solutions_])
        self.structural_indices_ = np.array(
            [ei.structural_index for ei in self.solutions_]
        )
        return self

    def fit_grid(
        self,
        grid,
        weights=(1, 0.1, 0.1, 0.05),
        data_names=("field", "deriv_east", "deriv_north", "deriv_up"),
        coordinate_names=("easting", "northing", "height"),
    ):
        """
        Perform Euler Deconvolution on data on a regular grid
        """
        shape = grid[data_names[0]].shape
        table = vd.grid_to_table(grid)
        coordinates = [np.asarray(table[n]) for n in coordinate_names]
        data = [np.asarray(table[n]) for n in data_names]
        self.fit(coordinates, data, weights)
        return self


def fit_window(
    window_coordinates, window_data, weights, max_variance, structural_indices, kwargs
):
    solution = None
    window_region = vd.get_region(window_coordinates)
    candidates = []
    for si in structural_indices:
        ei = EulerInversion(si, **kwargs)
        ei.fit(window_coordinates, window_data, weights)
        inside_window = vd.inside(ei.location_, window_region)
        if not inside_window:
            break
        candidates.append(ei)
    else:
        ei = candidates[np.argmin([ei.data_misfit_[-1] for ei in candidates])]
        small_variance = np.sqrt(ei.covariance_[2, 2]) < max_variance * np.abs(
            ei.location_[2]
        )
        inside_window = vd.inside(ei.location_, window_region)
        if small_variance and inside_window:
            solution = ei
    return solution
