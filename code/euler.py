"""
Quick implementation of Euler inversion
"""

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
        parameters = sp.linalg.solve(hessian, jacobian.T @ pseudo_data, assume_a="pos")
        pseudo_residuals = pseudo_data - jacobian @ parameters
        chi_squared = np.sum(pseudo_residuals**2) / (pseudo_data.size - parameters.size)
        self.covariance_ = chi_squared * sp.linalg.inv(hessian)
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
        parameters = self._initial_parameters(coordinates, data)
        # Data weights
        data_weights = np.empty_like(data_predicted)
        data_weights[:n_data] = weights[0]
        data_weights[n_data : 2 * n_data] = weights[1]
        data_weights[2 * n_data : 3 * n_data] = weights[2]
        data_weights[3 * n_data : 4 * n_data] = weights[3]
        Wd_inv = sp.sparse.diags(1 / data_weights, format="csc")
        euler = self._eulers_equation(coordinates, data_predicted, parameters)
        # Keep track of the way these things vary with iteration
        self.euler_misfit_ = [np.linalg.norm(euler)]
        self.data_misfit_ = [
            np.linalg.norm((data_observed - data_predicted) * data_weights)
        ]
        self.merit_ = [
            self.data_misfit_[-1] + self.euler_misfit_balance * self.euler_misfit_[-1]
        ]
        self.location_per_iteration_ = [parameters[:3].copy()]
        self.base_level_per_iteration_ = [parameters[3]]
        for iteration in range(self.max_iterations):
            parameter_step, data_step = self._newton_step(
                coordinates,
                data_observed,
                data_predicted,
                parameters,
                euler,
                Wd_inv,
            )
            parameters += parameter_step
            data_predicted += data_step
            euler = self._eulers_equation(coordinates, data_predicted, parameters)
            # Check stopping criteria
            new_euler_misfit = np.linalg.norm(euler)
            new_data_misfit = np.linalg.norm(
                (data_observed - data_predicted) * data_weights
            )
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
        parameter_step = sp.linalg.solve(ATQ @ A, -ATQ @ (euler + Br), assume_a="pos")
        data_step = residuals - BTQ @ Br - BTQ @ (euler + A @ parameter_step)
        return parameter_step, data_step

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
        return parameters

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
