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
        self.source_location_ = parameters[:3]
        self.base_level_ = parameters[-1]
        return self

    def fit_grid(
        self,
        field,
        deriv_east,
        deriv_north,
        deriv_up,
        coords_names=("easting", "northing", "height"),
    ):
        """
        Perform Euler Deconvolution on data on a regular grid
        """
        grid = xr.Dataset(
            dict(
                field=field,
                deriv_east=deriv_east,
                deriv_north=deriv_north,
                deriv_up=deriv_up,
            )
        )
        table = vd.grid_to_table(grid)
        self.fit(
            coordinates=[table[c] for c in coords_names],
            data=[table.field, table.deriv_east, table.deriv_north, table.deriv_up],
        )
        return self


class EulerInversion:

    def __init__(self, structural_index, max_iterations=20, tol=0.1):
        self.structural_index = structural_index
        self.max_iterations = max_iterations
        self.tol = tol

    def fit(self, coordinates, data):
        """
        Perform Euler Inversion on a single window spanning the entire data
        """
        for _ in self.fit_iterator(coordinates, data):
            continue
        return self

    def fit_iterator(self, coordinates, data):
        """
        Generator of the Euler Inversion iterations
        """
        balance = 0.1
        n_data = data[0].size
        # The data are organized into a single vector because of the maths
        data_observed = np.concatenate(data)
        data_predicted = self._initial_data(coordinates, data)
        parameters = self._initial_parameters(coordinates)
        # Data weights
        weights = np.ones_like(data_predicted)
        weights[:n_data] /= np.linalg.norm(data[0])
        weights[n_data:2 * n_data] /= np.linalg.norm(data[1]) * 10
        weights[2 * n_data:3 * n_data] /= np.linalg.norm(data[2]) * 10
        weights[3 * n_data:4 * n_data] /= np.linalg.norm(data[3]) * 20
        weights /= weights.max()
        Wd_inv = sp.sparse.diags(1 / weights, format="csc")
        euler = self._eulers_equation(coordinates, data_predicted, parameters)
        # Keep track of the way these three functions vary with iteration
        self.euler_misfit_ = [np.linalg.norm(euler)]
        self.data_misfit_ = [np.linalg.norm((data_observed - data_predicted) * weights)]
        self.merit_ = [self.data_misfit_[-1] + balance * self.euler_misfit_[-1]]
        self.predicted_field_ = data_predicted[:n_data]
        self.predicted_deriv_east_ = data_predicted[n_data : 2 * n_data]
        self.predicted_deriv_north_ = data_predicted[2 * n_data : 3 * n_data]
        self.predicted_deriv_up_ = data_predicted[3 * n_data :]
        self.source_location_ = parameters[:3]
        self.base_level_ = parameters[3]
        yield self
        for _ in range(self.max_iterations):
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
            self.euler_misfit_.append(np.linalg.norm(euler))
            self.data_misfit_.append(np.linalg.norm((data_observed - data_predicted) * weights))
            self.merit_.append(self.data_misfit_[-1] + balance * self.euler_misfit_[-1])
            merit_change = abs((self.merit_[-2] - self.merit_[-1]) / self.merit_[-2])
            if self.merit_[-1] > self.merit_[-2]:
                self.merit_.pop()
                self.data_misfit_.pop()
                self.euler_misfit_.pop()
                self.stopping_reason_ = "Merit increased"
                break
            self.predicted_field_ = data_predicted[:n_data]
            self.predicted_deriv_east_ = data_predicted[n_data : 2 * n_data]
            self.predicted_deriv_north_ = data_predicted[2 * n_data : 3 * n_data]
            self.predicted_deriv_up_ = data_predicted[3 * n_data :]
            self.source_location_ = parameters[:3]
            self.base_level_ = parameters[3]
            yield self
            if merit_change < self.tol:
                self.stopping_reason_ = f"Merit change ({merit_change:.2e}) below tolerance ({self.tol:.2e})"
                break

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

    def fit_grid(
        self,
        field,
        deriv_east,
        deriv_north,
        deriv_up,
        coords_names=("easting", "northing", "height"),
    ):
        """
        Perform Euler Deconvolution on data on a regular grid
        """
        for _ in self.fit_grid_iterator(
            field, deriv_east, deriv_north, deriv_up, coords_names
        ):
            continue
        return self

    def fit_grid_iterator(
        self,
        field,
        deriv_east,
        deriv_north,
        deriv_up,
        coords_names=("easting", "northing", "height"),
    ):
        """
        Perform Euler Deconvolution on data on a regular grid
        """
        grid = xr.Dataset(
            dict(
                field=field,
                deriv_east=deriv_east,
                deriv_north=deriv_north,
                deriv_up=deriv_up,
            )
        )
        table = vd.grid_to_table(grid)
        coordinates = [table[c] for c in coords_names]
        data = [table.field, table.deriv_east, table.deriv_north, table.deriv_up]
        for _ in self.fit_iterator(coordinates, data):
            self.predicted_grids_ = xr.full_like(
                grid,
                {
                    "field": self.predicted_field_.reshape(field.shape),
                    "deriv_east": self.predicted_deriv_east_.reshape(field.shape),
                    "deriv_north": self.predicted_deriv_north_.reshape(field.shape),
                    "deriv_up": self.predicted_deriv_up_.reshape(field.shape),
                },
            )
            yield self

    def _initial_data(self, coordinates, data):
        # Initial estimate for the predicted data is close to the observed
        # data. It can't be exactly the observed data, zero or any other single
        # value because those lead to singular matrices
        data_predicted = 0.9 * np.concatenate(data)
        return data_predicted

    def _initial_parameters(self, coordinates):
        # Make an initial estimate for the parameters that places the source in
        # the center of the data, at 10% of the height below the data, and with
        # 0 base level
        region = vd.get_region(coordinates)
        mean_height = np.mean(coordinates[2])
        parameters = np.array(
            [
                0.5 * (region[1] + region[0]),
                0.5 * (region[3] + region[2]),
                mean_height - np.abs(0.1 * mean_height),
                0,
            ]
        )
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
