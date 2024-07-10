"""
Quick implementation of Euler inversion
"""

import numpy as np
import scipy as sp
import verde as vd


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
        jacobian[:, 0] = -deriv_east
        jacobian[:, 1] = -deriv_north
        jacobian[:, 2] = -deriv_up
        jacobian[:, 3] = -self.structural_index
        pseudo_data = (
            east * deriv_east
            + north * deriv_north
            + up * deriv_up
            + self.structural_index * field
        )
        hessian = jacobian.T @ jacobian
        parameters = scipy.linalg.solve(
            hessian, jacobian.T @ pseudo_data, assume_a="pos"
        )
        pseudo_residuals = pseudo_data - jacobian @ parameters
        chi_squared = np.sum(pseudo_residuals**2) / (pseudo_data.size - parameters.size)
        self.covariance_ = chi_squared * sp.linalg.inv(hessian)
        self.source_location_ = parameters[:3]
        self.base_level_ = parameters[-1]
        return self


class EulerInversion:
    def __init__(self, structural_index, max_iterations=50, tol=1e-2):
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
        n_data = data[0].size
        # The data are organized into a single vector because of the maths
        data_observed = np.concatenate(data)
        data_predicted = self._initial_data(coordinates, data)
        parameters = self._initial_parameters(coordinates)
        euler = self._eulers_equation(coordinates, data_predicted, parameters)
        # Keep track of the way these three functions vary with iteration
        self.euler_misfit_ = [np.sum(np.abs(euler))]
        self.data_misfit_ = [np.linalg.norm(data_observed - data_predicted)]
        self.merit_ = [self.data_misfit_[-1] ** 2 + self.euler_misfit_[-1]]
        for _ in range(self.max_iterations):
            parameter_step, data_step = self._newton_step(
                coordinates,
                data_observed,
                data_predicted,
                parameters,
            )
            parameters += parameter_step
            data_predicted += data_step
            euler = self._eulers_equation(coordinates, data_predicted, parameters)
            self.euler_misfit_.append(np.sum(np.abs(euler)))
            self.data_misfit_.append(np.linalg.norm(data_observed - data_predicted))
            self.merit_.append(self.data_misfit_[-1] ** 2 + self.euler_misfit_[-1])
            merit_change = abs((self.merit_[-2] - self.merit_[-1]) / self.merit_[-2])
            if self.merit_[-1] > self.merit_[-2]:
                self.merit_.pop()
                self.data_misfit_.pop()
                self.euler_misfit_.pop()
                self.stopping_reason_ = "Merit increased"
                break
            self.predicted_field_ = data_predicted[:n_data]
            self.predicted_deriv_east_ = data_predicted[n_data:2 * n_data]
            self.predicted_deriv_north_ = data_predicted[2 * n_data:3 * n_data]
            self.predicted_deriv_up_ = data_predicted[3 * n_data:]
            self.source_location_ = parameters[:3]
            self.base_level_ = parameters[3]
            yield self
            if merit_change < tol:
                self.stopping_reason_ = (
                    f"Merit change ({merit_change:.2e}) below tolerance ({self.tol:.2e})"
                )
                break

    def _initial_data(self, coordinates, data):
        # Initial estimate for the predicted data is 90% of the observed data.
        # It can't be exactly the observed data, zero or any other single value
        # because those lead to singular matrices
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

    def _newton_step(self, coordinates, data_observed, data_predicted, parameters):
        """
        Calculate the step in parameters and data in the Gauss-Newton iteration
        """
        # Weights don't seem to make a difference
        east, north, up = coordinates
        field, dx, dy, dz = np.split(data, 4)
        xo, yo, zo, base_level = parameters
        A = jacobian_parameters(dx, dy, dz, structural_index)
        B = jacobian_data(x, y, z, xo, yo, zo, structural_index)
        r = data_observed - data
        f = eulers_equation(x, y, z, data, parameters, structural_index)
        Q = B @ B.T
        Q_inv = sparse.linalg.inv(Q)
        ATQ = A.T @ Q_inv
        BTQ = B.T @ Q_inv
        Br = B @ r
        deltap = np.linalg.solve(ATQ @ A, -ATQ @ (f + Br))
        deltad = r - BTQ @ Br - BTQ @ (f + A @ deltap)
        return deltap, deltad

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
        jacobian = sparse.hstack(
            [
                sparse.diags(np.full(nequations, self.structural_index)),
                sparse.diags(east - east_s),
                sparse.diags(north - north_s),
                sparse.diags(up - up_s),
            ],
            format="csc",
        )
        return jacobian

    def _eulers_equation(self, coordingates, data, parameters):
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
