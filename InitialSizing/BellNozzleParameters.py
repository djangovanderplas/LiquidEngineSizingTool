# Author: Django van der Plas
# This file contains a function to interpolate the empirical bell nozzle parameters by Rao

import pandas as pd
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt


def Interpolate(x, y):
    data = pd.read_csv('BellAngles.csv')

    X_train = data['expansion_ratio']
    Y_train = data['length_percentage']
    Z1_train = data['theta_n']
    Z2_train = data['theta_e']

    # interpolate the points
    Z1_interpolated = interp.griddata((X_train, Y_train), Z1_train, (x, y), method='cubic')
    Z2_interpolated = interp.griddata((X_train, Y_train), Z2_train, (x, y), method='cubic')

    return Z1_interpolated, Z2_interpolated


def getBellParameters(expansion_ratio: float, length_percentage: float) -> (float, float):
    """
    Get the bell nozzle parameters for a given expansion ratio and length percentage
    :param expansion_ratio:
    :param length_percentage:
    :return: theta_n, theta_e
    """
    length_percentage = length_percentage * 100
    expansion_ratio_interp, length_percentage_interp = np.meshgrid(expansion_ratio, length_percentage)
    theta_n_interp, theta_e_interp = Interpolate(expansion_ratio_interp, length_percentage_interp)
    return theta_n_interp[0, 0], theta_e_interp[0, 0]


if __name__ == '__main__':
    # create a grid of points where you want to interpolate
    X_interp, Y_interp = np.meshgrid(np.linspace(4, 100, 1000), np.linspace(60, 100, 1000))
    Z1_interp, Z2_interp = Interpolate(X_interp, Y_interp)

    # plot the interpolated surface
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X_interp, Y_interp, Z1_interp)
    ax.plot_surface(X_interp, Y_interp, Z2_interp)
    plt.show()
