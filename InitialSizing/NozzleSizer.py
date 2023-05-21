# Authors: ravi4ram & Django van der Plas
# This file contains classes that host different nozzle types
import numpy as np
import matplotlib.pyplot as plt
import BellNozzleParameters


# Nozzle types: https://www.researchgate.net/figure/Parabolic-nozzles-e-a-Rao-nozzle-b-Modified-Rao-nozzle-c-Dual-Bell-nozzle-d_fig4_331321763


class ConicalNozzle:
    """
    Class that hosts a conical nozzle
    """

    def __init__(self, exit_diameter: float, throat_diameter: float, conical_half_angle: float = 15,
                 nozzle_arc_scalar: float = 1.5) -> object:
        """
        Instantiate a conical nozzle object
        :param exit_diameter: exit diameter of the nozzle [m]
        :param throat_diameter: throat diameter of the nozzle [m]
        :param conical_half_angle: (half) angle of the conical nozzle [deg]
        :param nozzle_arc_scalar: Parameter that scales the circular arc radius (R in Huzel and Huang)
        """
        self.exit_radius = exit_diameter / 2
        self.throat_radius = throat_diameter / 2
        self.conical_half_angle = conical_half_angle
        self.nozzle_arc_scalar = nozzle_arc_scalar
        self.nozzle_arc_radius = None
        self.contraction_ratio = None
        self.nozzle_length = None
        self.nozzle_correction_factor = None

        self.size_nozzle()

    def size_nozzle(self) -> object:
        """
        Size the nozzle based on the parameters in the class
        :rtype: object
        """
        self.nozzle_arc_radius = self.nozzle_arc_scalar * self.throat_radius
        self.contraction_ratio = self.exit_radius ** 2 / self.throat_radius ** 2
        self.nozzle_length = (self.throat_radius * (np.sqrt(self.contraction_ratio) - 1) +
                              self.nozzle_arc_radius * (1 / np.cos(np.deg2rad(self.conical_half_angle)) - 1)) / \
                             np.tan(np.deg2rad(self.conical_half_angle))
        self.nozzle_correction_factor = 1 / 2 * (1 + np.cos(np.deg2rad(self.conical_half_angle)))


class BellNozzle:
    """
    Thrust optimised parabolic bell nozzle
    based on http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
    and code from https://github.com/ravi4ram/Bell-Nozzle/blob/master/bell_nozzle.py#L137.
    """

    def __init__(self, expansion_ratio: float, throat_radius: float, length_percentage: float = 0.8,
                 entrance_angle: float = 135) -> object:
        """
        Initialises the bell nozzle
        :rtype: object
        :param expansion_ratio: Area expansion ratio
        :param throat_radius: Radius of the throat
        :param length_percentage: Length of the nozzle as a percentage (fraction) of the length of
        the parabolic section (standard 0.8)
        At 85% a nozzle efficiency of 99% is reached, getting to 100% yields only 0.2% more efficiency
        Below 70% the nozzle efficiency suffers
        Most common value is 80%
        :param entrance_angle: Entrance angle of the nozzle (typical is 135 degrees)
        """
        self.ny_bell = None
        self.y_bell = None
        self.x_bell = None
        self.nye2 = None
        self.ye2 = None
        self.xe2 = None
        self.ye = None
        self.nye = None
        self.xe = None
        self.nozzle_length = None
        self.theta_e = None
        self.theta_n = None
        self.expansion_ratio = expansion_ratio
        self.throat_radius = throat_radius
        self.entrance_angle = entrance_angle
        self.length_percentage = length_percentage

        if length_percentage > 1 or length_percentage < 0.6:
            raise ValueError("Length percentage should be between 0.6 and 1, standard 0.8")

        if entrance_angle > 135 or entrance_angle < 90:
            raise ValueError("Entrance angle should be between 90 and 135 degrees, standard 135")

        self.size_nozzle()

    def size_nozzle(self, interval: int = 100) -> object:
        """
        Size the nozzle based on the parameters in the class. The bell nozzle is sized by creating
        a quadratic Bézier curve.
        :rtype: object
        :param interval: amount of points to use for the Bézier curve
        """
        # Bell nozzle parameters
        # θn: bell nozzle start angle, θe: bell nozzle exit angle
        # These parameters are 2D interpolated from empirical data made by Rao
        self.theta_n, self.theta_e = BellNozzleParameters.getBellParameters(self.expansion_ratio,
                                                                            self.length_percentage)
        self.nozzle_length = self.length_percentage * ((np.sqrt(self.expansion_ratio) - 1) *
                                                       self.throat_radius) / np.tan(np.deg2rad(15))
        # Throat entrant section
        # Constructed using a circle
        entrant_angles = np.linspace(np.deg2rad(-self.entrance_angle), np.deg2rad(-90), interval)  # -135 to -90 degrees
        self.xe = 1.5 * self.throat_radius * np.cos(entrant_angles)
        self.ye = 1.5 * self.throat_radius * np.sin(entrant_angles) + 2.5 * self.throat_radius

        # Add negative version
        self.nye = -self.ye

        # Throat exit section (from throat into nozzle)
        # Constructed using a circle
        exit_angles = np.linspace(np.deg2rad(-90), np.deg2rad(self.theta_n - 90), interval)  # -90 to theta_n-90
        self.xe2 = 0.382 * self.throat_radius * np.cos(exit_angles)
        self.xe2 = self.xe2 - self.xe2[0]
        self.ye2 = 0.382 * self.throat_radius * np.sin(exit_angles) + 1.382 * self.throat_radius

        # Add negative version
        self.nye2 = -self.ye2

        # Bell section
        # Drawn using a Quadratic Bézier curve
        # Start point of quadratic Bézier curve
        Nx = 0.382 * self.throat_radius * np.cos(np.deg2rad(self.theta_n - 90))
        Ny = 0.382 * self.throat_radius * np.sin(np.deg2rad(self.theta_n - 90)) + 1.382 * self.throat_radius

        # Exit point of quadratic bezier curve
        Ex = self.nozzle_length
        Ey = np.sqrt(self.expansion_ratio) * self.throat_radius

        # Gradient of the entrance and exit of the nozzle
        # gradient m1,m2 - [Eqn. 8]
        m1 = np.tan(np.deg2rad(self.theta_n))
        m2 = np.tan(np.deg2rad(self.theta_e))

        # Lines used to determine Q
        C1 = Ny - m1 * Nx
        C2 = Ey - m2 * Ex

        # Bezier point
        # Found as intercept between C1 and C2
        Qx = (C2 - C1) / (m1 - m2)
        Qy = (m1 * C2 - m2 * C1) / (m1 - m2)

        # Quadratic Bézier curve
        # The bell is a quadratic Bézier curve, which has equations:
        # x(t) = (1 − t)^2 * Nx + 2(1 − t)t * Qx + t^2 * Ex, 0≤t≤1
        # y(t) = (1 − t)^2 * Ny + 2(1 − t)t * Qy + t^2 * Ey, 0≤t≤1
        coord_list = np.linspace(0, 1, interval)

        self.x_bell = ((1 - coord_list) ** 2) * Nx + 2 * (1 - coord_list) * coord_list * Qx + (coord_list ** 2) * Ex
        self.y_bell = ((1 - coord_list) ** 2) * Ny + 2 * (1 - coord_list) * coord_list * Qy + (coord_list ** 2) * Ey

        # Add a negative version of the bell
        self.ny_bell = -self.y_bell

    def plotNozzle(self, ax):
        ax.set_aspect('equal')

        # throat entrance
        ax.plot(self.xe, self.ye, linewidth=2.5, color='green', label='Throat Entrance')
        ax.plot(self.xe, self.nye, linewidth=2.5, color='green')

        # throat inlet point
        ax.plot(self.xe[0], 0, '+')

        # throat exit
        ax.plot(self.xe2, self.ye2, linewidth=2.5, color='red', label='Throat Exit')
        ax.plot(self.xe2, self.nye2, linewidth=2.5, color='red')

        # bell
        ax.plot(self.x_bell, self.y_bell, linewidth=2.5, color='blue', label='Bell')
        ax.plot(self.x_bell, self.ny_bell, linewidth=2.5, color='blue')


if __name__ == '__main__':
    nozzle = BellNozzle(4, 0.2, length_percentage=0.8, entrance_angle=135)

    fig = plt.figure(figsize=(12, 9))
    ax1 = fig.add_subplot()
    nozzle.plotNozzle(ax1)
    plt.show()
