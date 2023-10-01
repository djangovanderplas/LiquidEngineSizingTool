# Authors: ravi4ram & Django van der Plas
# This file contains classes that host different nozzle types
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import BellNozzleParameters


# Nozzle types: https://www.researchgate.net/figure/Parabolic-nozzles-e-a-Rao-nozzle-b-Modified-Rao-nozzle-c-Dual-Bell-nozzle-d_fig4_331321763


class ConicalNozzle:
    """
    Class that hosts a conical nozzle
    """

    def __init__(self, expansion_ratio: float, throat_radius: float, length_percentage=1,
                 conical_half_angle: float = 15, nozzle_arc_scalar: float = 1.5, entrance_angle: float = 135) -> object:
        """
        Instantiate a conical nozzle object
        :param expansion_ratio: Area expansion ratio
        :param throat_radius: Radius of the throat
        :param conical_half_angle: (half) angle of the conical nozzle [deg]
        :param nozzle_arc_scalar: Parameter that scales the circular arc radius (R in Huzel and Huang),
        ranges from 0.5 to 1.5
        """
        self.expansion_ratio = expansion_ratio
        self.throat_radius = throat_radius
        self.throat_area = np.pi * self.throat_radius ** 2
        self.exit_area = self.expansion_ratio * self.throat_area
        self.exit_radius = np.sqrt(4 / np.pi * self.exit_area)
        self.conical_half_angle = conical_half_angle
        self.nozzle_arc_scalar = nozzle_arc_scalar
        self.entrance_angle = entrance_angle
        self.length_percentage = length_percentage
        self.nozzle_arc_radius = None
        self.contraction_ratio = None
        self.nozzle_length = None
        self.nozzle_efficiency = None
        self.neg_y_nozzle = None
        self.y_nozzle = None
        self.x_nozzle = None
        self.neg_y_throat_exit = None
        self.y_throat_exit = None
        self.x_throat_exit = None
        self.neg_y_throat_entrant = None
        self.y_throat_entrant = None
        self.x_throat_entrant = None

        self.size_nozzle()

    def size_nozzle(self, interval=100) -> object:
        """
        Size the nozzle based on the parameters in the class
        :rtype: object
        """
        self.nozzle_arc_radius = self.nozzle_arc_scalar * self.throat_radius
        self.contraction_ratio = self.exit_radius ** 2 / self.throat_radius ** 2
        self.nozzle_length = (self.throat_radius * (np.sqrt(self.contraction_ratio) - 1) +
                              self.nozzle_arc_radius * (1 / np.cos(np.deg2rad(self.conical_half_angle)) - 1)) / \
                             np.tan(np.deg2rad(self.conical_half_angle)) * self.length_percentage
        self.nozzle_efficiency = 1 / 2 * (1 + np.cos(np.deg2rad(self.conical_half_angle)))

        # Throat entrant section
        entrant_angles = np.linspace(np.deg2rad(-90), np.deg2rad(-self.entrance_angle), interval)
        self.x_throat_entrant = self.nozzle_arc_scalar * self.throat_radius * np.cos(entrant_angles)
        self.y_throat_entrant = self.nozzle_arc_scalar * self.throat_radius * np.sin(entrant_angles) \
                                + 2.5 * self.throat_radius
        self.neg_y_throat_entrant = -self.y_throat_entrant

        # Throat exit section
        exit_angles = np.linspace(np.deg2rad(-90), np.deg2rad(self.conical_half_angle-90), interval)
        self.x_throat_exit = self.nozzle_arc_scalar * self.throat_radius * np.cos(exit_angles)
        self.x_throat_exit = self.x_throat_exit - self.x_throat_exit[0]
        self.y_throat_exit = self.nozzle_arc_scalar * self.throat_radius * np.sin(exit_angles) \
                             + 2.5 * self.throat_radius
        self.neg_y_throat_exit = -self.y_throat_exit

        # Nozzle section
        slope = np.tan(np.deg2rad(self.conical_half_angle))
        self.x_nozzle = np.linspace(self.x_throat_exit[-1], self.nozzle_length, interval)
        self.y_nozzle = slope * (self.x_nozzle - self.x_throat_exit[-1]) + self.y_throat_exit[-1]
        self.neg_y_nozzle = -self.y_nozzle

    def plot2D(self, ax):
        ax.set_aspect('equal')

        # throat entrance
        ax.plot(self.x_throat_entrant, self.y_throat_entrant, linewidth=2.5, color='green', label='Throat Entrance')
        ax.plot(self.x_throat_entrant, self.neg_y_throat_entrant, linewidth=2.5, color='green')

        # throat inlet point
        ax.plot(self.x_throat_exit[0], 0, '+', label='Throat')

        # throat exit
        ax.plot(self.x_throat_exit, self.y_throat_exit, linewidth=2.5, color='red', label='Throat Exit')
        ax.plot(self.x_throat_exit, self.neg_y_throat_exit, linewidth=2.5, color='red')

        # bell
        ax.plot(self.x_nozzle, self.y_nozzle, linewidth=2.5, color='blue', label='Bell')
        ax.plot(self.x_nozzle, self.neg_y_nozzle, linewidth=2.5, color='blue')

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('2D Nozzle')
        ax.legend()

    def plot3D(self, ax):
        ax.set_aspect('equal')

        x_list = np.concatenate((self.x_throat_entrant[::-1], self.x_throat_exit, self.x_nozzle))
        y_list = np.concatenate((self.y_throat_entrant[::-1], self.y_throat_exit, self.y_nozzle))

        Theta = np.linspace(0, 2.1 * np.pi, 1000)
        X, Theta = np.meshgrid(x_list, Theta)
        Y = y_list * np.cos(Theta)
        Z = y_list * np.sin(Theta)

        ax.plot_surface(X, Y, Z, cmap='viridis')

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title('3D Nozzle')


class BellNozzle:
    """
    Parabolic approximation of the bell nozzle
    based on http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
    and code from https://github.com/ravi4ram/Bell-Nozzle/blob/master/bell_nozzle.py#L137.
    """

    def __init__(self, expansion_ratio: float, throat_radius: float, length_percentage: float = 0.8,
                 nozzle_arc_scalar: float = 1.5, entrance_angle: float = 135) -> object:
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
        :param nozzle_arc_scalar: Scalar of the radius of the entrance. Should be 1.5 in a standard bell nozzle
        :param entrance_angle: Entrance angle of the nozzle (typical is 135 degrees)
        """
        self.expansion_ratio = expansion_ratio
        self.throat_radius = throat_radius
        self.entrance_angle = entrance_angle
        self.length_percentage = length_percentage
        self.nozzle_arc_scalar = nozzle_arc_scalar
        self.neg_y_nozzle = None
        self.y_nozzle = None
        self.x_nozzle = None
        self.neg_y_throat_exit = None
        self.y_throat_exit = None
        self.x_throat_exit = None
        self.y_throat_entrant = None
        self.neg_y_throat_entrant = None
        self.x_throat_entrant = None
        self.nozzle_length = None
        self.theta_e = None
        self.theta_n = None

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

        entrant_angles = np.linspace(np.deg2rad(-90), np.deg2rad(-self.entrance_angle), interval)  # -135 to -90 degrees
        self.x_throat_entrant = self.nozzle_arc_scalar * self.throat_radius * np.cos(-entrant_angles)
        # self.ye = 2.5* self.throat_radius + self.nozzle_arc_scalar * self.throat_radius * np.sin(entrant_angles)
        self.y_throat_entrant = self.nozzle_arc_scalar * self.throat_radius * np.sin(entrant_angles) + \
                                self.throat_radius * (self.nozzle_arc_scalar + 1)
        self.entrant_length = self.x_throat_entrant[0] - self.x_throat_entrant[-1]

        # Add negative version
        self.neg_y_throat_entrant = -self.y_throat_entrant

        # Throat exit section (from throat into nozzle)
        # Constructed using a circle
        exit_angles = np.linspace(np.deg2rad(-90), np.deg2rad(self.theta_n - 90), interval)  # -90 to theta_n-90
        self.x_throat_exit = 0.382 * self.throat_radius * np.cos(exit_angles)
        self.x_throat_exit = self.x_throat_exit - self.x_throat_exit[0]
        self.y_throat_exit = 0.382 * self.throat_radius * np.sin(exit_angles) + 1.382 * self.throat_radius

        # Add negative version
        self.neg_y_throat_exit = -self.y_throat_exit

        # Bell section
        # Drawn using a Quadratic Bézier curve
        # Start point of quadratic Bézier curve, N
        Nx = 0.382 * self.throat_radius * np.cos(np.deg2rad(self.theta_n - 90))
        Ny = 0.382 * self.throat_radius * np.sin(np.deg2rad(self.theta_n - 90)) + 1.382 * self.throat_radius

        # Exit point of quadratic bezier curve, E
        Ex = self.nozzle_length
        Ey = np.sqrt(self.expansion_ratio) * self.throat_radius

        # Gradient of the entrance and exit of the nozzle
        # gradient slope1, slope2
        slope1 = np.tan(np.deg2rad(self.theta_n))
        slope2 = np.tan(np.deg2rad(self.theta_e))

        # Lines used to determine Q
        line1 = Ny - slope1 * Nx
        line2 = Ey - slope2 * Ex

        # Bezier point
        # Found as intercept between line1 and C2
        Qx = (line2 - line1) / (slope1 - slope2)
        Qy = (slope1 * line2 - slope2 * line1) / (slope1 - slope2)

        # Quadratic Bézier curve
        # The bell is a quadratic Bézier curve, which has equations:
        # x(t) = (1 − t)^2 * Nx + 2(1 − t)t * Qx + t^2 * Ex, 0≤t≤1
        # y(t) = (1 − t)^2 * Ny + 2(1 − t)t * Qy + t^2 * Ey, 0≤t≤1
        coord_list = np.linspace(0, 1, interval)

        self.x_nozzle = ((1 - coord_list) ** 2) * Nx + 2 * (1 - coord_list) * coord_list * Qx + (coord_list ** 2) * Ex
        self.y_nozzle = ((1 - coord_list) ** 2) * Ny + 2 * (1 - coord_list) * coord_list * Qy + (coord_list ** 2) * Ey

        # Add a negative version of the bell
        self.neg_y_nozzle = -self.y_nozzle

    def plot2D(self, ax):
        ax.set_aspect('equal')

        # throat entrance
        ax.plot(self.x_throat_entrant, self.y_throat_entrant, linewidth=2.5, color='green', label='Throat Entrance')
        ax.plot(self.x_throat_entrant, self.neg_y_throat_entrant, linewidth=2.5, color='green')

        # throat inlet point
        ax.plot(self.x_throat_exit[0], 0, '+', label='Throat')

        # throat exit
        ax.plot(self.x_throat_exit, self.y_throat_exit, linewidth=2.5, color='red', label='Throat Exit')
        ax.plot(self.x_throat_exit, self.neg_y_throat_exit, linewidth=2.5, color='red')

        # bell
        ax.plot(self.x_nozzle, self.y_nozzle, linewidth=2.5, color='blue', label='Bell')
        ax.plot(self.x_nozzle, self.neg_y_nozzle, linewidth=2.5, color='blue')

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('2D Nozzle')
        ax.legend()



    def plot3D(self, ax):
        ax.set_aspect('equal')

        x_list = np.concatenate((self.x_throat_entrant[::-1], self.x_throat_exit, self.x_nozzle))
        y_list = np.concatenate((self.y_throat_entrant[::-1], self.y_throat_exit, self.y_nozzle))

        Theta = np.linspace(0, 2.1 * np.pi, 1000)
        X, Theta = np.meshgrid(x_list, Theta)
        Y = y_list * np.cos(Theta)
        Z = y_list * np.sin(Theta)

        ax.plot_surface(X, Y, Z, cmap='viridis')

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title('3D Nozzle')



if __name__ == '__main__':
    #nozzle = BellNozzle(8, 0.2, length_percentage=1.0, nozzle_arc_scalar=1.5, entrance_angle=135)
    nozzle = ConicalNozzle(4, 0.1, length_percentage=1, conical_half_angle=15, nozzle_arc_scalar=1.5,entrance_angle=135)

    fig = plt.figure(figsize=(12, 9))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    nozzle.plot2D(ax1)
    nozzle.plot3D(ax2)
    plt.show()
