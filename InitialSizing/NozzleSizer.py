import numpy as np
import BellNozzleParameters


class ConicalNozzle:
    """
    Class that hosts a conical nozzle
    """
    def __init__(self, exit_diameter: float, throat_diameter: float, conical_half_angle: float = 15,
                 nozzle_arc_scalar: float = 1.5) -> object:
        """
        Instantiate a conical nozzle object
        :rtype: object
        :param exit_diameter: exit diameter of the nozzle [m]
        :param throat_diameter: throat diameter of the nozzle [m]
        :param conical_half_angle: (half) angle of the conical nozzle [deg]
        :param nozzle_arc_scalar: Parameter that scales the scales the circular arc radius (R in Huzel and Huang)
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


class BellNozzleV1:
    # Thrust optimised parabolic nozzle
    def __init__(self, expansion_ratio, throat_diameter, length_percentage=0.8, entrance_angle=135):
        self.expansion_ratio = expansion_ratio
        self.throat_radius = throat_diameter / 2
        self.entrance_angle = entrance_angle
        self.length_percentage = length_percentage

    def size_nozzle(self):
        self.theta_n, self.theta_e = BellNozzleParameters.get_bell_nozzle_parameters(self.expansion_ratio, self.length_percentage)
        self.nozzle_length = self.length_percentage * ((np.sqrt(self.expansion_ratio) - 1) * self.throat_radius)/ np.tan(np.deg2rad(15))
        # https://github.com/ravi4ram/Bell-Nozzle/blob/master/bell_nozzle.py#L137
