# Author: Django van der Plas
# This file contains a class that sizes an engine based on the input parameters
import sys
import yaml
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.path as mpath

sys.path.append('../')
from GeneralTools import propellant
import NozzleSizer

# Define constants
gas_constant = 8314.46261815324
g0 = 9.81


class Engine:
    """
    Class that hosts an engine
    """
    def __init__(self, path_to_config: str):
        """
        Instantiate an engine object
        :param path_to_config: path to the yaml config file
        """
        # Load config file
        with open(path_to_config) as f:
            config = yaml.load(f, Loader=yaml.FullLoader)

        # Instantiate propellant object
        self.thrust = float(config['Thrust'])
        self.chamber_pressure = float(config['ChamberPressure'])
        self.exit_pressure = float(config['DesignExitPressure'])
        self.mixture_ratio = float(config['MixtureRatio'])
        self.Lstar = float(config['L*'])
        self.fuel = config['Fuel']
        self.fuel_mass_fraction = config['FuelMassFraction']
        self.oxidizer = config['Oxidizer']
        self.oxidizer_mass_fraction = config['OxidizerMassFraction']

        # Engine Sizing Parameters
        self.pressure_ratio = None
        self.cea = None
        self.molecular_weight = None
        self.gamma = None
        self.design_isp = None
        self.combustion_temperature = None
        self.vandenkerckhove_constant = None
        self.exit_velocity = None
        self.mass_flow_rate = None
        self.throat_area = None
        self.exit_area = None
        self.throat_diameter = None
        self.expansion_ratio = None
        self.exit_diameter = None

        # Nozzle Sizing Parameters
        self.nozzle = None
        self.nozzle_type = config['NozzleType']
        self.nozzle_length_percentage = config['NozzleLengthPercentage']
        self.nozzle_conical_half_angle = config['NozzleHalfAngle']
        self.nozzle_arc_scalar = config['NozzleArcScalar']
        self.nozzle_entrance_angle = config['NozzleEntranceAngle']

        # Size the engine
        self.size_engine()

    def size_engine(self):
        # Calculate Pressure Ratio
        self.pressure_ratio = self.exit_pressure / self.chamber_pressure

        # Create CEA instance for this engine
        self.cea = propellant.genCEAObj(self.fuel, self.fuel_mass_fraction, self.oxidizer, self.oxidizer_mass_fraction)

        # Use CEA to calculate molecular weight, gamma, combustion temperature and estimate ambient Isp
        self.molecular_weight, self.gamma = self.cea.get_Chamber_MolWt_gamma(Pc=self.chamber_pressure,
                                                                             MR=self.mixture_ratio,
                                                                             eps=self.pressure_ratio)
        self.design_isp = self.cea.estimate_Ambient_Isp(Pc=self.chamber_pressure,
                                                        MR=self.mixture_ratio,
                                                        eps=self.pressure_ratio,
                                                        Pamb=self.exit_pressure)[0]
        self.combustion_temperature = self.cea.get_Tcomb(Pc=self.chamber_pressure,
                                                         MR=self.mixture_ratio)

        # [1] Sizing Throat
        # Calculate Vandenkerckhove constant
        self.vandenkerckhove_constant = np.sqrt(self.gamma * ((1+self.gamma)/2)**((1+self.gamma)/(1-self.gamma)))

        # Calculate Exit Velocity
        self.exit_velocity = self.design_isp * g0

        # Calculate Mass Flow Rate
        self.mass_flow_rate = self.thrust / self.exit_velocity

        # Calculate Throat Area & Throat Diameter
        self.throat_area = ((self.mass_flow_rate *
                             np.sqrt(gas_constant / self.molecular_weight * self.combustion_temperature))
                            / (self.vandenkerckhove_constant * self.chamber_pressure))

        self.throat_diameter = 2 * np.sqrt(self.throat_area / np.pi)

        # [2] Sizing Diverging Nozzle
        # Calculate optimum area expansion ratio
        self.expansion_ratio = 1 / (
                (((self.gamma + 1) / 2) ** (1 / (self.gamma - 1))) * (self.pressure_ratio ** (1 / self.gamma)) *
                np.sqrt(((self.gamma + 1) / (self.gamma - 1)) *
                        (1 - (self.pressure_ratio ** ((self.gamma - 1) / self.gamma)))))

        # Calculate exit area
        self.exit_area = self.throat_area * self.expansion_ratio
        self.exit_diameter = np.sqrt(4 / np.pi * self.exit_area)

        # Size the nozzle itself
        if self.nozzle_type == 'conical':
            self.nozzle = NozzleSizer.ConicalNozzle(expansion_ratio=self.expansion_ratio,
                                                    throat_radius=self.throat_diameter/2,
                                                    length_percentage=self.nozzle_length_percentage,
                                                    conical_half_angle=self.nozzle_conical_half_angle,
                                                    nozzle_arc_scalar=self.nozzle_arc_scalar,
                                                    entrance_angle=self.nozzle_entrance_angle)
        elif self.nozzle_type == 'bell':
            self.nozzle = NozzleSizer.BellNozzle(expansion_ratio=self.expansion_ratio,
                                                 throat_radius=self.throat_diameter/2,
                                                 length_percentage=self.nozzle_length_percentage,
                                                 nozzle_arc_scalar=self.nozzle_arc_scalar,
                                                 entrance_angle=self.nozzle_entrance_angle)
        else:
            raise ValueError('Nozzle type not recognized')

        # [3] Size Convergent Nozzle
        # Calculate Area Contraction Ratio
        self.contraction_ratio = 1.302 * self.throat_diameter ** -0.481 # Emperically determined, mildly sus,
        # I think it comes from Huzle & Huang
        self.chamber_diameter = self.throat_diameter * np.sqrt(self.contraction_ratio)

        # Calculate the radius of the intersection between chamber wall and convergent nozzle, needs to be less than 0.5
        self.contraction_radius = self.chamber_diameter * 0.3

        # [4] Size Chamber
        # Calculate Chamber Volume
        self.chamber_volume = self.Lstar * self.throat_area
        self.chamber_area = np.pi/4 * self.chamber_diameter**2
        self.chamber_length = self.chamber_volume / self.chamber_area

    def plot_geometry(self, ax):
        ax.set_aspect('equal')

        # [1] Plot Divergent Nozzle
        self.nozzle.plot2D(ax)

        # [2] Plot Convergent Nozzle
        # calculate intersect between chamber wall and convergent nozzle
        x_troat_entrant = self.nozzle.x_throat_entrant[-1]
        y_throat_entrant = self.nozzle.y_throat_entrant[-1]

        y_chamber_nozzle_intersect = self.chamber_diameter/2
        x_chamber_nozzle_intersect = (x_troat_entrant -(y_chamber_nozzle_intersect - y_throat_entrant)
                                      / np.arctan(np.deg2rad(self.nozzle_entrance_angle)))

        x_nozzle_div = np.linspace(x_troat_entrant, x_chamber_nozzle_intersect, 100)
        y_nozzle_div = y_throat_entrant - np.arctan(np.deg2rad(self.nozzle_entrance_angle)) * (x_nozzle_div-x_nozzle_div[0])

        ax.plot(x_nozzle_div, y_nozzle_div, 'purple', label='Diverging Nozzle', linewidth=2.5)
        ax.plot(x_nozzle_div, -y_nozzle_div, 'purple', label='Diverging Nozzle', linewidth=2.5)

        # This is WIP
        # # Plot transient
        # x_testing = np.linspace(x_chamber_nozzle_intersect, x_chamber_nozzle_intersect - 0.3*self.chamber_diameter, 100)
        # y_testing = y_chamber_nozzle_intersect - (x_testing-x_testing[0]) * np.tan(np.deg2rad(-self.nozzle_entrance_angle/2))
        # plt.plot(x_testing, y_testing)
        #
        # x_circle = x_chamber_nozzle_intersect - 0.3*self.chamber_diameter/2*np.cos(np.deg2rad(-self.nozzle_entrance_angle/2))
        # y_circle = y_chamber_nozzle_intersect - 0.3*self.chamber_diameter/2*np.sin(np.deg2rad(self.nozzle_entrance_angle/2))
        # pp1 = mpatches.Arc((x_circle, y_circle), 0.3*self.chamber_diameter/2, 0.3*self.chamber_diameter/2)
        # ax.add_patch(pp1)

        # [3] Plot Chamber Wall
        x_chamber = np.linspace(-self.chamber_length, x_chamber_nozzle_intersect, 100)
        y_chamber = self.chamber_diameter / 2 * np.ones(len(x_chamber))

        ax.plot(x_chamber, y_chamber, 'k', label='Chamber Wall', linewidth=2.5)
        ax.plot(x_chamber, -y_chamber, 'k', label='Chamber Wall', linewidth=2.5)



if __name__ == '__main__':
    test_engine = Engine('../config.yaml')

    fig = plt.figure(figsize=(12, 9))
    ax1 = fig.add_subplot(1, 1, 1)

    # test_engine.initial_size_geometry(nozzle_type='conical')
    test_engine.plot_geometry(ax1)

    plt.legend()
    plt.show()
