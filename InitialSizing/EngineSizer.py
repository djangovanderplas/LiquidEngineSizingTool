# Author: Django van der Plas
# This file contains a class that sizes an engine based on the input parameters
import sys
import yaml
import numpy as np

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
        self.pressure_ratio = None
        self.cea = None
        self.mw = None
        self.gamma = None
        self.area_ratio = None
        self.combustion_temperature = None
        self.throat_temperature = None
        self.exit_velocity = None
        self.mass_flow_rate = None
        self.isp = None
        self.throat_area = None
        self.exit_area = None
        self.chamber_volume = None
        self.chamber_length = None
        self.chamber_diameter = None
        self.chamber_area = None
        self.exit_diameter = None
        self.throat_diameter = None
        self.nozzle = None

        # Size the engine
        self.size_engine()

    def size_engine(self):
        """
        Size the engine with the given input parameters. This function will calculate the following parameters:
        - Pressure Ratio
        - Molecular Weight
        - Gamma
        - Area Ratio (Expansion Ratio)
        - Combustion Temperature
        - Throat Temperature
        - Exit Velocity
        - Mass Flow Rate
        - Isp
        - Throat Area
        - Exit Area
        The method features a lot of ugly formulae hence it's not too nice to read. The formulae are taken from Sutton
        and https://github.com/trenton-charlson/engine-dev.
        """
        # Calculate Pressure Ratio
        self.pressure_ratio = self.exit_pressure / self.chamber_pressure

        # Create CEA instance for this engine
        self.cea = propellant.genCEAObj(self.fuel, self.fuel_mass_fraction, self.oxidizer, self.oxidizer_mass_fraction)

        # Get molecular weight and gamma, uses Throat, Chamber can also be used
        self.mw, self.gamma = self.cea.get_Throat_MolWt_gamma(Pc=self.chamber_pressure, MR=self.mixture_ratio,
                                                              eps=self.pressure_ratio)

        # Calculate optimum area expansion ratio
        self.area_ratio = 1 / (
                (((self.gamma + 1) / 2) ** (1 / (self.gamma - 1))) * (self.pressure_ratio ** (1 / self.gamma)) *
                np.sqrt(((self.gamma + 1) / (self.gamma - 1)) *
                        (1 - (self.pressure_ratio ** ((self.gamma - 1) / self.gamma)))))

        # Calculate temperatures
        self.combustion_temperature = self.cea.get_Tcomb(Pc=self.chamber_pressure, MR=self.mixture_ratio)
        self.throat_temperature = 2 * self.combustion_temperature / (1 + self.gamma)

        # Calculate exit velocity
        self.exit_velocity = (2 * self.gamma / (self.gamma - 1) * gas_constant / self.mw * self.combustion_temperature *
                              (1 - self.pressure_ratio ** ((self.gamma - 1) / self.gamma))) ** 0.5

        # Calculate mass flow rate and Isp
        self.mass_flow_rate = self.thrust / self.exit_velocity
        self.isp = self.exit_velocity / g0

        # Calculate nozzle throat area and diameter
        self.throat_area = self.mass_flow_rate / self.chamber_pressure * np.sqrt(
            self.combustion_temperature * gas_constant / self.gamma) * \
                           (1 + (self.gamma - 1) / 2) ** ((self.gamma + 1) / (2 * (self.gamma - 1)))

        self.throat_diameter = np.sqrt(4 / np.pi * self.throat_area)

        # Calculate chamber volume
        self.chamber_volume = self.Lstar * self.throat_area

        # Calculate exit area
        self.exit_area = self.throat_area * self.area_ratio
        self.exit_diameter = np.sqrt(4 / np.pi * self.exit_area)

    def initial_size_geometry(self, nozzle_type='bell'):
        # Calculate chamber area and diameter
        self.chamber_area = self.throat_area * (8.0 * self.throat_diameter ** -0.6 + 1.25) # Emperically determined, mildly sus
        # https://static1.squarespace.com/static/5c639f2e11f78440d05623e3/t/5c6b9076e2c48352e13fa8ab/1550553213416/LCCLS+Rocket+Engine+Sizing+compressed.pdf
        self.chamber_diameter = np.sqrt(4 / np.pi * self.chamber_area)

        # Calculate chamber length
        self.chamber_length = self.chamber_volume / self.chamber_area

        # Size nozzle
        match nozzle_type:
            case 'bell':
                pass
            case 'conical':
                self.nozzle = NozzleSizer.ConicalNozzle(self.exit_diameter, self.throat_diameter,
                                                        conical_half_angle=15, nozzle_arc_scalar=1)
            case 'parabolic':
                pass



if __name__ == '__main__':
    test_engine = Engine('../config.yaml')
    print(vars(test_engine))
