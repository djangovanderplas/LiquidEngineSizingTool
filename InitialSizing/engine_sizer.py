# Author: Django van der Plas
import sys
import yaml
import numpy as np

sys.path.append('../')
from GeneralTools import propellant

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

        # Calculate nozzle throat area
        self.throat_area = self.mass_flow_rate / self.chamber_pressure * np.sqrt(
            self.combustion_temperature * gas_constant / self.gamma) * \
                           (1 + (self.gamma - 1) / 2) ** ((self.gamma + 1) / (2 * (self.gamma - 1)))

        # Calculate exit area
        self.exit_area = self.throat_area * self.area_ratio


if __name__ == '__main__':
    test_engine = Engine('../config.yaml')
    print(vars(test_engine))
