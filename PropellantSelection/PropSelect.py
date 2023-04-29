# Author: Django van der Plas
# This file contains the functions to select a propellant and to check the characteristic length of the propellant
import yaml
import os.path
import pandas as pd
import warnings
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import Propellant
matplotlib.use('Qt5Agg')


def check_characteristic_length(path_to_config: str) -> None:
    """
    This function returns the characteristic length of the selected propellant.
    :param path_to_config:
    :return:
    """
    # Load config file
    with open(path_to_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    # Load propellant data
    fuel = config['Fuel']
    oxidizer = config['Oxidizer']

    # Load propellant data
    propellant_data = pd.read_csv(os.path.join(os.path.dirname(__file__), 'CharacteristicLengths.csv'))

    # Check if the selected propellant is in the database
    comp_in_fuel_list = False
    fuel_list = []
    for comp in fuel:
        if comp in propellant_data['Fuel'].values:
            comp_in_fuel_list = True
            fuel_list.append(comp)
    if not comp_in_fuel_list:
        raise ValueError(f'PropellantSelectionWarning: Fuels {fuel} not in database, please provide your own L* value')
    if len(fuel_list) < len(fuel):
        warnings.warn(f'PropellantSelectionWarning: Not all fuels in database, please use L* value with care')

    comp_in_ox_list = False
    ox_list = []
    for comp in oxidizer:
        if comp in propellant_data['Oxidizer'].values:
            comp_in_oxlist = True
            ox_list.append(comp)
    if not comp_in_oxlist:
        raise ValueError(f'PropellantSelectionWarning: Oxidizers {oxidizer} not in database, please provide your own '
                         f'L* value')
    if len(ox_list) < len(oxidizer):
        warnings.warn(f'PropellantSelectionWarning: Not all oxidizers in database, please use L* value with care')

    # If the selected propellant is in the database, return the characteristic length
    if len(fuel_list) > 0 and len(ox_list) > 0:
        for fuel in fuel_list:
            for oxidizer in ox_list:
                l_star_max = propellant_data['L*max [m]'].where(propellant_data['Fuel'] == fuel).where(
                    propellant_data['Oxidizer'] == oxidizer).dropna().values[0]
                l_star_min = propellant_data['L*min [m]'].where(propellant_data['Fuel'] == fuel).where(
                    propellant_data['Oxidizer'] == oxidizer).dropna().values[0]

        print(f'Characteristic length range for {fuel} and {oxidizer}: {l_star_min} - {l_star_max} [m]')


def plot_OF_Temperature(path_to_config: str, start: float = 1, stop: float = 2, mode='Vac') -> None:
    """
    This function plots the Isp and the combustion temperature vs the mixture ratio.
    :param path_to_config:
    :param start:
    :param stop:
    :param mode:
    :return:
    """
    # Load config file
    with open(path_to_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    # Load engine data
    fuel = config['Fuel']
    fuel_mass_fraction = config['FuelMassFraction']
    oxidizer = config['Oxidizer']
    oxidizer_mass_fraction = config['OxidizerMassFraction']
    chamber_pressure = float(config['ChamberPressure'])
    expansion_ratio = float(config['ExpansionRatio'])

    # Generate CEA object
    cea = Propellant.genCEAObj(fuel, fuel_mass_fraction, oxidizer, oxidizer_mass_fraction)

    # Generate data
    mixture_ratios = np.arange(start, stop, 0.01)
    specific_impulses_vacuum = [cea.get_Isp(Pc=chamber_pressure, MR=mr, eps=expansion_ratio) for mr in mixture_ratios]
    specific_impulses_ambient = [cea.estimate_Ambient_Isp(Pc=chamber_pressure, MR=mr, eps=expansion_ratio,
                                                          Pamb=1.01325e5)[0] for mr in mixture_ratios]
    CombustionTemperatures = [cea.get_Tcomb(Pc=chamber_pressure, MR=mr) for mr in mixture_ratios]
    ThroatTemperatures = [cea.get_Temperatures(Pc=chamber_pressure, MR=mr,
                                               eps=expansion_ratio)[1] for mr in mixture_ratios]

    # Find optimal O/F and corresponding values
    maxIspIndex = np.argmax(specific_impulses_vacuum)
    print(f'O/F for Max Isp: {round(mixture_ratios[maxIspIndex], 2)}')
    print(f'Isp (Vac) for Max Isp: {round(specific_impulses_vacuum[maxIspIndex], 2)} [s]')
    print(f'Isp (Amb) for Max Isp: {round(specific_impulses_ambient[maxIspIndex], 2)} [s]')
    print(f'Combustion Temperature for Max Isp: {round(CombustionTemperatures[maxIspIndex], 2)} [K]')
    print(f'Throat Temperature for Max Isp: {round(ThroatTemperatures[maxIspIndex], 2)} [K]')

    # Plot data
    plt.rcParams["figure.autolayout"] = True
    ax1 = plt.subplot()

    if mode == 'Vac':
        ax1.plot(mixture_ratios, specific_impulses_vacuum, color='blue', label='Vacuum Isp')
    elif mode == 'Amb':
        ax1.plot(mixture_ratios, specific_impulses_ambient, color='blue', label='Ambient Isp')
    else:
        raise ValueError('Invalid mode')

    plt.legend()
    ax1.set_ylabel('Specific impulse [s]')
    ax2 = ax1.twinx()
    ax2.plot(mixture_ratios, CombustionTemperatures, color='red', label='Combustion Temperature')
    ax2.plot(mixture_ratios, ThroatTemperatures, color='orange', label='Throat Temperature')
    ax2.set_ylabel('Temperature [K]')

    plt.xlabel('Mixture ratio [-]')
    plt.legend()
    plt.title(f'{fuel} / {oxidizer} at {chamber_pressure/1e5} bar')
    plt.show()


if __name__ == '__main__':
    path_to_config = os.path.join(os.path.dirname(__file__), '../config.yaml')
    check_characteristic_length(path_to_config)
    plot_OF_Temperature(path_to_config)
