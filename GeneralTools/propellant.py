from pycea import CEA, blends
import yaml

def genCEAObj(fuel, fuel_mass_fraction, oxidizer, oxidizer_mass_fraction):
    """
    This function generates a CEA object based on the given propellant and mass fractions.
    :param fuel:
    :param fuel_mass_fraction:
    :param oxidizer:
    :param oxidizer_mass_fraction:
    :return:
    """
    # Create new blends
    Fuel = blends.newFuelBlend(fuelL=fuel, fuelPcentL=fuel_mass_fraction)
    Oxidizer = blends.newOxBlend(oxL=oxidizer, oxPcentL=oxidizer_mass_fraction)

    # Create CEA object
    cea = CEA(propName='', fuelName=Fuel, oxName=Oxidizer, units="metric")
    return cea

def genCEAfromConfig(path_to_config):
    with open(path_to_config) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    fuel = config['Fuel']
    fuel_mass_fraction = config['FuelMassFraction']
    oxidizer = config['Oxidizer']
    oxidizer_mass_fraction = config['OxidizerMassFraction']

    return genCEAObj(fuel, fuel_mass_fraction, oxidizer, oxidizer_mass_fraction)

if __name__ == '__main__':
    with open('../config.yaml') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    fuel = config['Fuel']
    fuel_mass_fraction = config['FuelMassFraction']
    oxidizer = config['Oxidizer']
    oxidizer_mass_fraction = config['OxidizerMassFraction']

    cea = genCEAObj(fuel, fuel_mass_fraction, oxidizer, oxidizer_mass_fraction)
    p_test = 305e5
    of_test = 1.02
    print(f'IspVac: {round(cea.get_Isp(Pc=p_test, MR=of_test, eps=10), 2)}/{round(cea.get_Isp(Pc=p_test, MR=1.49, eps=10), 2)}[s]')
    print(f'IspSL: {round(cea.estimate_Ambient_Isp(Pc=p_test, MR=of_test, eps=10, Pamb=1e5)[0], 2)}/{round(cea.estimate_Ambient_Isp(Pc=p_test, MR=1.49, eps=10, Pamb=1e5)[0], 2)}[s], {cea.estimate_Ambient_Isp(Pc=p_test, MR=of_test, eps=10, Pamb=1e5)[1]}')
    print(f'C*: {round(cea.get_Cstar(Pc=p_test, MR=of_test), 2)}/{round(cea.get_Cstar(Pc=p_test, MR=1.49), 2)}[m/s]')
    print(f'Tcomb: {round(cea.get_Tcomb(Pc=p_test, MR=of_test), 2)}/{round(cea.get_Tcomb(Pc=p_test, MR=1.49), 2)}[K]')
    print(f'Ct: {round(cea.get_PambCf(Pc=p_test, MR=of_test, eps=10, Pamb=1e5)[1], 4)}/{round(cea.get_PambCf(Pc=p_test, MR=1.49, eps=10, Pamb=1e5)[1], 4)}')