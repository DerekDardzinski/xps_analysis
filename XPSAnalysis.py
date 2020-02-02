import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import copy as cp


def import_database_to_df(core_level: str) -> pd.DataFrame:
    """
    This function imports the core level cross section and asymmetry parameter values at different
    energy values. The cross section and asymmetry parameters are used to calculate the differential
    cross section, needed for the calculation of the sensitivity factors

    :param core_level: str (Input the desired core level eg. 'Mo3d', 'C1s', 'W4f', ...)
    :return: df: pd.DataFrame (df of cross section and asymmetry parameters at different energies)
    """

    # DataBase name
    # data_base_path = 'DataBase'
    data_base_path = 'https://raw.githubusercontent.com/DerekDardzinski/xps_analysis/master/DataBase/'

    # File type
    file_type = 'txt'

    # Create file path
    file_path = f'{data_base_path}/{core_level}.{file_type}'

    # Read the data into a dataframe
    df = pd.read_csv(file_path, delimiter="\t")

    # Return the dataframe
    return df


def calculate_cross_section(x_ray_source_energy: float, df: pd.DataFrame) -> float:

    x_name = 'Energy [eV]'
    y_name = 'Cross Section'

    # Import the x and y data
    x = df[x_name].__array__(dtype=float)
    y = df[y_name].__array__(dtype=float)

    # If statement to find out whether the x_val is already in the data set
    if np.sum(x == x_ray_source_energy) != 0:

        # y_val is the value at the index where the x_val is
        y_val = y[np.argwhere(x == x_ray_source_energy)].flatten()[0]

    # Else is the x_val is not already in the data set
    else:

        # Insert the x_val into the array and sort it in ascending order
        sorted_array = np.sort(np.insert(x, 0, x_ray_source_energy))

        # Find the index of the value in the sorted array
        index = np.argwhere(sorted_array == x_ray_source_energy).flatten()[0]

        # The lower x and y values for linear interpolation is at the index - 1
        x1 = x[index - 1]
        y1 = y[index - 1]

        # The upper x and y values for linear interpolation is at the index
        x2 = x[index]
        y2 = y[index]

        # Use linear interpolation to solve for the y_val
        y_val = y1 + (x_ray_source_energy - x1) * ((y2 - y1) / (x2 - x1))

    # Return the y_val
    return y_val


def calculate_asymmetry_parameter(x_ray_source_energy: float, df: pd.DataFrame) -> float:

    x_name = 'Energy [eV]'
    y_name = 'Asymmetry Parameter'

    # Import the x and y data
    x = df[x_name].__array__(dtype=float)
    y = df[y_name].__array__(dtype=float)

    # If statement to find out whether the x_val is already in the data set
    if np.sum(x == x_ray_source_energy) != 0:

        # y_val is the value at the index where the x_val is
        y_val = y[np.argwhere(x == x_ray_source_energy)].flatten()[0]

    # Else is the x_val is not already in the data set
    else:

        # Insert the x_val into the array and sort it in ascending order
        sorted_array = np.sort(np.insert(x, 0, x_ray_source_energy))

        # Find the index of the value in the sorted array
        index = np.argwhere(sorted_array == x_ray_source_energy).flatten()[0]

        # The lower x and y values for linear interpolation is at the index - 1
        x1 = x[index - 1]
        y1 = y[index - 1]

        # The upper x and y values for linear interpolation is at the index
        x2 = x[index]
        y2 = y[index]

        # Use linear interpolation to solve for the y_val
        y_val = y1 + (x_ray_source_energy - x1) * ((y2 - y1) / (x2 - x1))

    # Return the y_val
    return y_val


def calculate_differential_cross_section(cross_section: float, asymmetry_parameter: float, analyzer_angle: float) -> float:
    """
    This function calculates the differential photoionization cross section of a given core level.
    The values for the photoionization cross section and asymmetry parameter are extracted from
    the database. The analyzer angle must be specified by the user.

    :param cross_section: float (Photoionization Cross Section)
    :param asymmetry_parameter: float (Asymmetry Parameter)
    :param analyzer_angle: float (Analyzer Angle)
    :return: differential_cross_section: float (Differential Photoionization Cross Section)
    """

    # Split the differential cross section equation up into four separate parts
    A = cross_section / (4 * np.pi)
    B = 1
    C = asymmetry_parameter / 2
    D = ((3 / 2) * (np.sin(analyzer_angle * (np.pi/180)) ** 2)) - 1

    # Calculate the differential photoionization cross section
    differential_cross_section = A * (B + (C * D))

    return differential_cross_section


def import_periodic_table_information(file_path: str) -> dict:
    """
    This function imports the atomic mass, atomic number, and element symbol from a .txt file

    :param file_path: str (File path to the .txt file containing the periodic table information)
    :return:
    """
    period_table_info = pd.read_csv(file_path, delimiter='\t')

    atomic_number = period_table_info['AtomicNumber'].__array__(dtype=int)
    element_symbol = period_table_info['Symbol'].__array__(dtype=str)
    atomic_mass = period_table_info['AtomicMass'].__array__(dtype=float)

    periodic_table_dict = {'Atomic Number': dict(zip(element_symbol, atomic_number)),
                           'Atomic Mass': dict(zip(element_symbol, atomic_mass))}

    return periodic_table_dict


def split_chemical_formula(chemical_formula: str) -> dict:
    """
    This function takes in a chemical formula breaks it to its individual constituents

    :param chemical_formula: str (Chemical Formula of species being analyzed)
    :return:
    """

    # This splits the formula at each Captital Letter eg. 'MoS2' -> ['Mo', 'S2']
    split_formula_by_caps = re.findall('[A-Z][^A-Z]*', chemical_formula)

    # Empty list for later extending
    chemical_formula_components_dict = {}

    # For loop to loop through each term in the split_formula_by_caps list
    for split_term in split_formula_by_caps:

        # Split each sub-component by number eg. 'S2' -> ['S', '2']
        split_component_by_number = re.findall(r'[A-Za-z]+|\d+', split_term)

        if len(split_component_by_number) == 1:
            split_component_by_number.extend('1')

        chemical_formula_components_dict[split_component_by_number[0]] = float(
            split_component_by_number[1])

        # Extend the chemical_formula_components_dict list
        # chemical_formula_components_dict.extend(split_component_by_number)

    # Return the dictionary
    return chemical_formula_components_dict


def calculate_average_atomic_mass(chemical_formula_components: dict, periodic_table_information_dict: dict) -> float:
    """

    :param chemical_formula_components:
    :param periodic_table_information_dict:
    :return:
    """
    elements_in_compound = list(chemical_formula_components.keys())

    number_of_each_element = np.array(list(chemical_formula_components.values()))

    total_number_of_elements = sum(number_of_each_element)

    atomic_masses = [periodic_table_information_dict['Atomic Mass'][element]
                     for element in elements_in_compound]

    average_atomic_mass = sum(atomic_masses * number_of_each_element) / total_number_of_elements

    return average_atomic_mass


def calculate_average_atomic_number(chemical_formula_components: dict, periodic_table_information_dict: dict) -> float:
    """

    :param chemical_formula_components:
    :param periodic_table_information_dict:
    :return:
    """
    elements_in_compound = list(chemical_formula_components.keys())

    number_of_each_element = np.array(list(chemical_formula_components.values()))

    total_number_of_elements = sum(number_of_each_element)

    atomic_numbers = [periodic_table_information_dict['Atomic Number'][element]
                      for element in elements_in_compound]

    average_atomic_number = sum(atomic_numbers * number_of_each_element) / total_number_of_elements

    return average_atomic_number


# chem_form = split_chemical_formula('MoS2')
# periodic_table = import_periodic_table_information('PeriodicTableData.txt')
# avg_mass = calculate_average_atomic_mass(chem_form, periodic_table)
#
# print(avg_mass)


def calculate_lattice_parameter(average_atomic_mass: float, density: float) -> float:
    """
    This function calculates the average lattice parameter of a material. Given the name of the
    material, the function will determine the molecular weight of the compound. Also from the name,
    the density will be looked up from an list of densities that is input by the user.

    :param average_atomic_mass: float (Average atomic mass of the chemical species)
    :param density: float (Density of material) [kg / m^3]
    :return: lattice_parameter_value: float (Lattice Parameter) [m]
    """

    avagadro_number = 6.0221409E23
    lattice_parameter_value = 1E8 * ((average_atomic_mass / (density * avagadro_number)) ** (1/3))

    return lattice_parameter_value


def calculate_attenuation_length(lattice_parameter: float, avg_atomic_number: float, kinetic_energy: float) -> float:
    """
    This function calculates the X-ray attenuation length of a material given the lattice parameter,
    average molecular weight, and kinetic energy. All values will be calculated without user input.

    :param lattice_parameter: float (lattice parameter of the material)
    :param avg_atomic_number: float (Average atomic number of the matrix)
    :param kinetic_energy: float (Kinetic energy)
    :return: attenuation_length_value: float (X-ray attenuation length)
    """
    A = 0.316 * lattice_parameter ** (3/2)
    B = kinetic_energy
    C = avg_atomic_number ** 0.45
    D = np.log(kinetic_energy / 27) + 3

    attenuation_length = A * ((B / (C * D)) + 4)

    return attenuation_length


def calculate_sensitivity_factor(differential_cross_section: float, attenuation_length: float) -> float:
    """
    This function calculates the sensitivity factor of a core level

    :param differential_cross_section: float (Differential cross section value)
    :param attenuation_length: float (X-ray attenuation length)
    :return: sensitivity_factor: float (Sensitivity factor)
    """

    sensitivity_factor = differential_cross_section * attenuation_length

    return sensitivity_factor


def calculate_kinetic_energy_dict(position_dict: dict, x_ray_source_energy: float,
                                  detector_work_function: float) -> dict:

    kinetic_energy_dict = cp.deepcopy(position_dict)

    for key in position_dict.keys():
        kinetic_energy_dict[key] = x_ray_source_energy - position_dict[key] - detector_work_function

    return kinetic_energy_dict

# def create_sensitivity_factor_dict(kinetic_energy_dict: dict, core_level: str, analyzer_angle: float,
#                                    density_dict: dict, period_table_information: dict,
#                                    x_ray_source_energy: float) -> dict:
#     """
#     This function creates a dictionary of sensitivity factors for each species in present in a core level
#
#     :param kinetic_energy_dict: dict (Dictionary that contains the kinetic energies of each peak)
#     :param core_level: str (Core level that is being analyzed. Necessary for calculating cross section and asymmetry
#     parameters)
#     :param analyzer_angle: float (Angle between the source and the analyzer. Used in the calculation of the
#     differential cross section)
#     :param density_dict: dict (Dictionary containing the densities of chemical species present on the sample. Used in
#     the calculation of the lattice parameter which is need to calculate the x-ray attenuation length)
#     :param period_table_information: dict (Dictionary containing atomic numbers, atomic masses, and element symbols
#     for each element)
#     :param x_ray_source_energy: float (Energy of your x-ray source)
#     :return: sensitivity_factor_dict: dict (Dictionary containing sensitivity factors for each species on the surface
#     of the sample for a given core level)
#     """
#
#     df = import_database_to_df(core_level)
#
#     sensitivity_factor_dict = cp.deepcopy(kinetic_energy_dict)
#
#     maximum_intensity_chemical_formula = max(sensitivity_factor_dict, key=sensitivity_factor_dict.get())
#
#     for chemical_formula in kinetic_energy_dict.keys():
#
#         density = density_dict[chemical_formula]
#         split_formula = split_chemical_formula(chemical_formula=chemical_formula)
#
#         average_atomic_mass = calculate_average_atomic_mass(
#             chemical_formula_components=split_formula,
#             periodic_table_information_dict=period_table_information
#         )
#
#         average_atomic_number = calculate_average_atomic_number(
#             chemical_formula_components=split_formula,
#             periodic_table_information_dict=period_table_information
#         )
#
#         cross_section = calculate_cross_section(
#             x_ray_source_energy=x_ray_source_energy,
#             df=df
#         )
#
#         asymmetry_parameter = calculate_asymmetry_parameter(
#             x_ray_source_energy=x_ray_source_energy,
#             df=df
#         )
#
#         differential_cross_section = calculate_differential_cross_section(
#             cross_section=cross_section,
#             asymmetry_parameter=asymmetry_parameter,
#             analyzer_angle=analyzer_angle
#         )
#
#         lattice_parameter = calculate_lattice_parameter(
#             average_atomic_mass=average_atomic_mass,
#             density=density
#         )
#
#         attenuation_length = calculate_attenuation_length(
#             lattice_parameter=lattice_parameter,
#             avg_atomic_number=average_atomic_number,
#             kinetic_energy=kinetic_energy_dict[chemical_formula]
#         )
#
#         sensitivity_factor = calculate_sensitivity_factor(
#             differential_cross_section=differential_cross_section,
#             attenuation_length=attenuation_length
#         )
#
#         sensitivity_factor_dict[chemical_formula] = sensitivity_factor
#
#     return sensitivity_factor_dict


def create_sensitivity_factor_dict(area_dict: dict, kinetic_energy_dict: dict, core_level: str, analyzer_angle: float,
                                   density_dict: dict, period_table_information: dict,
                                   x_ray_source_energy: float) -> dict:
    """
    This function creates a dictionary of sensitivity factors for each species in present in a core level

    :param kinetic_energy_dict: dict (Dictionary that contains the kinetic energies of each peak)
    :param core_level: str (Core level that is being analyzed. Necessary for calculating cross section and asymmetry
    parameters)
    :param analyzer_angle: float (Angle between the source and the analyzer. Used in the calculation of the
    differential cross section)
    :param density_dict: dict (Dictionary containing the densities of chemical species present on the sample. Used in
    the calculation of the lattice parameter which is need to calculate the x-ray attenuation length)
    :param period_table_information: dict (Dictionary containing atomic numbers, atomic masses, and element symbols
    for each element)
    :param x_ray_source_energy: float (Energy of your x-ray source)
    :return: sensitivity_factor_dict: dict (Dictionary containing sensitivity factors for each species on the surface
    of the sample for a given core level)
    """

    df = import_database_to_df(core_level)

    sensitivity_factor_dict = cp.deepcopy(kinetic_energy_dict)

    maximum_intensity_chemical_formula = max(area_dict, key=area_dict.get)
    density_chemical_formula = maximum_intensity_chemical_formula.split(" ")[0]

    density = density_dict[density_chemical_formula]
    split_formula = split_chemical_formula(chemical_formula=maximum_intensity_chemical_formula)

    average_atomic_mass = calculate_average_atomic_mass(
        chemical_formula_components=split_formula,
        periodic_table_information_dict=period_table_information
    )

    average_atomic_number = calculate_average_atomic_number(
        chemical_formula_components=split_formula,
        periodic_table_information_dict=period_table_information
    )

    cross_section = calculate_cross_section(
        x_ray_source_energy=x_ray_source_energy,
        df=df
    )

    asymmetry_parameter = calculate_asymmetry_parameter(
        x_ray_source_energy=x_ray_source_energy,
        df=df
    )

    differential_cross_section = calculate_differential_cross_section(
        cross_section=cross_section,
        asymmetry_parameter=asymmetry_parameter,
        analyzer_angle=analyzer_angle
    )

    lattice_parameter = calculate_lattice_parameter(
        average_atomic_mass=average_atomic_mass,
        density=density
    )

    attenuation_length = calculate_attenuation_length(
        lattice_parameter=lattice_parameter,
        avg_atomic_number=average_atomic_number,
        kinetic_energy=kinetic_energy_dict[maximum_intensity_chemical_formula]
    )

    print("sigma =", cross_section)
    print("beta =", asymmetry_parameter)
    print("average atomic mass=", average_atomic_mass)
    print("diff cross sec=", differential_cross_section)
    print("lattice_parameter=", lattice_parameter)
    print("attenuation length=", attenuation_length)

    sensitivity_factor = calculate_sensitivity_factor(
        differential_cross_section=differential_cross_section,
        attenuation_length=attenuation_length
    )

    sensitivity_factor_dict[maximum_intensity_chemical_formula] = sensitivity_factor

    return sensitivity_factor_dict


def import_to_df(folder_path: str) -> list:
    """
    This function takes in the file path to a folder which contains the exported CasaXPS data files for a single
    core level. It will pull the data from the last exported file, allowing you to track iterations of files into
    the same folder. From the data, it will extract the positions, FWHMs, and areas of each peak and return it in
    a list of pd.DataFrames

    :param folder_path: str (path to folder that contains the exported CasaXPS files for a certain for level)
    :return: [position_row, FWHM_row, area_row] list (list of three pd.DataFrames containing the information
             about peak positions, peak FWHM's, and peak areas)
    """

    # For windows only, replace the \\ with /
    folder_name = folder_path.replace(os.sep, '/')

    # List the file names in the folder ex. C1s.txt, O1s.txt, ...
    input_file_names = os.listdir(folder_name)

    # If there are multiple files, it takes the last one
    last_export_file_name = input_file_names[-1]

    # Generate the paths to each file
    file_path = folder_name + '/' + last_export_file_name

    # Read in the files
    specta_df = pd.read_csv(file_path, delimiter='\t', skiprows=0, nrows=4)

    # Drop NaN columns
    removed_NA_df = specta_df.dropna(axis=1, inplace=False)

    # Drop the name column
    removed_name_df = removed_NA_df.set_index(['Name'], inplace=False)

    # Extract the position, FWHM, and area rows
    position_row = removed_name_df.iloc[0]
    FWHM_row = removed_name_df.iloc[1]
    area_row = removed_name_df.iloc[2]

    # Return the position, FWHM, and area rows
    return [position_row, FWHM_row, area_row]


def import_xps_data(folder_path: str) -> list:
    """
    This function can take a file path that contains folders for each experiment included in the analysis

    :param folder_path: str (path to the folder that contains all of the experiments)
    :return: [position_dict, fwhm_dict, area_dict] list[dict[str:dict[str:float]]] (dictionaries that contain the
             position, fwhm, and area information for each experiment)
    """

    # Create empty dictionaries for Position, FWHM, and Area
    position_dict = {}
    fwhm_dict = {}
    area_dict = {}

    # For widows only, replace the \\ with /
    backslash_folder_path = folder_path.replace(os.sep, '/')

    # Collect the folder names that contain each experiment
    names_array = os.listdir(backslash_folder_path)

    # Create the paths to the subfolders
    sub_folder_path_list = [backslash_folder_path + '/' + Name for Name in names_array]

    # For loop to loop through each subfolder and make a subdict out of
    for (i, path) in enumerate(sub_folder_path_list):

        # Import the dataframes
        position, fwhm, area = import_to_df(path)

        # Create a dictionary key that is the experiment name, and make its value a sub dictionary
        # that contains all the  data from the experiment
        position_dict[f'{names_array[i]}'] = {column_name: value for (column_name, value)
                                              in zip(position.index, position.__array__(dtype=float))}
        fwhm_dict[f'{names_array[i]}'] = {column_name: value for (column_name, value)
                                          in zip(fwhm.index, fwhm.__array__(dtype=float))}
        area_dict[f'{names_array[i]}'] = {column_name: value for (column_name, value)
                                          in zip(area.index, area.__array__(dtype=float))}

    # Return the dictionaries
    return [position_dict, fwhm_dict, area_dict]


def Normalize_Areas_Dict(Area_Dict: dict, Position_Dict: dict, Analyzer_Angle, Unwanted_Species):

    a = 1
    Normalized_Area_Dict = {}
    Core_Level_List = Area_Dict.keys()
    for Core_Level in Core_Level_List:
        for Species in Unwanted_Species:
            if Species in Area_Dict[Core_Level].keys():
                Area_Dict[Core_Level].pop(Species)
        Normalized_Area_Dict[Core_Level] = {Key: (Value / Sensitivity_Factor_Dict[Analyzer_Angle][Core_Level])
                                            for (Key, Value) in
                                            zip(Area_Dict[Core_Level].keys(), Area_Dict[Core_Level].values())}
    return Normalized_Area_Dict


# file_path = "test_folder/Files_For_Plotting"
#
# position, FWHM, areas = import_xps_data(file_path)
#
# KE_dict = calculate_kinetic_energy_dict(
#     position_dict={'TiO2': 36.5},
#     x_ray_source_energy=1253.6,
#     detector_work_function=4.6
# )


# print(areas['W4f'])
# sensitivity_factor_dict = create_sensitivity_factor_dict(
#     area_dict=areas['W4f'],
#     kinetic_energy_dict=KE_dict,
#     core_level='W4f',
#     analyzer_angle=60,
#     density_dict={'WS2': 7500},
#     period_table_information=periodic_table,
#     x_ray_source_energy=1487
# )

# sensitivity_factor_dict = create_sensitivity_factor_dict(
#     area_dict={'TiO2': 1000, 'something else': 100},
#     kinetic_energy_dict=KE_dict,
#     core_level='Ti3p',
#     analyzer_angle=55,
#     density_dict={'TiO2': 4240},
#     period_table_information=periodic_table,
#     x_ray_source_energy=1253.6
# )
# print(sensitivity_factor_dict)
