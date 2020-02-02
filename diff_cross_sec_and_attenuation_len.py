import streamlit as st
import XPSAnalysis as xps
import matplotlib.pyplot as plt
import pandas as pd


def differential_cross_section_equation():
    st.header('Differential Cross Section Calculation')

    st.markdown('''
    Differential cross section are calculated using the below equation. $^{\\text{1}}$ Cross section and asymmetry parameters have been acquired from: $^{\\text{3-5}}$ https://vuo.elettra.eu/services/elements/WebElements.html
    ''')

    st.latex(
        "\\frac{d\\sigma_{i}}{d\\Omega_{i}} = \\frac{\\sigma_{i}}{4\\pi} \\left[1 + \\frac{\\beta_{i}}{2} \\left(\\frac{3}{2} \\sin^{2}(\\alpha) - 1\\right)\\right]")

    st.latex(
        "\\frac{d\\sigma_{i}}{d\\Omega_{i}} = \\text{ Differential cross section}"
    )

    st.latex(
        "\\sigma_{i} = \\text{ Cross section}"
    )

    st.latex(
        "\\beta_{i} = \\text{ Asymmetry parameter}"
    )

    st.latex(
        "\\alpha = \\text{ Analyzer angle}"
    )

    return ''


def x_ray_attenuation_length_equation():

    st.header('Effective Attenuation Length Calculation')

    st.markdown("""
        Effective attenuation lengths are calculated using the below equation: $^{2}$
    """)

    st.latex(
        "\\Lambda_{i} = 0.316a^{3/2}\\left[\\frac{E}{Z^{0.45}\\left[\\ln\\left(\\frac{E}{27}\\right) + 3\\right]} + 4\\right]")

    st.latex(
        "\\Lambda_{i} = \\text{ Effective attenuation length  [nm]}"
    )

    st.latex(
        "a = \\text{ Lattice parameter  [nm]}"
    )

    st.latex(
        "E = \\text{ Kinetic energy  [eV]}"
    )

    st.latex(
        "Z = \\text{ Average atomic number}"
    )

    st.markdown(' ')
    st.markdown(' ')
    st.markdown("The lattice parameter can be calculated using the following equation: $^{2}$")

    st.latex(
        "a = 10^{8} \\left(\\frac{\\mu}{\\rho N_{A}}\\right)^{1/3}"
    )

    st.latex(
        "a = \\text{ Lattice parameter  [nm]}"
    )

    st.latex(
        "\\mu = \\text{ Average atomic mass  [g/mol]}"
    )

    st.latex(
        "\\rho = \\text{ Density  [kg/m}^{3}\\text{]}]"
    )

    st.latex(
        "N_{A} = \\text{ Avagadro's number}"
    )

    return ''


def info_section():

    st.title('About the App')

    st.markdown("""

        **Author: ** Derek Dardzinski  
        **Contact: ** dardzinski.derek@gmail.com


        This web application is a tool to assist X-ray photoelectron spectroscopy (XPS) users in the analysis of their spectra.

        The intent of this tool is to make XPS related information easier to access, and to reduce the work required to calculate important parameters such as the differential cross section and the effective attenuation length in a material.

        This application has been created in Python using [Streamlit](https://www.streamlit.io/) and deployed through [Heroku](https://www.heroku.com/apps). All code is open source and can be found on [GitHub](https://github.com/DerekDardzinski/xps_analysis).

        ### References

        [1] Bertrand, Alain, et al. “Atomic Calculation of Photoionization Cross-Sections and Asymmetry Parameters.” WebCrossSections, vuo.elettra.eu/services/elements/WebElements.html.

        [2] Cumpson, Peter J., and Martin P. Seah. "Elastic scattering corrections in AES and XPS. II. Estimating attenuation lengths and conditions required for their valid use in overlayer/substrate experiments." Surface and Interface Analysis: An International Journal devoted to the development and application of techniques for the analysis of surfaces, interfaces and thin films 25.6 (1997): 430-446.

        [3] Fadley, C. S. "Basic concepts of X-ray photoelectron spectroscopy." Electron spectroscopy: theory, techniques and applications 2 (1978): 1-156.

        [4] Yeh, J. J., and I. Lindau. "Atomic subshell photoionization cross sections and asymmetry parameters: 1⩽ Z⩽ 103." Atomic data and nuclear data tables 32.1 (1985): 1-155.

        [5] Yeh, J. J. "Atomic Calculation of Photoionization Cross-Section and Asymmetry Parameters Gordon and Breach Science Publishers." Langhorne, PE (USA) (1993).

    """)

    return ''


def import_and_plot_data(element_core_level, df):

    st.dataframe(df)

    font_size = 9

    fig = plt.figure(dpi=400, figsize=(6, 3))
    ax1 = fig.add_subplot(121)
    ax1.tick_params(axis='both', labelsize=font_size)
    plt.ylabel('Cross Section, $\\alpha$', fontsize=font_size)
    plt.xlabel('Energy [eV]', fontsize=font_size)

    plt.minorticks_on()
    plt.grid(
        True,
        which='major',
        axis="both",
        alpha=0.3
    )

    plt.grid(
        True,
        which='minor',
        axis="both",
        alpha=0.1,
        linestyle="--"
    )

    ax2 = fig.add_subplot(122)
    ax2.tick_params(axis='both', labelsize=font_size)
    plt.ylabel('Asymmetry Parameter, $\\beta$', fontsize=font_size)
    plt.xlabel('Energy [eV]', fontsize=font_size)

    plt.minorticks_on()
    plt.grid(
        True,
        which='major',
        axis="both",
        alpha=0.3
    )

    plt.grid(
        True,
        which='minor',
        axis="both",
        alpha=0.1,
        linestyle="--"
    )

    ax1.plot(
        df['Energy [eV]'].__array__(),
        df['Cross Section'],
        color='black',
        linestyle='',
        linewidth=0.5,
        alpha=1,
        marker='x',
        markersize=3
    )
    ax1.plot(
        df['Energy [eV]'].__array__(),
        df['Cross Section'],
        color='black',
        linestyle='-',
        linewidth=0.5,
        alpha=0.6,
    )
    ax2.plot(
        df['Energy [eV]'].__array__(),
        df['Asymmetry Parameter'],
        color='black',
        linestyle='',
        linewidth=0.5,
        alpha=1,
        marker='x',
        markersize=3
    )
    ax2.plot(
        df['Energy [eV]'].__array__(),
        df['Asymmetry Parameter'],
        color='black',
        linestyle='-',
        linewidth=0.5,
        alpha=0.6,
    )
    plt.tight_layout()
    st.pyplot()

    return ''


def import_data_and_calculate_differential_cross_section(element_core_level, x_ray_source_energy, analyzer_angle):
    df = xps.import_database_to_df(element_core_level)

    crosss_section = xps.calculate_cross_section(x_ray_source_energy, df)

    asymmetry_parameter = xps.calculate_asymmetry_parameter(x_ray_source_energy, df)

    differential_cross_section = xps.calculate_differential_cross_section(
        cross_section=crosss_section,
        asymmetry_parameter=asymmetry_parameter,
        analyzer_angle=analyzer_angle
    )

    return df, crosss_section, asymmetry_parameter, differential_cross_section


def calculate_x_ray_attenuation_length_from_formula(chemical_formula, period_table_info, material_density, x_ray_source_energy, binding_energy, analyzer_work_function):
    split_formula = xps.split_chemical_formula(chemical_formula)

    avg_atomic_mass = xps.calculate_average_atomic_mass(
        chemical_formula_components=split_formula,
        periodic_table_information_dict=period_table_info
    )

    avg_atomic_number = xps.calculate_average_atomic_number(
        chemical_formula_components=split_formula,
        periodic_table_information_dict=period_table_info
    )

    lattice_parameter = xps.calculate_lattice_parameter(
        average_atomic_mass=avg_atomic_mass,
        density=material_density
    )

    attenuation_length = xps.calculate_attenuation_length(
        lattice_parameter=lattice_parameter,
        avg_atomic_number=avg_atomic_number,
        kinetic_energy=x_ray_source_energy - binding_energy - analyzer_work_function
    )

    return attenuation_length, lattice_parameter, avg_atomic_number, avg_atomic_mass


differential_cross_section_title = "Differential Cross Section"
x_ray_attenuation_length_title = "Effective Attentuation Length"
info_title = "About the App / References"

option = st.sidebar.radio(
    "Select an Option", [differential_cross_section_title, x_ray_attenuation_length_title, info_title], index=0)

# ========================================
# ============= Info Page ================
# ========================================

if option == info_title:
    st.write(info_section())

# ========================================================
# ========== Differential Cross Section Page =============
# ========================================================

if option == differential_cross_section_title:

    st.title('Differential Photoionization Cross Section Calculator')

    show_differential_cross_section_equation = st.checkbox('Show Equation')

    if show_differential_cross_section_equation is True:
        st.write(differential_cross_section_equation())

    element_core_level = st.text_input(
        label='Enter element symbol and core level, eg. C1s, O1s, Sn3d, Mo3d',
        value=''
    )

    analyzer_angle = st.number_input(
        'Input analyzer angle in degrees',
        min_value=0.0,
        max_value=180.0,
        value=60.0
    )

    x_ray_source_energy = st.number_input(
        'Input X-ray source energy [eV]',
        min_value=0.0,
        value=1486.7
    )

    if element_core_level != '':

        df, cross_section, asymmetry_parameter, differential_cross_section = import_data_and_calculate_differential_cross_section(
            element_core_level=element_core_level,
            x_ray_source_energy=x_ray_source_energy,
            analyzer_angle=analyzer_angle
        )

        results_dict = {'Cross Section': cross_section, 'Asymmetry Parameter': asymmetry_parameter,
                        'Differential Cross Section': differential_cross_section}

        st.markdown(f"""
            | Variable | Value | Units |
            | :------: | :---: | :---: |
            | Differential Cross Section | {round(differential_cross_section, 5)} | barns |
            | Cross Section | {round(cross_section, 5)} | barns |
            | Asymmetry Parameter | {round(asymmetry_parameter, 5)} | - |
        """)

        st.markdown(' ')

    show_data = st.checkbox('Show cross section and asymmetry parameter tables and charts')

    if element_core_level != '' and show_data is True:
        st.write(import_and_plot_data(element_core_level, df))

# ========================================================
# ========== X-ray Attenuation Length Page ===============
# ========================================================

if option == x_ray_attenuation_length_title:

    st.title('Effective Attenuation Length Calculator')

    show_equation = st.checkbox('Show Equations')

    if show_equation is True:
        st.write(x_ray_attenuation_length_equation())

    period_table_info = xps.import_periodic_table_information(
        file_path='https://raw.githubusercontent.com/DerekDardzinski/xps_analysis/master/PeriodicTableData.txt')

    chemical_formula = st.text_input(
        'Enter the chemical formula of the compound using whole numbers for, eg. SnO2, MoS2, AlO3',
        value=''
    )

    material_density = st.number_input(
        'Enter the density of the material [kg/m^3]',
        min_value=0.0,
        value=1000.0
    )

    x_ray_source_energy = st.number_input(
        'Input X-ray source energy [eV]',
        min_value=0.0,
        value=1486.7
    )

    binding_energy = st.number_input(
        'Binding Energy of the XPS peak [eV]'
    )

    analyzer_work_function = st.number_input(
        'Work function of XPS analyzer [eV]',
        value=4.6
    )

    if chemical_formula != '':
        attenuation_length, lattice_parameter, avg_atomic_number, avg_atomic_mass = calculate_x_ray_attenuation_length_from_formula(
            chemical_formula, period_table_info, material_density, x_ray_source_energy, binding_energy, analyzer_work_function)

        st.markdown(f"""
            | Variable | Value | Units |
            | :------: | :---: | :---: |
            | Attenuation Length | {round(attenuation_length, 3)} | nm |
            | Lattice Parameter | {round(lattice_parameter, 3)} | nm |
            | Average Atomic Mass | {round(avg_atomic_mass, 3)} | g/mol |
            | Average Atomic Number | {round(avg_atomic_number, 3)} | - |
        """)
