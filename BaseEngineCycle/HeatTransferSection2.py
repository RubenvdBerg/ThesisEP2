import itertools

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from BaseEngineCycle.ThrustChamber import ThrustChamber
from dataclasses import dataclass, replace, field
from BaseEngineCycle.FlowState import FlowState, ManualFlowState, CoolantFlowState
from typing import Optional, Callable
from math import pi
from CoolProp import CoolProp
import warnings
from functools import cached_property

@dataclass
class HeatExchanger:
    thrust_chamber: ThrustChamber
    combustion_chamber_flow_state: ManualFlowState
    coolant_inlet_flow_state: FlowState

    number_of_coolant_channels: float
    radiative_factor: float
    chamber_wall_conductivity: float
    chamber_wall_thickness: float

    max_iterations: float = 10
    amount_of_sections: float = 100
    counter_flow: bool = True
    post_injection_build_up_ratio: float = .25  # [-]
    coolant_coefficient_roughness_correction = 1
    coolant_coefficient_fin_correction = 1
    _coolant_channel_diameter: Optional[float] = None
    _recovery_factor: Optional[float] = None  # [-]
    _initial_flow_speed: float = 10  # [m/s]
    verbose: bool = True
    iteration_accuracy: float = 1e-5
    data: dict = field(init=False)

    # Iterative section variables, not required at init
    section_hot_side_wall_temp: float = field(init=False, repr=False)
    section_cold_side_wall_temp: float = field(init=False, repr=False)
    section_distance_from_throat: float = field(init=False, repr=False)
    next_section_distance_from_throat: float = field(init=False, repr=False)
    section_coolant_total_temp: float = field(init=False, repr=False)
    section_coolant_total_pressure: float = field(init=False, repr=False)
    section_coolant_state: CoolantFlowState = field(init=False, repr=False)
    section_heat_flux: float = field(init=False, repr=False)
    section_hot_gas_adiabatic_wall_temp: float = field(init=False, repr=False)
    section_hot_gas_film_temp: float = field(init=False, repr=False)
    section_hot_gas_static_temp: float = field(init=False, repr=False)

    def __post_init__(self):
        # Set values at start of coolant channel
        self.section_hot_side_wall_temp = 400 if self.counter_flow else self.coolant_inlet_flow_state.temperature
        self.section_cold_side_wall_temp = self.coolant_inlet_flow_state.temperature
        self.section_coolant_total_temp = self.coolant_inlet_flow_state.temperature
        self.section_coolant_total_pressure = self.coolant_inlet_flow_state.pressure

        self.init_data()
        self.iterate()

    def init_data(self):
        self.data = {'Distance from Throat [m]': [],
                     'Temperature [K]': {'CoolantBulk': [],
                                         'HotSideWall': [],
                                         'ColdSideWall': [],
                                         'HotGasAdiabatic': [],
                                         'HotGasStatic': [],
                                         'HotGasFilm': [], },
                     'Heat-Transfer Coefficient [W/(K*m2]': {'HotSide': [],
                                                             'ColdSide': [],
                                                             'HotSideConv.': [], },
                     'Heat Flux [W/m2]': {'Total': [],
                                          'HotSide': [],
                                          'ColdSide': [],
                                          'HotSideConv.': [],
                                          'HotSideRad.': [], },
                     'Coolant State': {'T': [],
                                       'T0': [],
                                       'p': [],
                                       'p0': [],
                                       'rho': [],
                                       'cp': [],
                                       'T_b': []}
                     }

    def write_data(self):
        self.data['Distance from Throat [m]'].append(self.section_distance_from_throat)
        self.data['Temperature [K]']['CoolantBulk'].append(self.section_coolant_bulk_temp)
        self.data['Temperature [K]']['HotSideWall'].append(self.section_hot_side_wall_temp)
        self.data['Temperature [K]']['ColdSideWall'].append(self.section_cold_side_wall_temp)
        self.data['Temperature [K]']['HotGasAdiabatic'].append(self.section_hot_gas_adiabatic_wall_temp)
        self.data['Temperature [K]']['HotGasStatic'].append(self.section_hot_gas_static_temp)
        self.data['Temperature [K]']['HotGasFilm'].append(self.section_hot_gas_film_temp)
        self.data['Heat-Transfer Coefficient [W/(K*m2]']['HotSide'].append(
            self.section_hot_gas_total_heat_transfer_coefficient)
        self.data['Heat-Transfer Coefficient [W/(K*m2]']['ColdSide'].append(
            self.section_coolant_convective_heat_transfer_coefficient)
        self.data['Heat-Transfer Coefficient [W/(K*m2]']['HotSideConv.'].append(
            self.section_hot_gas_convective_heat_transfer_coefficient)
        self.data['Heat Flux [W/m2]']['Total'].append(self.section_heat_flux)
        self.data['Heat Flux [W/m2]']['HotSide'].append(self.section_hot_gas_total_heat_flux)
        self.data['Heat Flux [W/m2]']['ColdSide'].append(self.section_coolant_heat_flux)
        self.data['Heat Flux [W/m2]']['HotSideConv.'].append(self.section_hot_gas_convective_heat_flux)
        self.data['Heat Flux [W/m2]']['HotSideRad.'].append(self.section_hot_gas_radiative_heat_flux)
        self.data['Coolant State']['T_b'].append(self.section_coolant_bulk_temp)
        self.data['Coolant State']['T'].append(self.section_coolant_state.static_temperature)
        self.data['Coolant State']['T0'].append(self.section_coolant_total_temp)
        self.data['Coolant State']['p'].append(self.section_coolant_state.static_pressure)
        self.data['Coolant State']['p0'].append(self.section_coolant_total_pressure)

    @cached_property
    def section_axial_length(self):
        return self.thrust_chamber.length / self.amount_of_sections

    @cached_property
    def coolant_channel_area(self):
        return (pi / 4) * self.coolant_channel_diameter ** 2

    @cached_property
    def mass_flow_per_coolant_channel(self):
        return self.coolant_inlet_flow_state.mass_flow / self.number_of_coolant_channels

    @cached_property
    def coolant_channel_mass_flux(self):
        return self.mass_flow_per_coolant_channel / self.coolant_channel_area

    @cached_property
    def combustion_temp(self):
        return self.combustion_chamber_flow_state.temperature

    @cached_property
    def recovery_factor(self):
        if self._recovery_factor is None:
            # Zandbergen 2017 p.160, experimental relationship for turbulent boundary layers
            return self.combustion_chamber_flow_state.prandtl_number ** (1 / 3)
        else:
            return self._recovery_factor

    @cached_property
    def section_pressure_change(self):
        return -1 * self.combustion_chamber_flow_state.pressure * .2 / self.amount_of_sections

    @property
    def percentage_iterated(self):
        return f'{100 * (self.section_distance_from_throat - self.thrust_chamber.min_distance_from_throat) / self.thrust_chamber.length:.2f}%'

    @cached_property
    def coolant_channel_diameter(self):
        if self._coolant_channel_diameter is None:
            chosen_flow_speed = self._initial_flow_speed
            initial_density = self.coolant_inlet_flow_state.density
            mass_flow = self.coolant_inlet_flow_state.mass_flow
            volume_flow = mass_flow / initial_density
            required_area = volume_flow / chosen_flow_speed
            channel_area = required_area / self.number_of_coolant_channels
            return 2 * (channel_area / pi) ** .5
        else:
            return self._coolant_channel_diameter

    @cached_property
    def total_coolant_flow_area(self):
        return self.coolant_channel_area * self.number_of_coolant_channels

    def iterate(self):
        """
        Loop over the nozzle contour, while calculating the heat transfer variables at every section and adding them
        to the data dictionary.

        The process is started at the injector face normally, if coolant flows in a direction counter to
        the combustion gas the loop is started from the nozzle exit side instead.
        """
        if self.counter_flow:
            linspace_input = self.thrust_chamber.throat_distance_tuple[::-1]
        else:
            linspace_input = self.thrust_chamber.throat_distance_tuple

        axial_distances = np.linspace(*linspace_input, int(self.amount_of_sections), endpoint=False)

        for i, (x, x_next) in enumerate(zip(axial_distances, axial_distances[1:])):
            self.section_distance_from_throat = x
            self.next_section_distance_from_throat = x_next
            self.iterate_coolant_dynamic_state()
            print(self.percentage_iterated)
            self.iterate_wall_temps()
            self.write_data()
            self.set_next_state()

        if self.verbose:
            print('Heat Transfer Iteration Complete')

    def iterate_coolant_dynamic_state(self):
        """Dynamic state iteration happens internally in CoolantFlowState"""
        self.section_coolant_state = CoolantFlowState(propellant_name=self.coolant_inlet_flow_state.propellant_name,
                                                      total_temperature=self.section_coolant_total_temp,
                                                      total_pressure=self.section_coolant_total_pressure,
                                                      mass_flow=self.coolant_inlet_flow_state.mass_flow,
                                                      type=self.coolant_inlet_flow_state.type,
                                                      _iteration_accuracy=self.iteration_accuracy,
                                                      verbose=self.verbose,
                                                      total_flow_area=self.total_coolant_flow_area)

    def iterate_wall_temps(self):
        """Loop until hot and cold side wall temperatures converge"""
        th_w = self.chamber_wall_thickness
        k_w = self.chamber_wall_conductivity
        t_b = self.section_coolant_bulk_temp
        self.set_hot_gas_temps()
        t_ad = self.section_hot_gas_adiabatic_wall_temp

        def get_q(h_hg, h_cool):
            """Get heat flux between hot gas and coolant"""
            return (t_ad - t_b) / (1 / h_hg + th_w / k_w + 1 / h_cool)

        def get_wall_temp_hot(q_hg, h_hg):
            return t_ad - q_hg / h_hg if h_hg else self.section_hot_side_wall_temp

        def get_wall_temp_cold(q_cool, h_cool):
            return t_b + q_cool / h_cool if h_cool else self.section_cold_side_wall_temp

        h_g = self.section_hot_gas_total_heat_transfer_coefficient
        h_c = self.section_coolant_convective_heat_transfer_coefficient
        q = get_q(h_g, h_c)
        new_wall_temp_hot = get_wall_temp_hot(q, h_g)
        new_wall_temp_cold = get_wall_temp_cold(q, h_c)

        not_ran = True
        while (self.error_too_large(self.section_hot_side_wall_temp, new_wall_temp_hot)
               or self.error_too_large(self.section_cold_side_wall_temp, new_wall_temp_cold)
               or not_ran):
            not_ran = False
            self.section_hot_side_wall_temp = new_wall_temp_hot
            self.section_cold_side_wall_temp = new_wall_temp_cold
            self.set_film_temp()
            new_h_g = self.section_hot_gas_total_heat_transfer_coefficient
            new_h_c = self.section_coolant_convective_heat_transfer_coefficient
            new_q = get_q(new_h_g, new_h_c)
            new_wall_temp_hot = get_wall_temp_hot(new_q, new_h_g)
            new_wall_temp_cold = get_wall_temp_cold(new_q, new_h_c)

            if self.verbose:
                print(f'Hot Wall temp [K]-> \n'
                      f'Current: {self.section_hot_side_wall_temp:.4e}, \n'
                      f'Expected:{new_wall_temp_hot:.4e}')
                print(f'Cold Wall temp [K]-> \n'
                      f'Current: {self.section_cold_side_wall_temp:.4e}, \n'
                      f'Expected:{new_wall_temp_cold:.4e}\n\n')

            self.section_hot_side_wall_temp = new_wall_temp_hot
            self.section_cold_side_wall_temp = new_wall_temp_cold

        self.section_heat_flux = new_q
        self.set_hot_gas_temps()

    def error_too_large(self, current: float, expected: float):
        error = abs((current - expected) / expected)
        return error > self.iteration_accuracy

    def set_hot_gas_temps(self):
        y = self.combustion_chamber_flow_state.heat_capacity_ratio
        m = self.thrust_chamber.get_mach(self.section_distance_from_throat)
        t_c = self.combustion_temp
        f_r = self.recovery_factor
        f_red = self.section_reduction_factor
        t_hg_st = f_red * t_c / (1 + (y - 1) / 2 * m ** 2)
        t_hg_ad = t_hg_st * (1 + (y - 1) / 2 * m ** 2 * f_r)
        self.section_hot_gas_static_temp = t_hg_st
        self.section_hot_gas_adiabatic_wall_temp = t_hg_ad
        self.set_film_temp()

    def set_film_temp(self):
        t_w = self.section_hot_side_wall_temp
        t_ad = self.section_hot_gas_adiabatic_wall_temp
        t_st = self.section_hot_gas_static_temp
        self.section_hot_gas_film_temp = t_st * .28 + t_w * .5 + t_ad * .22

    @property
    def section_reduction_factor(self):
        """Factor that reduces the heat flow until maximum heat flux build up after the injector is achieved"""
        len_cc = self.thrust_chamber.chamber.length
        dist_min = self.thrust_chamber.min_distance_from_throat
        r_build_up = self.post_injection_build_up_ratio
        # Distance from injector divided by total chamber length
        inj_distance_ratio = (self.section_distance_from_throat - dist_min) / len_cc
        return 1
        if inj_distance_ratio < r_build_up:
            return (self.section_distance_from_throat - dist_min) / (r_build_up * len_cc)
        else:
            return 1

    @property
    def section_hot_gas_convective_heat_transfer_coefficient(self):
        mf = self.combustion_chamber_flow_state.mass_flow
        di = self.thrust_chamber.get_radius(self.section_distance_from_throat) * 2
        mu = self.combustion_chamber_flow_state.dynamic_viscosity
        cp = self.combustion_chamber_flow_state.specific_heat_capacity
        pr = self.combustion_chamber_flow_state.prandtl_number
        t0 = self.section_hot_gas_static_temp
        tf = self.section_hot_gas_film_temp
        return 0.026 * 1.213 * mf ** .8 * di ** -1.8 * mu ** .2 * cp * pr ** -.6 * (t0 / tf) ** .68

    @property
    def section_hot_gas_total_heat_transfer_coefficient(self):
        return self.section_hot_gas_convective_heat_transfer_coefficient * (1 + self.radiative_factor)

    @property
    def section_hot_gas_convective_heat_flux(self):
        return self.section_hot_gas_convective_heat_transfer_coefficient * (self.section_hot_gas_adiabatic_wall_temp
                                                                            - self.section_hot_side_wall_temp)

    @property
    def section_hot_gas_radiative_heat_flux(self):
        return self.section_hot_gas_convective_heat_flux * self.radiative_factor

    @property
    def section_hot_gas_total_heat_flux(self):
        return self.section_hot_gas_convective_heat_flux + self.section_hot_gas_radiative_heat_flux

    @property
    def section_coolant_convective_heat_transfer_coefficient(self):
        k_c = self.section_coolant_state.conductivity
        D = self.coolant_channel_diameter
        re = self.section_coolant_state.get_reynolds(
            linear_dimension=self.coolant_channel_diameter,
            flow_speed=self.coolant_channel_mass_flux / self.section_coolant_state.density)
        pr = self.section_coolant_state.prandtl_number
        t_b = self.section_coolant_bulk_temp
        t_w = self.section_cold_side_wall_temp
        ksi = self.coolant_coefficient_roughness_correction
        ksi2 = self.coolant_coefficient_fin_correction
        # print(f'kc:{k_c}, D:{D}, re:{re}, pr:{pr}, tb:{t_b}, tw:{t_w}')
        return 0.025 * k_c / D * re ** .8 * pr ** .4 * (t_b / t_w) ** .55 * ksi * ksi2

    @property
    def section_coolant_bulk_temp(self):
        return self.section_coolant_state.static_temperature

    @property
    def section_coolant_heat_flux(self):
        return self.section_coolant_convective_heat_transfer_coefficient * (self.section_cold_side_wall_temp
                                                                            - self.section_coolant_bulk_temp)

    def get_section_nozzle_surface(self):
        # noinspection PyTupleAssignmentBalance
        section_surface, _ = integrate.quad(func=lambda x: 2 * pi * self.thrust_chamber.get_radius(x),
                                            a=self.section_distance_from_throat,
                                            b=self.next_section_distance_from_throat)
        return abs(section_surface)

    def set_next_state(self):
        """Calculate change in temperature and pressure and update """
        # noinspection PyTupleAssignmentBalance
        A = self.get_section_nozzle_surface()
        q = self.section_heat_flux
        m_dot = self.coolant_inlet_flow_state.mass_flow
        h_in = self.section_coolant_state.total_mass_specific_enthalpy
        t_in = self.section_coolant_state.total_temperature

        Q = A * q  # Heat Transfer
        delta_h = Q / m_dot  # Enthalpy Increase

        h_out = h_in + delta_h
        t_out = CoolProp.PropsSI('T',
                                 'H', h_out,
                                 'P', self.section_coolant_state.total_pressure,
                                 self.section_coolant_state.coolprop_name)

        delta_t = t_out - t_in
        self.section_coolant_total_temp += delta_t
        delta_p = self.section_pressure_change
        self.section_coolant_total_pressure += delta_p

    # Plot Stuff Below
    @property
    def distances(self):
        return self.data['Distance from Throat [m]']

    def plot_values(self, variable: str, extra_func: Optional[Callable] = None):
        fig, ax = plt.subplots()

        self.plot_nozzle_contour_background(ax=ax)

        if extra_func is not None:
            extra_func(ax)

        for key, value in self.data[variable].items():
            ax.plot(self.distances, value, label=key)

        ax.legend()
        ax.set_ylabel(variable)
        ax.set_xlabel('Distance from Throat [m]')
        ax.set_title(variable.split('[')[0].strip(' '))
        plt.show()

    def plot_temps(self):
        self.plot_values('Temperature [K]')

    def plot_coeffs(self):
        self.plot_values('Heat-Transfer Coefficient [W/(K*m2]')

    def plot_flux(self, **kwargs):
        self.plot_values('Heat Flux [W/m2]', **kwargs)

    def plot_coolant(self):
        c_data = self.data['Coolant State']
        fig, ax = plt.subplots()

        self.plot_nozzle_contour_background(ax=ax)

        ax2 = ax.twinx()
        for variable, axis, color in zip(['T', 'T0', 'p', 'p0'], [ax, ax, ax2, ax2],
                                         ['orange', 'red', 'blue', 'lightblue']):
            axis.plot(self.distances, c_data[variable], label=variable, color=color)

        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc=0)
        ax.set_ylabel('Temperature [K]')
        ax2.set_ylabel('Pressure [Pa]')
        ax.set_xlabel('Distance from Throat [m]')
        ax.set_title(variable.split('[')[0].strip(' '))
        plt.savefig(rf'savefigs\plot_{self.max_iterations}')
        plt.show()

    def plot_perakis_data(self, ax: plt.Axes):
        import pandas as pd
        filename = r'..\BaseEngineCycle\Data\Perakis2021_HeatTransfer'
        file1 = np.genfromtxt(fname=filename + '.txt', delimiter=',')
        distances_exp = file1[1:, 0] * 1e-3 + self.thrust_chamber.min_distance_from_throat
        values_exp = file1[1:, 1] * 1e6
        ax.plot(distances_exp, values_exp, label='q_Perakis', color='pink', linestyle='--')

    def plot_nozzle_contour_background(self, ax: plt.Axes):
        ax3 = ax.twinx()
        ax3.plot(self.distances, [self.thrust_chamber.get_radius(x) for x in self.distances], color='grey',
                 linestyle='--')
        ax3.get_yaxis().set_visible(False)

    def plot_all(self):
        self.plot_coeffs()
        self.plot_temps()
        self.plot_flux()
        self.plot_coolant()
