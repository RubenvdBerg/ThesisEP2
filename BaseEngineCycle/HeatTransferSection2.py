import itertools

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
from scipy import integrate
from BaseEngineCycle.ThrustChamber import ThrustChamber
from dataclasses import dataclass, replace, field
from BaseEngineCycle.FlowState import FlowState, ManualFlowState, CoolantFlowState
from typing import Optional, Callable
from math import pi, log, tanh
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
    amount_of_sections: float = 286
    counter_flow: bool = False
    post_injection_build_up_ratio: float = .25  # [-]
    coolant_channel_roughness_height: float = 6e-6
    hot_gas_convective_heat_transfer_coefficient_mode: str = 'ModifiedBartz'
    _coolant_coefficient_roughness_correction: float = 1
    _coolant_coefficient_fin_correction: float = 1
    _coolant_channel_diameter: Optional[float] = None
    _recovery_factor: Optional[float] = None  # [-]
    _initial_flow_speed: float = 10  # [m/s]
    verbose: bool = True
    iteration_accuracy: float = 1e-3
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
    section_hot_gas_static_temp: float = field(init=False, repr=False)

    def __post_init__(self):
        # Set values at start of coolant channel
        self.section_hot_side_wall_temp = 400 if self.counter_flow else self.coolant_inlet_flow_state.temperature
        self.section_cold_side_wall_temp = self.coolant_inlet_flow_state.temperature
        self.section_coolant_total_temp = self.coolant_inlet_flow_state.temperature
        self.section_coolant_total_pressure = self.coolant_inlet_flow_state.pressure

        self.init_data()
        self.iterate()
        print(f'Max. T_hg_wall: {max(self.data["Temperature [K]"]["HotSideWall"]):.2f} K')
        print(f'Out. T_b      : {self.data["Temperature [K]"]["CoolantBulk"][-1]:.2f} K')
        print(
            f'DeltaP        : {(self.data["Coolant State"]["p"][-1] - self.data["Coolant State"]["p"][0]) * 1e-5 :.2f} bar')
        print(
            f'DeltaP0       : {(self.data["Coolant State"]["p0"][-1] - self.data["Coolant State"]["p0"][0]) * 1e-5 :.2f} bar')

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
            self.iterate_wall_temps()

            # if float(self.percentage_iterated.split('.')[0]) < 72:
            #     pass
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

        iterations = 0
        while True:
            h_g_n = self.section_hot_gas_convective_heat_transfer_coefficient
            h_g = self.section_hot_gas_total_heat_transfer_coefficient
            h_c = self.section_coolant_convective_heat_transfer_coefficient
            q = (t_ad - t_b) / (1 / h_g + th_w / k_w + 1 / h_c)
            self.section_heat_flux = q

            new_wall_temp_hot = t_ad - q / h_g if h_g else self.section_hot_side_wall_temp
            new_wall_temp_cold = t_b + q / h_c if h_c else self.section_cold_side_wall_temp

            if self.verbose:
                print(f'Hot Wall temp [K]-> \n'
                      f'Current: {self.section_hot_side_wall_temp:.4e}, \n'
                      f'Expected:{new_wall_temp_hot:.4e}')
                print(f'Cold Wall temp [K]-> \n'
                      f'Current: {self.section_cold_side_wall_temp:.4e}, \n'
                      f'Expected:{new_wall_temp_cold:.4e}\n\n')

            # Check if iteration is done
            iterations += 1
            error_hot = self.error_is_small(self.section_hot_side_wall_temp, new_wall_temp_hot)
            error_cold = self.error_is_small(self.section_cold_side_wall_temp, new_wall_temp_cold)
            if iterations > self.max_iterations or (error_hot and error_cold):
                break
            # Update variables
            self.section_hot_side_wall_temp = float(new_wall_temp_hot)
            self.section_cold_side_wall_temp = float(new_wall_temp_cold)

    def error_is_large(self, current: float, expected: float):
        error = abs((current - expected) / expected)
        return error > self.iteration_accuracy

    def error_is_small(self, current: float, expected: float):
        error = abs((current - expected) / expected)
        return error < self.iteration_accuracy

    def set_hot_gas_temps(self):
        y = self.combustion_chamber_flow_state.heat_capacity_ratio
        m = self.thrust_chamber.get_mach(self.section_distance_from_throat)
        # TODO: Remove
        m = 3
        t_c = self.combustion_temp
        f_r = self.recovery_factor
        f_red = self.section_reduction_factor
        t_hg_st = f_red * t_c / (1 + (y - 1) / 2 * m ** 2)
        t_hg_ad = t_hg_st * (1 + f_r * (y - 1) / 2 * m ** 2)
        self.section_hot_gas_static_temp = t_hg_st
        self.section_hot_gas_adiabatic_wall_temp = t_hg_ad

    @property
    def section_hot_gas_film_temp(self):
        t_w = self.section_hot_side_wall_temp
        t_ad = self.section_hot_gas_adiabatic_wall_temp
        t_st = self.section_hot_gas_static_temp
        t_f = t_st * .28 + t_w * .5 + t_ad * .22
        return t_f

    @property
    def section_reduction_factor(self):
        """Factor that reduces the heat flow until maximum heat flux build up after the injector is achieved"""
        len_cc = self.thrust_chamber.chamber.length
        dist_min = self.thrust_chamber.min_distance_from_throat
        r_build_up = self.post_injection_build_up_ratio
        # Distance from injector divided by total chamber length
        inj_distance_ratio = (self.section_distance_from_throat - dist_min) / len_cc
        return 1
        # TODO: Remove

        # if inj_distance_ratio < r_build_up:
        #     return (self.section_distance_from_throat - dist_min) / (r_build_up * len_cc)
        # else:
        #     return 1

    @property
    def section_radius(self):
        return self.thrust_chamber.get_radius(self.section_distance_from_throat)

    @property
    def section_hot_gas_convective_heat_transfer_coefficient(self):
        mf = self.combustion_chamber_flow_state.mass_flow
        di = self.section_radius * 2
        mu = self.combustion_chamber_flow_state.dynamic_viscosity
        cp = self.combustion_chamber_flow_state.specific_heat_capacity
        pr = self.combustion_chamber_flow_state.prandtl_number
        t0 = self.section_hot_gas_adiabatic_wall_temp
        tf = self.section_hot_gas_film_temp
        if self.hot_gas_convective_heat_transfer_coefficient_mode == 'Bartz':
            m = self.thrust_chamber.get_mach(self.section_distance_from_throat)
            tw = self.section_hot_side_wall_temp
            y = self.combustion_chamber_flow_state.heat_capacity_ratio
            dt = self.thrust_chamber.get_radius(0) * 2
            rc = self.thrust_chamber.nozzle.div_longi_throat_radius
            p0 = self.combustion_chamber_flow_state.pressure
            ai_at = di ** 2 / dt ** 2
            cstar = 1847.1
            s = (1 + m ** 2 * (y - 1) / 2) ** -.12 / (.5 + .5 * (tw / t0) * (1 + m ** 2 * (y - 1) / 2)) ** .68
            h_hg = (0.026 * dt ** -0.2 * mu ** 0.2 * cp * pr ** -0.6 * (p0 / cstar) ** 0.8
                    * (dt / rc) ** .1 * (dt / di) ** 1.8 * s)
        elif self.hot_gas_convective_heat_transfer_coefficient_mode == 'Bartz2':
            t0 = self.combustion_temp  # Overwrite
            m = self.thrust_chamber.get_mach(self.section_distance_from_throat)
            # TODO: Remove
            m = 3
            tw = self.section_hot_side_wall_temp
            y = self.combustion_chamber_flow_state.heat_capacity_ratio
            dt = self.thrust_chamber.get_radius(0) * 2
            p0 = self.combustion_chamber_flow_state.pressure
            cstar = 1847.1
            s = (1 + m ** 2 * (y - 1) / 2) ** -.12 / (.5 + .5 * (tw / t0) * (1 + m ** 2 * (y - 1) / 2)) ** .68
            h_hg = (0.026 * dt ** -0.2 * mu ** 0.2 * cp * pr ** -0.6 * (p0 / cstar) ** 0.8 * (dt / di) ** 1.8 * s)
        elif self.hot_gas_convective_heat_transfer_coefficient_mode == 'ModifiedBartz':
            h_hg = 0.026 * 1.213 * mf ** .8 * di ** -1.8 * mu ** .2 * cp * pr ** -.6 * (t0 / tf) ** .68
        else:
            raise ValueError('Select proper mode for estimation of the hot gas convective heat transfer coefficient')
        return h_hg / 0.026 * 0.0195

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
        re = self.section_coolant_state.get_reynolds(linear_dimension=D)
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

    @property
    def section_wall_length(self):
        r = self.section_radius
        r2 = self.thrust_chamber.get_radius(self.next_section_distance_from_throat)
        dy = abs(r - r2)
        dx = self.section_axial_length
        return (dx ** 2 + dy ** 2) ** .5

    @property
    def section_pressure_change(self):
        rho = self.section_coolant_state.density
        v = self.section_coolant_state.flow_speed
        l = self.section_wall_length
        fd = self.section_friction_factor
        print(self.section_coolant_bulk_temp, rho, l)
        D = self.coolant_channel_diameter
        dp = fd * l / D * .5 * rho * v ** 2
        return -1 * dp
        # return -1 * self.combustion_chamber_flow_state.pressure * .2 / self.amount_of_sections

    def get_friction_factor(self, roughness_height: float):
        D = self.coolant_channel_diameter
        re = self.section_coolant_state.get_reynolds(linear_dimension=D)
        hs = roughness_height

        def colebrook(darcy_weisbach_friction_factor, rough_height, diameter, reynolds):
            f = darcy_weisbach_friction_factor
            h = rough_height
            d = diameter
            r = reynolds
            return - 2 * log(h / (3.7 * d) + 2.51 / (r * f ** .5), 10) - (1 / (f ** .5))

        answer = scipy.optimize.fsolve(colebrook, 1e-5, args=(hs, D, re))
        return float(answer)

    @property
    def section_friction_factor(self):
        return self.get_friction_factor(roughness_height=self.coolant_channel_roughness_height)

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
        # print(f'{delta_p * -1:.0f} Pa')
        self.section_coolant_total_pressure += delta_p

    # Plot Stuff Below
    @property
    def distances(self):
        return self.data['Distance from Throat [m]']

    def plot_values(self, variable: str, extra_func: Optional[Callable] = None):
        fig, ax = plt.subplots()

        self.plot_nozzle_contour_background(ax=ax)
        data = self.data[variable]
        if extra_func is not None:
            extra_func(ax, data)

        for key, value in data.items():
            ax.plot(self.distances, value, label=key)

        ax.legend()
        ax.set_ylabel(variable)
        ax.set_xlabel('Distance from Throat [m]')
        ax.set_title(variable.split('[')[0].strip(' '))
        plt.show()

    def plot_temps(self):
        self.plot_values('Temperature [K]')

    def plot_wall_temps(self):
        def extra_func(ax: plt.Axes, data: dict):
            keys = list(data)
            for key in keys:
                if 'HotGas' in key:
                    del data[key]

        self.plot_values('Temperature [K]', extra_func=extra_func)

    def plot_coeffs(self):
        def extra_func(ax: plt.Axes, data: dict):
            ax2 = ax.twinx()
            ax2.plot(self.distances, data['ColdSide'], label='ColdSide', linestyle='-.')
            ax2.set_ylabel('Cold Hc')
            ax2.legend()
            del data['ColdSide']
        self.plot_values('Heat-Transfer Coefficient [W/(K*m2]', extra_func=extra_func)

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

    def plot_coolant_properties(self):
        self.plot_values(['Reynolds'])

    def plot_all(self):
        self.plot_coeffs()
        self.plot_temps()
        self.plot_wall_temps()
        self.plot_flux()
        self.plot_coolant()


@dataclass
class DetailedHeatExchanger(HeatExchanger):
    coolant_channel_fin_thickness: float = 1e-3
    coolant_channel_height_input: float = 1e-3
    section_length_to_start_channel: float = field(init=False, repr=False, default=.01)
    coolant_heat_transfer_coefficient_mode: str = 'Taylor'

    def set_next_state(self):
        super().set_next_state()
        self.section_length_to_start_channel += self.section_wall_length

    @cached_property
    def coolant_channel_height(self):
        return self.coolant_channel_height_input

    @property
    def section_coolant_channel_width(self):
        section_circumference = self.section_radius * pi * 2
        return section_circumference / self.number_of_coolant_channels - self.coolant_channel_fin_thickness

    @property
    def coolant_channel_area(self):
        return self.coolant_channel_height * self.section_coolant_channel_width

    @property
    def coolant_channel_diameter(self):
        h = self.coolant_channel_height
        w = self.section_coolant_channel_width
        return 2 * h * w / (h + w)
        # return 2 * (self.coolant_channel_area / pi) ** .5

    @property
    def coolant_channel_equivalent_diameter(self):
        h = self.coolant_channel_height
        w = self.section_coolant_channel_width
        return 2 * h * w / (h + w)

    @property
    def total_coolant_flow_area(self):
        return self.coolant_channel_area * self.number_of_coolant_channels

    @property
    def section_coolant_convective_heat_transfer_coefficient(self):

        k_c = self.section_coolant_state.conductivity
        D = self.coolant_channel_diameter
        re = self.section_coolant_state.get_reynolds(linear_dimension=D)
        pr = self.section_coolant_state.prandtl_number
        t_b = self.section_coolant_bulk_temp
        t_w = self.section_cold_side_wall_temp

        if self.coolant_heat_transfer_coefficient_mode == 'Taylor':
            x = self.section_length_to_start_channel
            a = 0.023
            exponent = 0.57 - 1.59 * D / x
        elif self.coolant_heat_transfer_coefficient_mode == 'SiederTate':
            a = 0.023
            exponent = .55
        elif self.coolant_heat_transfer_coefficient_mode == 'DittusBoelter':
            a = 0.023
            exponent = 0
        elif self.coolant_heat_transfer_coefficient_mode == 'HessKunz':
            raise NotImplementedError

        h_c_base = a * k_c / D * re ** .8 * pr ** .4 * (t_b / t_w) ** exponent
        self.hc_base = h_c_base
        k_r = self.coolant_coefficient_roughness_correction
        h_c = h_c_base * k_r

        # Fin Correction: Luka Denies - Regenerative cooling analysis of oxygen/methane rocket engines 2015
        th_fin = self.coolant_channel_fin_thickness
        ht_c = self.coolant_channel_height
        w_c = self.section_coolant_channel_width
        k = self.chamber_wall_conductivity
        b = (2 * h_c * th_fin / k) ** .5 / th_fin * ht_c
        eta = tanh(b) / b
        c_fin = (w_c + eta * 2 * ht_c) / (w_c + th_fin)
        self.coolant_coefficient_fin_correction = c_fin
        h_c = h_c * c_fin
        return h_c
    
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
                                       'T_b': []},
                     'Geometry [mm]': {'H':[],
                                  # 'W':[],
                                  'D':[],}
                     }

    def write_data(self):
        super().write_data()
        self.data['Geometry [mm]']['H'].append(self.coolant_channel_height*1e3)
        # self.data['Geometry [mm]']['W'].append(self.section_coolant_channel_width*1e3)
        self.data['Geometry [mm]']['D'].append(self.coolant_channel_diameter*1e3)
    
    def plot_geometry(self):
        self.plot_values('Geometry [mm]')

    @property
    def coolant_coefficient_roughness_correction(self):
        k = self.section_friction_factor / self.get_friction_factor(roughness_height=0)
        pr = self.section_coolant_state.prandtl_number
        re = self.section_coolant_state.get_reynolds(linear_dimension=self.coolant_channel_diameter)
        b = 1.5 * pr ** (-1. / 6.) * re ** (-1. / 8.)
        c_rough = ((1 + b * (pr - 1)) / (1 + b * (pr * k - 1)) * k)
        return c_rough
