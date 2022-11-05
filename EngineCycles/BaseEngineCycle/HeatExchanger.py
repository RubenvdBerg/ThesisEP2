from dataclasses import dataclass, field
from typing import Optional, Callable
from functools import cached_property
from numpy import linspace
from scipy import integrate
from math import pi
from EngineCycles.BaseEngineCycle.ThrustChamber import ThrustChamberSection
from EngineCycles.BaseEngineCycle.FlowState import FlowState, ManualFlowState, CoolantFlowState
import EngineCycles.Functions.EmpiricalRelations as empirical
from CoolProp import CoolProp

from matplotlib import pyplot as plt
from EngineCycles.Functions.BaseFunctions import multi_legend


@dataclass
class HeatExchanger:
    thrust_chamber_section: ThrustChamberSection
    combustion_chamber_flow_state: ManualFlowState
    coolant_inlet_flow_state: FlowState
    number_of_coolant_channels: float
    radiative_factor: float
    maximum_wall_temp: float
    is_counter_flow: bool = False
    max_iterations: float = 10
    amount_of_sections: float = 500
    post_injection_build_up_ratio: float = .25  # [-]
    coolant_channel_roughness_height: float = 6e-6

    _coolant_channel_diameter: Optional[float] = None
    _recovery_factor: Optional[float] = None  # [-]
    _initial_flow_speed: float = 10  # [m/s]
    _save_data: bool = False
    verbose: bool = True
    iteration_accuracy: float = 1e-3
    data: dict = field(init=False)

    # Iterative section variables, not required at init
    section_distance_from_throat: float = field(init=False, repr=False)
    next_section_distance_from_throat: float = field(init=False, repr=False)
    section_coolant_total_temp: float = field(init=False, repr=False)
    section_coolant_total_pressure: float = field(init=False, repr=False)
    section_coolant_state: CoolantFlowState = field(init=False, repr=False)
    section_heat_flux: float = field(init=False, repr=False)
    section_local_mach: float = field(init=False, repr=False)
    section_hot_gas_adiabatic_wall_temp: float = field(init=False, repr=False)
    section_hot_gas_static_temp: float = field(init=False, repr=False)
    section_static_guess_tuple: tuple = field(init=False, repr=False, default=(None, None))

    def __post_init__(self):
        # Set values at start of coolant channel
        self.section_coolant_total_temp = self.coolant_inlet_flow_state.temperature
        self.section_coolant_total_pressure = self.coolant_inlet_flow_state.pressure
        if self._save_data:
            self.init_data()
        self.iterate()
        self.output_message()

    # Cached properties below remain constant throughout iteration
    @cached_property
    def axial_length_one_section(self):
        return self.thrust_chamber_section.length / self.amount_of_sections

    @cached_property
    def mass_flow_per_coolant_channel(self):
        return self.coolant_inlet_flow_state.mass_flow / self.number_of_coolant_channels

    @cached_property
    def combustion_temp(self):
        return self.combustion_chamber_flow_state.temperature

    @cached_property
    def hot_wall_temp(self):
        """Pass along maximum wall temp, method required for child overwrite."""
        return self.maximum_wall_temp

    @cached_property
    def recovery_factor(self):
        if self._recovery_factor is None:
            # Zandbergen 2017 p.160, experimental relationship for turbulent boundary layers
            return self.combustion_chamber_flow_state.prandtl_number ** (1 / 3)
        else:
            return self._recovery_factor

    @cached_property
    def axial_distances_inputs(self):
        """Return distance from throat to start of each section, pairwise.

        Normally process starts at injector face, for counter flow it start at the nozzle exit instead.
        """
        if self.is_counter_flow:
            linspace_input = reversed(self.thrust_chamber_section.throat_distance_tuple)
        else:
            linspace_input = self.thrust_chamber_section.throat_distance_tuple
        axial_distances = linspace(*linspace_input, int(self.amount_of_sections))
        return list(zip(axial_distances, axial_distances[1:]))

    def iterate(self):
        """Calculate (and store) the heat transfer variables at every section."""
        i_tot = len(self.axial_distances_inputs)
        for i, (x, x_next) in enumerate(self.axial_distances_inputs):
            if self.verbose:
                print(f'{i / i_tot * 100:.1f}% done')

            self.section_distance_from_throat = x
            self.next_section_distance_from_throat = x_next
            self.iterate_coolant_dynamic_state()
            self.determine_heat_flux()

            if self._save_data:
                self.write_data()

            self.set_next_state()

        if self.verbose:
            print('Heat Transfer Iteration Complete')

    def iterate_coolant_dynamic_state(self):
        """Iterate dynamic state, which happens internally in CoolantFlowState."""
        self.section_coolant_state = CoolantFlowState(propellant_name=self.coolant_inlet_flow_state.propellant_name,
                                                      total_temperature=self.section_coolant_total_temp,
                                                      total_pressure=self.section_coolant_total_pressure,
                                                      mass_flow=self.coolant_inlet_flow_state.mass_flow,
                                                      type=self.coolant_inlet_flow_state.type,
                                                      total_flow_area=self.section_total_coolant_flow_area,
                                                      _iteration_accuracy=self.iteration_accuracy,
                                                      verbose=self.verbose,
                                                      _static_temperature_initial_guess=self.section_static_guess_tuple[0],
                                                      _static_pressure_initial_guess=self.section_static_guess_tuple[1],)

    def determine_heat_flux(self):
        self.set_local_mach()
        self.set_hot_gas_temps()
        hg = self.section_hot_gas_total_heat_transfer_coefficient
        tw = self.hot_wall_temp
        tr = self.section_hot_gas_adiabatic_wall_temp
        self.section_heat_flux = hg * (tr - tw)

    def set_next_state(self):
        """Calculate change in temperature and pressure and update."""
        t_guess = self.section_coolant_state._static_temperature
        p_guess = self.section_coolant_state._static_pressure
        self.section_static_guess_tuple = (t_guess, p_guess)
        self.section_coolant_total_temp += self.section_temperature_change
        self.section_coolant_total_pressure += self.section_pressure_change

    def output_message(self):
        p0_1 = self.coolant_inlet_flow_state.pressure
        p0_2 = self.section_coolant_total_pressure
        T0_1 = self.coolant_inlet_flow_state.temperature
        T0_2 = self.section_coolant_total_temp
        print(f'Coolant Outlet Total Temp.   : {T0_2:.2f} K')
        print(f'Coolant Outlet Total Pressure: {p0_2*1e-5:.3f} bar')
        print(f'Coolant Change in Temp.      : {T0_2 - T0_1:.2f} K')
        print(f'Coolant Pressure Drop        : {(p0_1 - p0_2) * 1e-5:.3f} bar')
        if self._save_data:
            print(
                f'DeltaP        : {(self.data["Coolant State"]["p"][-1] - self.data["Coolant State"]["p"][0]) * 1e-5 :.2f} bar')
            print(
                f'DeltaP0       : {(self.data["Coolant State"]["p0"][-1] - self.data["Coolant State"]["p0"][0]) * 1e-5 :.2f} bar')

    def error_is_small(self, current: float, expected: float):
        error = abs((current - expected) / expected)
        return error < self.iteration_accuracy

    @property
    def section_coolant_channel_diameter(self):
        if self._coolant_channel_diameter is None:
            # Choose channel diameter based on an assumed initial flow speed and number of channels
            chosen_flow_speed = self._initial_flow_speed
            initial_density = self.coolant_inlet_flow_state.density
            mass_flow = self.coolant_inlet_flow_state.mass_flow
            volume_flow = mass_flow / initial_density
            required_area = volume_flow / chosen_flow_speed
            channel_area = required_area / self.number_of_coolant_channels
            return 2 * (channel_area / pi) ** .5
        else:
            return self._coolant_channel_diameter

    @property
    def section_coolant_channel_area(self):
        return (pi / 4) * self.section_coolant_channel_diameter ** 2

    @property
    def section_total_coolant_flow_area(self):
        return self.section_coolant_channel_area * self.number_of_coolant_channels

    @property
    def section_hot_gas_film_temp(self):
        return empirical.get_film_temp(wall_temp=self.hot_wall_temp,
                                       adiabatic_wall_temp=self.section_hot_gas_adiabatic_wall_temp,
                                       static_temp=self.section_hot_gas_static_temp)

    @property
    def section_coolant_bulk_temp(self):
        return self.section_coolant_state.static_temperature

    @property
    def section_reduction_factor(self):
        """Reduces heat transfer variables at the injector face to 0 and builds up to 1"""
        len_cc = self.thrust_chamber_section.chamber.length
        dist_min = self.thrust_chamber_section.min_distance_from_throat
        r_build_up = self.post_injection_build_up_ratio

        # Distance from injector divided by total chamber length
        inj_distance_ratio = (self.section_distance_from_throat - dist_min) / len_cc
        if inj_distance_ratio < r_build_up:
            return (self.section_distance_from_throat - dist_min) / (r_build_up * len_cc)
        else:
            return 1

    def set_local_mach(self):
        self.section_local_mach = self.thrust_chamber_section.get_mach(self.section_distance_from_throat)

    def set_hot_gas_temps(self):
        y = self.combustion_chamber_flow_state.heat_capacity_ratio
        m = self.section_local_mach
        t_c = self.combustion_temp
        f_r = self.recovery_factor
        t_hg_st = t_c / (1 + (y - 1) / 2 * m ** 2)
        t_hg_ad = t_hg_st * (1 + f_r * (y - 1) / 2 * m ** 2)
        self.section_hot_gas_static_temp = t_hg_st
        self.section_hot_gas_adiabatic_wall_temp = t_hg_ad

    @property
    def section_hot_gas_convective_heat_transfer_coefficient(self):
        h_g = empirical.get_hot_gas_convective_heat_transfer_coefficient(
            mass_flow=self.coolant_inlet_flow_state.mass_flow,
            local_diameter=self.section_radius * 2,
            dynamic_viscosity=self.combustion_chamber_flow_state.dynamic_viscosity,
            specific_heat_capacity=self.combustion_chamber_flow_state.specific_heat_capacity,
            prandtl_number=self.combustion_chamber_flow_state.prandtl_number,
            stagnation_temp=self.section_hot_gas_adiabatic_wall_temp,
            film_temp=self.section_hot_gas_film_temp,
            mode='ModifiedBartz',
        )
        return h_g * self.section_reduction_factor

    @property
    def section_hot_gas_total_heat_transfer_coefficient(self):
        return self.section_hot_gas_convective_heat_transfer_coefficient * (1 + self.radiative_factor)

    @property
    def section_radius(self):
        return self.thrust_chamber_section.get_radius(self.section_distance_from_throat)

    @property
    def section_wall_length(self):
        r = self.section_radius
        r2 = self.thrust_chamber_section.get_radius(self.next_section_distance_from_throat)
        dy = abs(r - r2)
        dx = self.axial_length_one_section
        return (dx ** 2 + dy ** 2) ** .5

    @property
    def section_wall_surface(self):
        # noinspection PyTupleAssignmentBalance
        section_surface, _ = integrate.quad(func=lambda x: 2 * pi * self.thrust_chamber_section.get_radius(x),
                                            a=self.section_distance_from_throat,
                                            b=self.next_section_distance_from_throat)
        return abs(section_surface)

    @property
    def section_temperature_change(self):
        A = self.section_wall_surface
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
        return delta_t

    @property
    def section_friction_factor(self):
        e = self.coolant_channel_roughness_height
        Dh = self.section_coolant_channel_diameter
        re = self.section_coolant_state.get_reynolds(linear_dimension=Dh)
        return empirical.get_friction_factor(roughness_height=e,
                                             diameter=Dh,
                                             reynolds_number=re)

    @property
    def section_pressure_change(self):
        rho = self.section_coolant_state.density
        v = self.section_coolant_state.flow_speed
        l = self.section_wall_length
        fd = self.section_friction_factor
        Dh = self.section_coolant_channel_diameter
        dp = fd * l / Dh * .5 * rho * v ** 2
        return -1 * dp

    # Properties below are only for writing to data, not required for any other calculations
    @property
    def section_distance_from_injector(self):
        return self.section_distance_from_throat - self.thrust_chamber_section.min_distance_from_throat

    @property
    def section_hot_gas_convective_heat_flux(self):
        return self.section_hot_gas_convective_heat_transfer_coefficient * (self.section_hot_gas_adiabatic_wall_temp
                                                                            - self.hot_wall_temp)

    @property
    def section_hot_gas_radiative_heat_flux(self):
        return self.section_hot_gas_convective_heat_flux * self.radiative_factor

    @property
    def section_hot_gas_total_heat_flux(self):
        return self.section_hot_gas_convective_heat_flux + self.section_hot_gas_radiative_heat_flux

    def init_data(self):
        self.data = {'Distance from Throat [m]': [],
                     'Temperature [K]': {'CoolantBulk': [],
                                         'HotGasAdiabatic': [],
                                         'HotGasStatic': [],
                                         'HotGasFilm': [], },
                     'Heat-Transfer Coefficient [W/(K*m2]': {'Hot Gas': [],
                                                             'Hot SideConv.': [], },
                     'Heat Flux [W/m2]': {'Total': [],
                                          'Hot SideConv.': [],
                                          'Hot SideRad.': [], },
                     'Coolant State': {'T': [],
                                       'T0': [],
                                       'p': [],
                                       'p0': [],
                                       'rho': [],
                                       'cp': [], },
                     'Channel Geometry [mm]': {'Diameter': []},
                     'Plot Vals [cm]': {'x': [],
                                        'y': [], }
                     }

    def write_data(self):
        self.data['Distance from Throat [m]'].append(self.section_distance_from_throat)
        self.data['Temperature [K]']['CoolantBulk'].append(self.section_coolant_bulk_temp)
        self.data['Temperature [K]']['HotGasAdiabatic'].append(self.section_hot_gas_adiabatic_wall_temp)
        self.data['Temperature [K]']['HotGasStatic'].append(self.section_hot_gas_static_temp)
        self.data['Temperature [K]']['HotGasFilm'].append(self.section_hot_gas_film_temp)
        self.data['Heat-Transfer Coefficient [W/(K*m2]']['Hot Gas'].append(
            self.section_hot_gas_total_heat_transfer_coefficient)
        self.data['Heat-Transfer Coefficient [W/(K*m2]']['Hot SideConv.'].append(
            self.section_hot_gas_convective_heat_transfer_coefficient)
        self.data['Heat Flux [W/m2]']['Total'].append(self.section_heat_flux)
        self.data['Heat Flux [W/m2]']['Hot SideConv.'].append(self.section_hot_gas_convective_heat_flux)
        self.data['Heat Flux [W/m2]']['Hot SideRad.'].append(self.section_hot_gas_radiative_heat_flux)
        self.data['Coolant State']['T'].append(self.section_coolant_state.static_temperature)
        self.data['Coolant State']['T0'].append(self.section_coolant_total_temp)
        self.data['Coolant State']['p'].append(self.section_coolant_state.static_pressure)
        self.data['Coolant State']['p0'].append(self.section_coolant_total_pressure)
        self.data['Channel Geometry [mm]']['Diameter'].append(self.section_coolant_channel_diameter * 1e3)
        self.data['Plot Vals [cm]']['x'].append(self.section_distance_from_injector * 1e2)
        self.data['Plot Vals [cm]']['y'].append(self.section_radius * 1e2)


@dataclass
class HeatTransferPlots:
    heattransferdata: dict

    @property
    def distances(self):
        return self.heattransferdata['Distance from Throat [m]']

    def plot_values(self, variable: str, extra_func: Optional[Callable] = None):
        fig, ax = plt.subplots()

        self.plot_nozzle_contour_background(ax=ax)
        data = self.heattransferdata[variable]

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

    def plot_flux(self, **kwargs):
        self.plot_values('Heat Flux [W/m2]', **kwargs)

    def plot_geometry(self):
        self.plot_values('Channel Geometry [mm]')

    def plot_wall_temps(self):
        fig, ax = plt.subplots()

        self.plot_nozzle_contour_background(ax=ax)
        hot_tws = self.heattransferdata['Temperature [K]']['Hot SideWall']
        cold_tws = self.heattransferdata['Temperature [K]']['Cold SideWall']
        heat_fluxs = [q * 1e-6 for q in self.heattransferdata['Heat Flux [W/m2]']['Total']]

        ax2 = ax.twinx()
        ax2.set_ylabel(r'Heat Flux [$MW/m^2$]')

        ax2.plot(self.distances, heat_fluxs, color='darkorange', label='Heat Flux', linestyle='-.')
        ax.plot(self.distances, hot_tws, color='r', label='Hot Side')
        ax.plot(self.distances, cold_tws, color='b', label='Cold Side')

        multi_legend([ax, ax2])
        ax.set_ylim([min(cold_tws) * .7, max(hot_tws) * 1.2])
        ax.set_ylabel(r'Temperature [$K$]')
        ax.set_xlabel(r'Distance from Throat [$m$]')
        ax.set_title('Wall Temperatures')
        plt.show()

    def plot_coeffs(self):
        fig, ax = plt.subplots()

        self.plot_nozzle_contour_background(ax=ax)
        variable_name = 'Heat-Transfer Coefficient [W/(K*m2]'
        ax2 = ax.twinx()
        axes = [ax, ax2]
        for name, color, axis in zip(['Hot Gas', 'Coolant'], ['r', 'b'], axes):
            values = [hc * 1e-3 for hc in self.heattransferdata[variable_name][name]]
            axis.plot(self.distances, values, color=color, label=name)
            axis.set_ylabel(name + r' Heat-Transfer Coefficient [$kW/(K\cdot m^2$)]')

        multi_legend((ax, ax2))
        ax.set_xlabel(r'Distance from Throat [$m$]')
        ax.set_title(variable_name.split('[')[0])
        plt.show()

    def plot_coolant(self):
        c_data = self.heattransferdata['Coolant State']
        fig, ax = plt.subplots()

        self.plot_nozzle_contour_background(ax=ax)

        ax2 = ax.twinx()
        for variable, axis, color, style in zip(['T', 'T0', 'p', 'p0'], [ax, ax, ax2, ax2],
                                                ['orange', 'red', 'lightblue', 'blue'], ['-', '-.', '-', '-.']):
            axis.plot(self.distances, c_data[variable], label=variable, color=color, linestyle=style)

        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc=0)
        ax.set_ylabel(r'Temperature [$K$]')
        ax2.set_ylabel(r'Pressure [$Pa$]')
        ax.set_xlabel(r'Distance from Throat [$m$]')
        ax.set_title('Coolant Bulk State')
        plt.show()

    def plot_nozzle_contour_background(self, ax: plt.Axes):
        ax3 = ax.twinx()
        radii = self.heattransferdata['Plot Vals [cm]']['y']
        ax3.plot(self.distances, radii, color='grey',
                 linestyle='--')
        ax3.set_ylim([min(radii) * .9, max(radii) * 1.4])
        ax3.get_yaxis().set_visible(False)

    def plot_all(self):
        self.plot_coeffs()
        self.plot_temps()
        self.plot_wall_temps()
        self.plot_flux()
        self.plot_coolant()
        self.plot_geometry()
