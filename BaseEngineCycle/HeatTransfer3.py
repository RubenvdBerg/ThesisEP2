import numpy as np
from BaseEngineCycle.ThrustChamber import ThrustChamber
from dataclasses import dataclass, replace, field
from BaseEngineCycle.FlowState import FlowState, ManualFlowState
from typing import Optional
from math import pi
from CoolProp import CoolProp
import warnings


@dataclass
class HotSideHeatTransferSection:
    wall_temperature: float  # [K]
    distance_from_throat: float  # [m]
    radiative_factor: float
    thrust_chamber: ThrustChamber
    combustion_chamber_flow_state: ManualFlowState

    post_injection_build_up_ratio: Optional[float] = None  # [-]
    recovery_factor: Optional[float] = None  # [-]
    verbose: bool = True

    def __post_init__(self):
        if self.recovery_factor is None:
            self.recovery_factor = self.turbulent_recovery_factor
        if self.post_injection_build_up_ratio is None:
            self.post_injection_build_up_ratio = 0.25

    @property
    def min_distance_from_throat(self):
        return self.thrust_chamber.min_distance_from_throat

    @property
    def max_distance_from_throat(self):
        return self.thrust_chamber.max_distance_from_throat

    @property
    def turbulent_recovery_factor(self):
        # Zandbergen 2017 p.160
        return self.combustion_chamber_flow_state.prandtl_number ** (1 / 3)

    @property
    def radius(self):
        return self.thrust_chamber.get_radius(self.distance_from_throat)

    @property
    def mach(self):
        return self.thrust_chamber.get_mach(self.distance_from_throat)

    @property
    def heat_capacity_ratio(self):
        return self.combustion_chamber_flow_state.heat_capacity_ratio

    @property
    def combustion_temperature(self):
        return self.combustion_chamber_flow_state.temperature

    @property
    def static_temp(self):
        return self.reduction_factor * self.combustion_temperature / (1 + (self.heat_capacity_ratio - 1)
                                                                      / 2 * self.mach ** 2)

    @property
    def adiabatic_wall_temp(self):
        return self.static_temp * (1 + (self.heat_capacity_ratio - 1) / 2 * self.mach ** 2 * self.recovery_factor)

    @property
    def film_temperature(self):
        return self.static_temp*.28 + self.wall_temperature*.5 + self.adiabatic_wall_temp*.22

    @property
    def convective_heat_transfer_coefficient(self):
        mf = self.combustion_chamber_flow_state.mass_flow
        di = self.radius * 2
        mu = self.combustion_chamber_flow_state.dynamic_viscosity
        cp = self.combustion_chamber_flow_state.specific_heat_capacity
        pr = self.combustion_chamber_flow_state.prandtl_number
        t0 = self.static_temp
        tf = self.film_temperature
        return 0.026 * 1.213 * mf ** .8 * di ** -1.8 * mu ** .2 * cp * pr ** -.6 * (t0 / tf) ** .68

    @property
    def reduction_factor(self):
        """Factor that reduces the heat flow until maximum heat flux build up after the injector is achieved"""
        len_cc = self.thrust_chamber.chamber.length
        dist_min = self.thrust_chamber.min_distance_from_throat
        r_build_up = self.post_injection_build_up_ratio
        # Distance from injector divided by total chamber length
        inj_distance_ratio = (self.distance_from_throat - dist_min) / len_cc
        if inj_distance_ratio < r_build_up:
            return (self.distance_from_throat - dist_min) / (r_build_up * len_cc)
        else:
            return 1


    @property
    def convective_heat_flux(self):
        return self.convective_heat_transfer_coefficient * (self.adiabatic_wall_temp - self.wall_temperature)

    @property
    def radiative_heat_flux(self):
        return self.convective_heat_flux * self.radiative_factor

    @property
    def total_heat_transfer_coefficient(self):
        return self.convective_heat_transfer_coefficient * (1 + self.radiative_factor)

    @property
    def total_heat_flux(self):
        """Should be the same as check in comments"""
        check = self.total_heat_transfer_coefficient * (self.adiabatic_wall_temp - self.wall_temperature)
        value = self.convective_heat_flux + self.radiative_heat_flux
        if not check == value:
            warnings.warn(f'Total Heat Flux [W/m2] -> Value:{value:.4e} Check:{check:.4e}')
        return value

    @property
    def heat_transfer_per_axial_meter(self):
        return 2 * pi * self.radius * self.total_heat_flux


@dataclass
class CoolantSideHeatTransferSection:
    wall_temperature: float  # [K]
    coolant_flow_state: FlowState
    channel_diameter: float  # [m]
    number_of_channels: float  # [-]

    @property
    def bulk_temperature(self):
        return self.coolant_flow_state.temperature

    @property
    def channel_area(self):
        return (pi / 4) * self.channel_diameter ** 2

    @property
    def mass_flow_per_channel(self):
        return self.coolant_flow_state.mass_flow / self.number_of_channels

    @property
    def channel_mass_flux(self):
        return self.mass_flow_per_channel / self.channel_area

    @property
    def channel_flow_speed(self):
        return self.channel_mass_flux / self.coolant_flow_state.density

    @property
    def reynolds(self):
        return self.coolant_flow_state.get_reynolds(flow_speed=self.channel_flow_speed,
                                                    linear_dimension=self.channel_diameter)

    @property
    def convective_heat_transfer_coefficient(self):
        k_c = self.coolant_flow_state.conductivity
        D = self.channel_diameter
        re = self.reynolds
        pr = self.coolant_flow_state.prandtl_number
        t_b = self.bulk_temperature
        t_w = self.wall_temperature
        # print(f'kc:{k_c}, D:{D}, re:{re}, pr:{pr}, tb:{t_b}, tw:{t_w}')
        return 0.025 * k_c / D * re ** .8 * pr ** .4 * (t_b / t_w) ** .55

    @property
    def heat_flux(self):
        return self.convective_heat_transfer_coefficient * (self.wall_temperature - self.bulk_temperature)


@dataclass
class HeatTransferSection:
    coolant_inlet_flow_state: FlowState
    coolant_channel_diameter: float
    number_of_coolant_channels: float
    radiative_factor: float
    amount_of_sections: float
    thrust_chamber: ThrustChamber
    combustion_chamber_flow_state: ManualFlowState
    data: dict = field(init=False)

    wall_temperature: Optional[float] = None
    distance_from_throat: Optional[float] = None
    bulk_temperature: Optional[float] = None

    post_injection_build_up_ratio: Optional[float] = None  # [-]
    recovery_factor: Optional[float] = None  # [-]
    verbose: bool = True
    iteration_accuracy: float = 0.0001

    def __post_init__(self):
        if self.wall_temperature is None:
            self.wall_temperature = self.coolant_inlet_flow_state.temperature
        if self.distance_from_throat is None:
            self.distance_from_throat = self.thrust_chamber.min_distance_from_throat + 1e-9
        if self.bulk_temperature is None:
            self.bulk_temperature = self.coolant_inlet_flow_state.temperature
        self.init_data()
        self.iterate()

    def iterate_section(self):
        while abs(
                self.wall_temperature - self.expected_wall_temperature) / self.expected_wall_temperature > self.iteration_accuracy:
            if self.verbose:
                print(
                    f'Wall Temperature [K]-> Current:{self.wall_temperature:.1f}, Expected:{self.expected_wall_temperature:.1f}')
            self.wall_temperature = self.expected_wall_temperature

    def iterate(self):
        for x in np.linspace(self.thrust_chamber.max_distance_from_throat, self.thrust_chamber.min_distance_from_throat,
                             self.amount_of_sections):
            self.distance_from_throat = x
            self.iterate_section()
            self.bulk_temperature = self.ideal_outlet_temperature
            self.write_data()
        if self.verbose:
            print('Heat Transfer Iteration Complete')

    def init_data(self):
        self.data = {'Temperature [K]': [],
                     'Heat-Transfer Coefficient [W/(K*m2]': [],
                     'Heat Flux [W/m2]': [],
                     'Distance From Throat [m]': []}

    def write_data(self):
        self.data['Temperature [K]'].append([self.bulk_temperature,
                                              self.wall_temperature,
                                              self.hot_side_heat_transfer_section.adiabatic_wall_temp,
                                              self.hot_side_heat_transfer_section.static_temp,
                                              self.hot_side_heat_transfer_section.film_temperature,])
        self.data['Heat-Transfer Coefficient [W/(K*m2]'].append([self.hot_side_heat_transfer_section.convective_heat_transfer_coefficient,
                                              self.hot_side_heat_transfer_section.total_heat_transfer_coefficient,
                                              self.coolant_side_heat_transfer_section.convective_heat_transfer_coefficient, ])
        self.data['Heat Flux [W/m2]'].append([self.heat_flux,
                                                self.coolant_side_heat_transfer_section.heat_flux,
                                                self.hot_side_heat_transfer_section.total_heat_flux,
                                                self.hot_side_heat_transfer_section.convective_heat_flux,
                                                self.hot_side_heat_transfer_section.radiative_heat_flux, ])
        self.data['Distance From Throat [m]'].append([self.distance_from_throat])

    @property
    def section_length(self):
        return (self.thrust_chamber.max_distance_from_throat
                - self.thrust_chamber.min_distance_from_throat) / (self.amount_of_sections - 1)

    @property
    def coolant_flow_state(self):
        return replace(self.coolant_inlet_flow_state,
                       temperature=self.bulk_temperature, )

    @property
    def hot_side_heat_transfer_section(self):
        return HotSideHeatTransferSection(wall_temperature=self.wall_temperature,
                                          distance_from_throat=self.distance_from_throat,
                                          radiative_factor=self.radiative_factor,
                                          thrust_chamber=self.thrust_chamber,
                                          combustion_chamber_flow_state=self.combustion_chamber_flow_state,
                                          post_injection_build_up_ratio=self.post_injection_build_up_ratio,
                                          recovery_factor=self.recovery_factor,
                                          verbose=self.verbose)

    @property
    def coolant_side_heat_transfer_section(self):
        return CoolantSideHeatTransferSection(wall_temperature=self.wall_temperature,
                                              coolant_flow_state=self.coolant_flow_state,
                                              channel_diameter=self.coolant_channel_diameter,
                                              number_of_channels=self.number_of_coolant_channels,
                                              )

    @property
    def heat_flux(self):
        h_g_eff = self.hot_side_heat_transfer_section.total_heat_transfer_coefficient
        h_c_eff = self.coolant_side_heat_transfer_section.convective_heat_transfer_coefficient
        t_w_ad = self.hot_side_heat_transfer_section.adiabatic_wall_temp
        t_b = self.coolant_side_heat_transfer_section.bulk_temperature
        return (t_w_ad - t_b) / (1 / h_g_eff + 1 / h_c_eff)

    @property
    def heat_transfer(self):
        return 2 * pi * self.thrust_chamber.get_radius(self.distance_from_throat) * self.heat_flux * self.section_length

    @property
    def hot_side_heat_transfer(self):
        return self.hot_side_heat_transfer_section.heat_transfer_per_axial_meter * self.section_length

    @property
    def estimated_wal_temperature2(self):
        return (self.hot_side_heat_transfer_section.adiabatic_wall_temp
                - self.heat_flux / self.hot_side_heat_transfer_section.total_heat_transfer_coefficient)

    @property
    def increase_mass_specific_enthalpy(self):
        return self.heat_transfer / self.coolant_inlet_flow_state.mass_flow

    @property
    def outlet_mass_specific_enthalpy(self):
        return self.coolant_flow_state.mass_specific_enthalpy + self.increase_mass_specific_enthalpy

    @property
    def ideal_outlet_temperature(self):
        return CoolProp.PropsSI('T', 'H', self.outlet_mass_specific_enthalpy, 'P',
                                self.coolant_inlet_flow_state.pressure,
                                self.coolant_inlet_flow_state.coolprop_name)

    @property
    def expected_bulk_temperature(self):
        return (self.coolant_inlet_flow_state.temperature + self.ideal_outlet_temperature) / 2

    @property
    def expected_wall_temperature(self):
        h_g_eff = self.hot_side_heat_transfer_section.total_heat_transfer_coefficient
        h_c_eff = self.coolant_side_heat_transfer_section.convective_heat_transfer_coefficient
        t_w_ad = self.hot_side_heat_transfer_section.adiabatic_wall_temp
        t_b = self.coolant_side_heat_transfer_section.bulk_temperature
        return (h_g_eff * t_w_ad + h_c_eff * t_b) / (h_g_eff + h_c_eff)

    @property
    def distances(self):
        return [x for x in zip(*self.data['Distance From Throat [m]'])][0]

    def plot_values(self, variable: str, line_names: list[str,...]):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        for temp_data, name in zip(zip(*self.data[variable]),line_names):
            ax.plot(self.distances, temp_data, label=name)
        ax2 = ax.twinx()
        ax2.plot(self.distances, [self.thrust_chamber.get_radius(x) for x in self.distances], color='grey', linestyle='--')
        ax2.get_yaxis().set_visible(False)
        ax.legend()
        ax.set_ylabel(variable)
        ax.set_xlabel('Distance from Throat [m]')
        ax.set_title(variable.split('[')[0].strip(' '))
        plt.show()

    def plot_temps(self):
        self.plot_values('Temperature [K]', ['Bulk', 'Wall', 'Adiabatic', 'Static', 'Film'])

    def plot_coeffs(self):
        self.plot_values('Heat-Transfer Coefficient [W/(K*m2]', ['Hot-Gas Conv.', 'Hot-Gas Total', 'Coolant'])

    def plot_flux(self):
        self.plot_values('Heat Flux [W/m2]', ['q', 'q_c', 'q_hg', 'q_hg_c', 'q_hg_r'])

    def plot_all(self):
        self.plot_coeffs()
        self.plot_temps()
        self.plot_flux()

@dataclass
class HeatTransferPlotter:
    data: dict

    @property
    def distances(self):
        return [x for x in zip(*self.data['Distance From Throat [m]'])][0]

    def plot_temps(self):
        import matplotlib.pyplot as plt
        for temp_data, name in zip(zip(*self.data['Temperature [K]']),['Bulk', 'Wall', 'Adiabatic', 'Static']):
            fig, ax = plt.subplots()
            ax.plot(self.distances, temp_data, label=name)
            ax2 = ax.secondary_yaxis(location='right')
            ax2.plot()
        plt.legend()
        plt.show()



if __name__ == '__main__':
    from BaseEngineCycle.EngineCycle import EngineCycle
    from BaseEngineCycle.Nozzle import ConicalNozzle
    from BaseEngineCycle.Injector import Injector
    from BaseEngineCycle.CombustionChamber import CombustionChamber
    import arguments as args
    from math import radians

    main_kwargs = args.change_to_conical_nozzle(args.tcd1_kwargs, throat_half_angle=radians(25))
    engine = EngineCycle(**main_kwargs, **args.duel_pump_kwargs)
    injector = Injector(material_density=0, safety_factor=0, yield_strength=0)
    nozzle = ConicalNozzle(throat_area=engine.throat_area,
                           expansion_ratio=main_kwargs['expansion_ratio'],
                           conv_chamber_bend_ratio=main_kwargs['convergent_chamber_bend_ratio'],
                           conv_throat_bend_ratio=main_kwargs['convergent_throat_bend_ratio'],
                           conv_half_angle=main_kwargs['convergent_half_angle'],
                           div_throat_half_angle=main_kwargs['divergent_throat_half_angle'],
                           area_ratio_chamber_throat=main_kwargs['area_ratio_chamber_throat'])
    chamber = CombustionChamber(material_density=0, safety_factor=0, yield_strength=0,
                                throat_area=engine.throat_area,
                                combustion_chamber_pressure=engine.combustion_chamber_pressure,
                                convergent_volume_estimate=nozzle.conv_volume_estimate,
                                area_ratio_chamber_throat=main_kwargs['area_ratio_chamber_throat'],
                                propellant_mix=engine.propellant_mix_name,
                                characteristic_length=main_kwargs['chamber_characteristic_length'],
                                verbose=True)
    manual_cc_flow_state = ManualFlowState(propellant_name='Combustion_Chamber_Gas',
                                           temperature=engine.combustion_temperature,
                                           pressure=engine.combustion_chamber_pressure,
                                           mass_flow=engine.chamber_mass_flow,
                                           type='combusted',
                                           _specific_heat_capacity=engine.cc_hot_gas_specific_heat_capacity,
                                           _heat_capacity_ratio=engine.cc_hot_gas_heat_capacity_ratio,
                                           _prandtl_number=engine.cc_hot_gas_prandtl_number,
                                           _dynamic_viscosity=engine.cc_hot_gas_dynamic_viscosity,
                                           )

    thrustchamber = ThrustChamber(injector=injector, chamber=chamber, nozzle=nozzle,
                                  heat_capacity_ratio=manual_cc_flow_state.heat_capacity_ratio)

    heattransfer = HeatTransferSection(distance_from_throat=thrustchamber.min_distance_from_throat + 0.0001,
                                       coolant_inlet_flow_state=engine.cooling_inlet_flow_state,
                                       coolant_channel_diameter=.01,
                                       number_of_coolant_channels=150,
                                       radiative_factor=engine.radiative_heat_transfer.radiative_factor,
                                       thrust_chamber=thrustchamber,
                                       combustion_chamber_flow_state=manual_cc_flow_state,
                                       amount_of_sections=300,
                                       verbose=False
                                       )

    print(heattransfer.data)
    heattransfer.plot_all()