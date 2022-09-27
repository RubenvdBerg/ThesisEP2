import warnings
from dataclasses import dataclass, field
from functools import cached_property
from math import sqrt, log
from typing import Optional, Literal

from scipy import constants as constants

from BaseEngineCycle.BaseFunctions import get_propellant_mix_name
from BaseEngineCycle.CombustionChamber import CombustionChamber
from BaseEngineCycle.Injector import Injector
from BaseEngineCycle.Cooling import CoolingChannelSection
from BaseEngineCycle.HeatTransferSection import HeatTransferSection, ConvectiveHeatTransfer, RadiativeHeatTransfer
from BaseEngineCycle.Nozzle import BellNozzle, ConicalNozzle
from BaseEngineCycle.Pressurant import Pressurant, PressurantTank
from BaseEngineCycle.Propellant import Propellant
from BaseEngineCycle.Tank import Tank
from BaseEngineCycle.ThrustChamber import ThrustChamber
from BaseEngineCycle.Pump import Pump
from BaseEngineCycle.FlowState import FlowState, DefaultFlowState
from cea import get_cea_values_dict
from cea_new import get_cea_dict, get_cea_chamber_dict
from irt import get_kerckhove, get_expansion_ratio_from_p_ratio, get_pressure_ratio_fsolve


@dataclass
class EngineCycle:
    thrust: float  # [N]
    burn_time: float  # [s]
    combustion_chamber_pressure: float  # [Pa]
    max_acceleration: float  # [m/s]
    mass_mixture_ratio: float  # [-]
    fuel_initial_pressure: float  # [Pa]
    fuel_initial_temperature: float  # [k]
    fuel_pump_pressure_factor: float  # [Pa]
    fuel_pump_specific_power: float  # [W]
    fuel_pump_efficiency: float  # [-]
    oxidizer_initial_pressure: float  # [Pa]
    oxidizer_initial_temperature: float  # [k]
    oxidizer_pump_pressure_factor: float  # [Pa]
    oxidizer_pump_specific_power: float  # [W]
    oxidizer_pump_efficiency: float  # [-]
    propellant_margin_factor: float  # [-]
    ullage_volume_factor: float  # [-]
    tanks_structural_factor: float  # [-]
    tanks_material_density: float  # [kg/m3]
    tanks_yield_strength: float  # [Pa]
    pressurant_heat_capacity_ratio: float  # [-]
    pressurant_molar_mass: float  # [kg/mol]
    pressurant_initial_temperature: float  # [K]
    pressurant_margin_factor: float  # [-]
    pressurant_initial_pressure: float  # [Pa]
    pressurant_final_pressure: float  # [Pa]
    pressurant_tank_material_density: float  # [kg/m3]
    pressurant_tank_yield_strength: float  # [Pa]
    pressurant_tank_safety_factor: float  # [-]
    combustion_chamber_material_density: float  # [kg/m3]
    combustion_chamber_yield_strength: float  # [Pa]
    combustion_chamber_safety_factor: float  # [-]
    injector_material_density: float  # [kg/m3]
    injector_yield_strength: float  # [Pa]
    injector_safety_factor: float  # [-]
    injector_propellant_is_gas: bool
    convergent_half_angle: float  # [rad]
    convergent_throat_bend_ratio: float  # [-]
    convergent_chamber_bend_ratio: float  # [-]
    divergent_throat_half_angle: float  # [rad]
    divergent_exit_half_angle: float  # [rad]
    maximum_wall_temperature: float  # [K]
    thrust_chamber_wall_emissivity: float  # [-]
    hot_gas_emissivity: float  # [-]

    nozzle_type: str
    convective_coefficient_mode: str
    oxidizer_name: str
    fuel_name: str
    is_frozen: bool

    # Values that override other inputs (one of them is required)
    expansion_ratio: Optional[float] = None  # [-]
    pressure_ratio: Optional[float] = None  # [-]
    exit_pressure_forced: Optional[float] = None  # [Pa]

    # Values that can be estimated or are not necessarily required
    recovery_factor: Optional[float] = None  # [-]
    area_ratio_chamber_throat: Optional[float] = None  # [-]
    chamber_characteristic_length: Optional[float] = None  # [m]
    coolant_inlet_temperature: Optional[float] = None  # [K]
    expansion_ratio_end_cooling: Optional[float] = None  # [-]
    distance_from_throat_end_cooling: Optional[float] = None  # [m]
    distance_from_throat_start_cooling: Optional[float] = None  # [m]
    cooling_section_pressure_drop: Optional[float] = None  # [Pa]
    ambient_pressure: Optional[float] = None  # [Pa]

    # Values that can be estimated by CEA
    characteristic_velocity: Optional[float] = None  # [m/s]
    ideal_thrust_coefficient: Optional[float] = None  # [-]
    combustion_temperature: Optional[float] = None  # [K]
    cc_hot_gas_molar_mass: Optional[float] = None  # [kg/mol]
    cc_hot_gas_heat_capacity_ratio: Optional[float] = None  # [-]
    cc_hot_gas_dynamic_viscosity: Optional[float] = None  # [Pa*s]
    cc_hot_gas_prandtl_number: Optional[float] = None  # [-]
    cc_hot_gas_specific_heat_capacity: Optional[float] = None  # [J/(kg*K)]

    _injector_pressure_drop_factor: float = .3
    iteration_accuracy = 0.0001

    verbose: bool = True
    fast_init: bool = False  # If True needs to make less calls to rocketCEA, assumes expansion ratio is provided, ignores pressure_ratio and exit_pressure_input
    iterate: bool = True

    def __post_init__(self):
        # Setting of internal variables
        self._cea_frozen, self._cea_frozenAtThroat = (1, 1) if self.is_frozen else (0, 0)
        if self.fast_init:
            assert self.expansion_ratio is not None
            self.set_cea()
            self.pressure_ratio = get_pressure_ratio_fsolve(self.expansion_ratio, self.cc_hot_gas_heat_capacity_ratio)
        else:
            self.resolve_expansion_pressure_ratio_choice()
            # Initiate CEA call if any of the CEA values is not provided
            if any(getattr(self, key) is None for key in self.cea_dict):
                self.set_cea()

    def resolve_expansion_pressure_ratio_choice(self):
        # Estimating the expansion ratio based on the pressure ratio or vice versa, exit_pressure_forced overrides both
        # TODO: This breaks updating of either variable as the other will not update
        if self.cc_hot_gas_heat_capacity_ratio is None:
            self.cc_hot_gas_heat_capacity_ratio = self.get_heat_capacity_ratio()  # Required for following calculations
        if self.exit_pressure_forced is not None:
            if self.verbose:
                warnings.warn('Exit pressure is given, pressure- and expansion ratio are ignored if provided')
            self.pressure_ratio = self.combustion_chamber_pressure / self.exit_pressure_forced
            self.expansion_ratio = get_expansion_ratio_from_p_ratio(self.pressure_ratio,
                                                                    self.cc_hot_gas_heat_capacity_ratio)
        elif not (self.pressure_ratio is None) ^ (self.expansion_ratio is None):
            raise ValueError(
                'Neither or both the pressure_ratio and expansion_ratio are given. Provide one and only one')
        elif self.pressure_ratio is None:
            self.pressure_ratio = get_pressure_ratio_fsolve(self.expansion_ratio,
                                                            self.cc_hot_gas_heat_capacity_ratio)
        elif self.expansion_ratio is None:
            self.expansion_ratio = get_expansion_ratio_from_p_ratio(self.pressure_ratio,
                                                                    self.cc_hot_gas_heat_capacity_ratio)

    @cached_property
    def cea_dict(self):
        return {'characteristic_velocity': 'c_star',
                'ideal_thrust_coefficient': 'C_F',
                'combustion_temperature': 'T_C',
                'cc_hot_gas_molar_mass': 'mm_cc',
                'cc_hot_gas_heat_capacity_ratio': 'y_cc',
                'cc_hot_gas_dynamic_viscosity': 'mu_cc',
                'cc_hot_gas_prandtl_number': 'pr_cc',
                'cc_hot_gas_specific_heat_capacity': 'cp_cc'}

    @property
    def cea_kwargs(self):
        # CEA values always calculated from pressure ratio
        # TODO: Give option to calculate CEA values with given expansion ratio (eps)
        return {'Pc': self.combustion_chamber_pressure,
                'MR': self.mass_mixture_ratio,
                'eps': None,
                'PcOvPe': self.pressure_ratio,
                'fuelName': self.fuel_name,
                'oxName': self.oxidizer_name,
                'frozen': self._cea_frozen,
                'frozenAtThroat': self._cea_frozenAtThroat}

    def set_cea(self):
        # Checking if value is given, if not: assign value found by CEA. Despite this check for each attribute, it is
        # recommended to either provide none of the CEA properties or all of them
        cea_attributes = self.cea_dict.keys()
        cea_values = get_cea_dict(**self.cea_kwargs)
        for attribute in cea_attributes:
            cea_name = self.cea_dict[attribute]
            if getattr(self, attribute) is None:
                setattr(self, attribute, cea_values[cea_name])

    def update_cea(self):
        cea_values = self.get_cea()
        for attribute in self.cea_dict.keys():
            cea_name = self.cea_dict[attribute]
            setattr(self, attribute, cea_values[cea_name])

    def get_heat_capacity_ratio(self):
        # Get the heat_capacity_ratio before an expansion ratio or pressure ratio is provided,
        #  required for setting all CEA values at once
        kwargs = {key: value for key, value in self.cea_kwargs.items() if key not in ('eps', 'PcOvPe')}
        return get_cea_chamber_dict(**kwargs)['y_cc']

    def do_pressure_check(self):
        self.get_pressure_check()

    def get_pressure_check(self, extra_pumps: tuple[Pump, ...] = ()):
        if self.verbose:
            pumps = (self.fuel_pump, self.oxidizer_pump, *extra_pumps)
            expected_outlet_pressures = (self.fuel_pump_expected_pressure, self.oxidizer_pump_expected_pressure)
            for pump, p_expect in zip(pumps, expected_outlet_pressures):
                if pump.outlet_pressure != p_expect:
                    warnings.warn(f'For the {pump.inlet_flow_state.type}-pump the used outlet pressure '
                                  f'[{pump.outlet_pressure:.4e} Pa] is not equal to estimated expected outlet pressure '
                                  f'[{p_expect:.4e} Pa]\n'
                                  f'{pump}\n')

    @property
    def fuel_pump_expected_pressure(self):
        return (self.combustion_chamber_pressure
                - self.injector.pressure_change
                - self.cooling_channel_section.pressure_change)

    @property
    def oxidizer_pump_expected_pressure(self):
        return (self.combustion_chamber_pressure
                - self.injector.pressure_change)

    @property
    def cstar_cf(self):
        return self.characteristic_velocity, self.ideal_thrust_coefficient

    @property
    def propellant_mix_name(self):
        return get_propellant_mix_name(fuel_name=self.fuel_name, oxidizer_name=self.oxidizer_name)

    @property
    def base_mass_flow(self):
        return self.thrust / (self.characteristic_velocity * self.ideal_thrust_coefficient)

    @property
    def chamber_mass_flow(self):
        return self.base_mass_flow

    @property
    def throat_area(self):
        return self.chamber_mass_flow * sqrt(constants.R / self.cc_hot_gas_molar_mass * self.combustion_temperature) / (
                get_kerckhove(self.cc_hot_gas_heat_capacity_ratio) * self.combustion_chamber_pressure)

    @property
    def exit_area(self):
        return self.throat_area * self.expansion_ratio

    @property
    def exit_pressure(self):
        return self.combustion_chamber_pressure / self.pressure_ratio

    @property
    def simple_specific_impulse(self):
        return self.thrust / self.chamber_mass_flow / g

    @property
    def chamber_fuel_flow(self):
        return 1 / (self.mass_mixture_ratio + 1) * self.chamber_mass_flow

    @property
    def chamber_oxidizer_flow(self):
        return self.mass_mixture_ratio / (self.mass_mixture_ratio + 1) * self.chamber_mass_flow

    @property
    def main_fuel_flow(self):  # Default to chamber flow, overriden in child classes/cycles
        return self.chamber_fuel_flow

    @property
    def main_oxidizer_flow(self):  # Default to chamber flow, overriden in child classes/cycles
        return self.chamber_oxidizer_flow

    @property
    def delta_p_oxidizer_pump(self):
        return self.combustion_chamber_pressure * self.oxidizer_pump_pressure_factor - self.oxidizer_initial_pressure

    @property
    def delta_p_fuel_pump(self):
        return self.combustion_chamber_pressure * self.fuel_pump_pressure_factor - self.fuel_initial_pressure

    @property
    def oxidizer_initial_flow_state(self):
        return FlowState(propellant_name=self.oxidizer_name,
                         temperature=self.oxidizer_initial_temperature,
                         pressure=self.oxidizer_initial_pressure,
                         mass_flow=self.main_oxidizer_flow,
                         type='oxidizer',)

    @property
    def fuel_initial_flow_state(self):
        return FlowState(propellant_name=self.fuel_name,
                         temperature=self.fuel_initial_temperature,
                         pressure=self.fuel_initial_pressure,
                         mass_flow=self.main_fuel_flow,
                         type='fuel',)

    @property
    def oxidizer(self):
        return Propellant(initial_flow_state=self.oxidizer_initial_flow_state,
                          burn_time=self.burn_time,
                          margin_factor=self.propellant_margin_factor)

    @property
    def fuel(self):
        return Propellant(initial_flow_state=self.fuel_initial_flow_state,
                          burn_time=self.burn_time,
                          margin_factor=self.propellant_margin_factor)

    @property
    def pressurant(self):
        return Pressurant(oxidizer_volume=self.oxidizer.volume,
                          fuel_volume=self.fuel.volume,
                          fuel_tank_initial_pressure=self.fuel_initial_pressure,
                          oxidizer_tank_initial_pressure=self.oxidizer_initial_pressure,
                          margin_factor=self.pressurant_margin_factor,
                          initial_pressure=self.pressurant_initial_pressure,
                          final_pressure=self.pressurant_final_pressure,
                          heat_capacity_ratio=self.pressurant_heat_capacity_ratio,
                          molar_mass=self.pressurant_molar_mass,
                          initial_temperature=self.pressurant_initial_temperature,
                          propellant_tanks_ullage_factor=self.ullage_volume_factor)

    @property
    def pressurant_tank(self):
        return PressurantTank(material_density=self.pressurant_tank_material_density,
                              safety_factor=self.pressurant_tank_safety_factor,
                              yield_strength=self.pressurant_tank_yield_strength,
                              pressurant=self.pressurant)

    @property
    def oxidizer_tank(self):
        return Tank(inlet_flow_state=self.oxidizer_initial_flow_state,
                    propellant_volume=self.oxidizer.volume,
                    max_acceleration=self.max_acceleration,
                    ullage_factor=self.ullage_volume_factor,
                    pressurant_tank_volume=self.pressurant_tank.volume,
                    material_density=self.tanks_material_density,
                    yield_strength=self.tanks_yield_strength,
                    safety_factor=self.tanks_structural_factor)

    @property
    def fuel_tank(self):
        return Tank(inlet_flow_state=self.fuel_initial_flow_state,
                    propellant_volume=self.fuel.volume,
                    max_acceleration=self.max_acceleration,
                    ullage_factor=self.ullage_volume_factor,
                    pressurant_tank_volume=None,
                    material_density=self.tanks_material_density,
                    yield_strength=self.tanks_yield_strength,
                    safety_factor=self.tanks_structural_factor, )

    @property
    def oxidizer_pump(self):
        return Pump(inlet_flow_state=self.oxidizer_tank.outlet_flow_state,
                    pressure_increase=self.delta_p_oxidizer_pump,
                    efficiency=self.oxidizer_pump_efficiency,
                    specific_power=self.oxidizer_pump_specific_power,)

    @property
    def fuel_pump(self):
        return Pump(inlet_flow_state=self.fuel_tank.outlet_flow_state,
                    pressure_increase=self.delta_p_fuel_pump,
                    efficiency=self.fuel_pump_efficiency,
                    specific_power=self.fuel_pump_specific_power,)

    @property
    def injector_inlet_flow_states(self):
        """Ugly hack to prevent recursion (CoolingChannelSection indirectly requires ThrustChamber, which requires
                Injector, which would require CoolingChannelSection without this hack, creating an infinite loop)"""
        if CoolingChannelSection._instance_created:
            return get_injector_inlet_flow_state_fuel, self.oxidizer_pump.outlet_flow_state
        else:
            return DefaultFlowState(), self.oxidizer_pump.outlet_flow_state

    @property
    def get_injector_inlet_flow_state_fuel(self) -> FlowState:
        """part of ugly hack"""
        return self.cooling_channel_section.outlet_flow_state

    @property
    def injector(self):
        # noinspection PyArgumentList
        return Injector(inlet_flow_states=self.injector_inlet_flow_states,
                        combustion_chamber_pressure=self.combustion_chamber_pressure,
                        combustion_chamber_area=self.combustion_chamber.area,
                        material_density=self.injector_material_density,
                        safety_factor=self.injector_safety_factor,
                        yield_strength=self.injector_yield_strength,
                        propellant_is_gas=self.injector_propellant_is_gas,
                        _pressure_drop_factor=.3)

    @property
    def nozzle(self):
        kwargs = {'throat_area': self.throat_area,
                  'expansion_ratio': self.expansion_ratio,
                  'area_ratio_chamber_throat': self.area_ratio_chamber_throat,
                  'conv_half_angle': self.convergent_half_angle,
                  'conv_throat_bend_ratio': self.convergent_throat_bend_ratio,
                  'conv_chamber_bend_ratio': self.convergent_chamber_bend_ratio,
                  'div_throat_half_angle': self.divergent_throat_half_angle}
        if self.nozzle_type == 'conical':
            return ConicalNozzle(**kwargs)
        elif self.nozzle_type == 'bell':
            return BellNozzle(**kwargs, div_exit_half_angle=self.divergent_exit_half_angle)
        else:
            raise ValueError('No proper nozzle_type was provided, choose from [bell, conical]')

    @property
    def combustion_chamber(self):
        return CombustionChamber(throat_area=self.throat_area,
                                 combustion_chamber_pressure=self.combustion_chamber_pressure,
                                 propellant_mix=self.propellant_mix_name,
                                 area_ratio_chamber_throat=self.area_ratio_chamber_throat,
                                 characteristic_length=self.chamber_characteristic_length,
                                 convergent_volume_estimate=self.nozzle.conv_volume_estimate,
                                 safety_factor=self.combustion_chamber_safety_factor,
                                 yield_strength=self.combustion_chamber_yield_strength,
                                 material_density=self.combustion_chamber_material_density,
                                 verbose=self.verbose)

    @property
    def thrust_chamber(self):
        return ThrustChamber(nozzle=self.nozzle,
                             chamber=self.combustion_chamber,
                             injector=self.injector,
                             heat_capacity_ratio=self.cc_hot_gas_heat_capacity_ratio)

    @property
    def max_distance_from_throat_heat_transfer_section(self):
        if self.expansion_ratio_end_cooling:
            if self.verbose:
                warnings.warn(
                    'Expansion_ratio_end_cooling is given, distance_from_throat_end_cooling is ignored if provided')
            return self.thrust_chamber.get_distance_for_divergent_expansion_ratio(self.expansion_ratio_end_cooling)
        else:
            return self.distance_from_throat_end_cooling

    @property
    def min_distance_from_throat_heat_transfer_section(self):
        return self.distance_from_throat_start_cooling

    @property
    def convective_heat_transfer_args(self):
        return {'combustion_temperature': self.combustion_temperature,
                'combustion_chamber_pressure': self.combustion_chamber_pressure,
                'dynamic_viscosity': self.cc_hot_gas_dynamic_viscosity,
                'specific_heat_capacity': self.cc_hot_gas_specific_heat_capacity,
                'mass_flow': self.chamber_mass_flow,
                'maximum_wall_temperature': self.maximum_wall_temperature,
                'thrust_chamber_wall_emissivity': self.thrust_chamber_wall_emissivity,
                'hot_gas_emissivity': self.hot_gas_emissivity,
                'heat_capacity_ratio': self.cc_hot_gas_heat_capacity_ratio,
                'convective_coefficient_mode': self.convective_coefficient_mode,
                'thrust_chamber': self.thrust_chamber,
                'recovery_factor': self.recovery_factor,
                'prandtl_number': self.cc_hot_gas_prandtl_number,
                'verbose': self.verbose, }

    @cached_property
    def theoretical_convective_heat_transfer(self):
        return ConvectiveHeatTransfer(**self.convective_heat_transfer_args)

    @cached_property
    def radiative_heat_transfer(self):
        return RadiativeHeatTransfer(
            thrust_chamber=self.thrust_chamber,
            combustion_temperature=self.combustion_temperature,
            maximum_wall_temperature=self.maximum_wall_temperature,
            thrust_chamber_wall_emissivity=self.thrust_chamber_wall_emissivity,
            hot_gas_emissivity=self.hot_gas_emissivity,
            theoretical_total_convective_heat_transfer=self.theoretical_convective_heat_transfer.total_convective_heat_transfer, )

    @cached_property
    def heat_transfer_section(self):
        return HeatTransferSection(**self.convective_heat_transfer_args,
                                   max_distance_section=self.max_distance_from_throat_heat_transfer_section,
                                   min_distance_section=self.min_distance_from_throat_heat_transfer_section,
                                   radiative_heat_transfer_factor=self.radiative_heat_transfer.radiative_factor)

    @property
    def cooling_inlet_flow_state(self):
        return self.fuel_pump.outlet_flow_state

    @property
    def cooling_channel_section(self):
        return CoolingChannelSection(inlet_flow_state=self.cooling_inlet_flow_state,
                                     heat_transfer_section=self.heat_transfer_section,
                                     _total_heat_transfer=None,
                                     combustion_chamber_pressure=self.combustion_chamber_pressure,
                                     pressure_drop=self.cooling_section_pressure_drop,
                                     verbose=self.verbose,)

    @property
    def pump_power_required(self):
        return self.fuel_pump.power_required + self.oxidizer_pump.power_required

    @property
    def pumps_mass(self):
        return self.fuel_pump.mass + self.oxidizer_pump.mass

    @property
    def tanks_mass(self):
        return self.fuel_tank.mass + self.oxidizer_tank.mass + self.pressurant_tank.mass

    @property
    def props_mass(self):
        return self.oxidizer.mass + self.fuel.mass

    @property
    def mass(self):
        return self.props_mass + self.tanks_mass + self.pumps_mass

    @property
    def mass_ratio(self):
        return (self.mass - self.props_mass) / self.mass

    @property
    def ideal_delta_v(self):
        return self.simple_specific_impulse * log(1 / self.mass_ratio) * g

    @property
    def gravity_delta_v(self, vertical_fraction: float = 0.2):
        return self.ideal_delta_v - constants.g * self.burn_time * vertical_fraction

    @property
    def chamber_ideal_specific_impulse(self):
        return self.ideal_thrust_coefficient * self.characteristic_velocity / constants.g

    @property
    def chamber_vacuum_specific_impulse(self):
        return self.characteristic_velocity / constants.g * (self.ideal_thrust_coefficient + self.exit_pressure / self.combustion_chamber_pressure * self.expansion_ratio)

    @property
    def payload_delta_v(self):
        payload = 10  # [kg]
        return self.simple_specific_impulse * log(self.mass + payload / (self.mass - self.props_mass + payload)) * g

    @property
    def payload_mass_ratio(self):
        if self.mass_u is None:
            ValueError('Payload mass not given, impossible to calculate mass ratio with payload mass')
        return (self.mass + self.mass_u - self.props_mass) / (self.mass + self.mass_u)
