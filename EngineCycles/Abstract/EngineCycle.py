import warnings
from dataclasses import dataclass, field
from functools import cached_property
from math import log, exp
from typing import Optional
from CoolProp import CoolProp

from scipy import constants as constants

from EngineComponents.Base.CombustionChamber import CombustionChamber
from EngineComponents.Base.Injector import Injector
from EngineComponents.Base.Cooling import CoolingChannelSection
from EngineComponents.Base.HeatTransferSection import HeatTransferSection, ConvectiveHeatTransfer, \
    RadiativeHeatTransfer
from EngineComponents.Base.Nozzle import BellNozzle, ConicalNozzle
from EngineComponents.Base.Pressurant import Pressurant, PressurantTank
from EngineComponents.Base.Propellant import Propellant
from EngineComponents.Base.Tank import Tank
from EngineComponents.Base.ThrustChamber import ThrustChamber
from EngineComponents.Base.Pump import Pump
from EngineComponents.Base.Merger import Merger
from EngineComponents.Abstract.FlowState import FlowState
from EngineComponents.Abstract.Material import Material
from EngineFunctions.CEAFunctions import get_cea_dict, get_cea_chamber_dict
from EngineFunctions.IRTFunctions import get_expansion_ratio_from_p_ratio, \
    get_pressure_ratio_fsolve, get_throat_area, get_thrust_coefficient_from_ideal
from EngineFunctions.AssumeValueFunctions import get_characteristic_length, get_initial_propellant_temperature


@dataclass
class EngineCycle:
    thrust: float  # [N]
    burn_time: float  # [s]
    combustion_chamber_pressure: float  # [Pa]
    oxidizer_name: str
    fuel_name: str
    max_acceleration: float  # [m/s]
    mass_mixture_ratio: float  # [-]
    fuel_initial_pressure: float  # [Pa]
    fuel_pump_specific_power: float  # [W]
    fuel_pump_efficiency: float  # [-]
    oxidizer_initial_pressure: float  # [Pa]
    oxidizer_pump_specific_power: float  # [W]
    oxidizer_pump_efficiency: float  # [-]
    propellant_margin_factor: float  # [-]
    ullage_volume_factor: float  # [-]
    tanks_structural_factor: float  # [-]
    fuel_tank_material: Material
    oxidizer_tank_material: Material
    pressurant_heat_capacity_ratio: float  # [-]
    pressurant_molar_mass: float  # [kg/mol]
    pressurant_initial_temperature: float  # [K]
    pressurant_margin_factor: float  # [-]
    pressurant_initial_pressure: float  # [Pa]
    pressurant_final_pressure: float  # [Pa]
    pressurant_tank_material: Material
    pressurant_tank_safety_factor: float  # [-]
    combustion_chamber_material: Material
    combustion_chamber_safety_factor: float  # [-]
    injector_material: Material
    injector_safety_factor: float  # [-]
    cooling_pressure_drop_factor: float  # [-]
    injector_pressure_drop_factor: float  # [-]
    convergent_half_angle: float  # [rad]
    convergent_throat_bend_ratio: float  # [-]
    convergent_chamber_bend_ratio: float  # [-]
    divergent_throat_half_angle: float  # [rad]
    nozzle_material: Material
    nozzle_safety_factor: float  # [-]
    nozzle_type: str
    maximum_wall_temperature: float  # [K]
    thrust_chamber_wall_emissivity: float  # [-]
    hot_gas_emissivity: float  # [-]
    specific_impulse_correction_factor: float  # [-]
    shaft_mechanical_efficiency: float  # [-]

    is_frozen: bool = True
    # Values that override other inputs (one of them is required)
    expansion_ratio: Optional[float] = None  # [-]
    pressure_ratio: Optional[float] = None  # [-]
    exit_pressure_forced: Optional[float] = None  # [Pa]

    # Values that can be estimated or are not necessarily required
    oxidizer_initial_temperature: Optional[float] = None  # [k]
    fuel_initial_temperature: Optional[float] = None  # [K]
    divergent_exit_half_angle: Optional[float] = None  # [rad]
    recovery_factor: Optional[float] = None  # [-]
    area_ratio_chamber_throat: Optional[float] = None  # [-]
    chamber_characteristic_length: Optional[float] = None  # [m]
    coolant_inlet_temperature: Optional[float] = None  # [K]
    expansion_ratio_start_cooling: Optional[float] = None
    expansion_ratio_end_cooling: Optional[float] = None  # [-]
    distance_from_throat_end_cooling: Optional[float] = None  # [m]
    distance_from_throat_start_cooling: Optional[float] = None  # [m]
    cooling_section_pressure_drop: Optional[float] = None  # [Pa]
    ambient_pressure: Optional[float] = None  # [Pa]
    fuel_density: Optional[float] = None  # [kg/m3]
    oxidizer_density: Optional[float] = None  # [kg/m3]

    # Values that can be estimated by CEA
    characteristic_velocity: Optional[float] = None  # [m/s]
    ideal_thrust_coefficient: Optional[float] = None  # [-]
    combustion_temperature: Optional[float] = None  # [K]
    cc_hot_gas_molar_mass: Optional[float] = None  # [kg/mol]
    cc_hot_gas_heat_capacity_ratio: Optional[float] = None  # [-]
    cc_hot_gas_dynamic_viscosity: Optional[float] = None  # [Pa*s]
    cc_hot_gas_prandtl_number: Optional[float] = None  # [-]
    cc_hot_gas_specific_heat_capacity: Optional[float] = None  # [J/(kg*K)]

    _oxidizer_pump_pressure_factor_first_guess: float = 1.15  # [Pa]
    _fuel_pump_pressure_factor_first_guess: float = 1.55  # [Pa]
    _is_temp_calc_needed: bool = True
    _ignore_cooling: bool = False
    iteration_accuracy: float = 0.0001

    # Values always calculated by program, but need to be saved as attributes, not properties
    total_heat_transfer: float = field(init=False, repr=False, default=0)
    minimum_required_coolant_mass_flow: float = field(init=False, repr=False, default=0)
    _fuel_pump_outlet_pressure: float = field(init=False, repr=False, default=None)
    _oxidizer_pump_outlet_pressure: float = field(init=False, repr=False, default=None)
    _expansion_ratio_end_cooling: float = field(init=False, repr=False, default=None)
    verbose: bool = True
    iterate: bool = True

    def __post_init__(self):
        self.initialize_cea()
        self.set_initial_values()
        if self._ignore_cooling:
            self._is_temp_calc_needed = False
            warnings.warn('!!_ignore_cooling flag has been set to True!!')
            self.set_pump_outlet_pressures()
        else:
            self.set_heat_transfer()
            self.set_pump_outlet_pressures()
            self.calc_minimum_required_coolant_mass_flow()

    def initialize_cea(self):
        # Setting of internal variables
        self._cea_frozen, self._cea_frozenAtThroat = (1, 1) if self.is_frozen else (0, 0)
        self.resolve_expansion_choice()
        self.set_cea()
        if self.expansion_ratio is None:
            self.expansion_ratio = get_expansion_ratio_from_p_ratio(self.pressure_ratio,
                                                                    self.cc_hot_gas_heat_capacity_ratio)

    def resolve_expansion_choice(self):
        if self.exit_pressure_forced is not None:
            if self.verbose:
                warnings.warn('Exit pressure is given, pressure- and expansion ratio are ignored if provided')
            self.pressure_ratio = self.combustion_chamber_pressure / self.exit_pressure_forced
        elif not ((self.pressure_ratio is None) ^ (self.expansion_ratio is None)):
            raise ValueError(
                'Neither or both the pressure_ratio and expansion_ratio are given. Provide one and only one')
        elif self.pressure_ratio is None:
            self.cc_hot_gas_heat_capacity_ratio = self.get_heat_capacity_ratio()
            self.pressure_ratio = get_pressure_ratio_fsolve(self.expansion_ratio,
                                                            self.cc_hot_gas_heat_capacity_ratio)

    def set_initial_values(self):
        """Used by child classes as well to init values directly after CEA is set, but before other operations"""
        if self.fuel_initial_temperature is None:
            self.fuel_initial_temperature = get_initial_propellant_temperature(
                propellant_name=self.fuel_name
            )
        if self.oxidizer_initial_temperature is None:
            self.oxidizer_initial_temperature = get_initial_propellant_temperature(
                propellant_name=self.oxidizer_name
            )
        if self.chamber_characteristic_length is None:
            self.chamber_characteristic_length = get_characteristic_length(
                fuel_name=self.fuel_name,
                oxidizer_name=self.oxidizer_name
            )

    def set_heat_transfer(self):
        self.total_heat_transfer = self.heat_transfer_section.total_heat_transfer
        self.heat_transfer_func = self.heat_transfer_section.init_heat_transfer()

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

    def set_pump_outlet_pressures(self):
        Merger._warn_pressure = False
        self.calc_pump_outlet_pressures()
        Merger._warn_pressure = True

    def calc_pump_outlet_pressures(self):
        self._fuel_pump_outlet_pressure = self.fuel_pump_expected_pressure
        self._oxidizer_pump_outlet_pressure = self.oxidizer_pump_expected_pressure

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
    def propellant_mix_name(self):
        return get_propellant_mix_name(fuel_name=self.fuel_name, oxidizer_name=self.oxidizer_name)
    
    @property
    def thrust_coefficient(self):
        return get_thrust_coefficient_from_ideal(ideal_thrust_coefficient=self.ideal_thrust_coefficient,
                                                 chamber_pressure=self.combustion_chamber_pressure,
                                                 exit_pressure=self.exit_pressure,
                                                 expansion_ratio=self.expansion_ratio,
                                                 ambient_pressure=self.ambient_pressure,)

    @cached_property
    def chamber_specific_impulse(self):
        return self.thrust_coefficient * self.characteristic_velocity * self.specific_impulse_correction_factor / constants.g

    @property
    def base_mass_flow(self):
        return self.thrust / (self.chamber_specific_impulse * constants.g)

    @property
    def chamber_mass_flow(self):
        return self.chamber_thrust / (self.chamber_specific_impulse * constants.g)

    @property
    def throat_area(self):
        return get_throat_area(molar_mass=self.cc_hot_gas_molar_mass,
                               heat_capacity_ratio=self.cc_hot_gas_heat_capacity_ratio,
                               chamber_temperature=self.combustion_temperature,
                               mass_flow=self.chamber_mass_flow,
                               chamber_pressure=self.combustion_chamber_pressure)

    @property
    def exit_area(self):
        return self.throat_area * self.expansion_ratio

    @property
    def exit_pressure(self):
        return self.combustion_chamber_pressure / self.pressure_ratio

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
    def total_mass_flow(self):
        return self.main_oxidizer_flow + self.main_fuel_flow

    @property
    def oxidizer_pump_outlet_pressure(self):
        # Initial engine run with estimate using pressure_factor, after which outlet_pressure will be known and used instead
        if self._oxidizer_pump_outlet_pressure:
            return self._oxidizer_pump_outlet_pressure
        else:
            return self.combustion_chamber_pressure * self._oxidizer_pump_pressure_factor_first_guess

    @property
    def fuel_pump_outlet_pressure(self):
        # Initial engine run with estimate using pressure_factor, after which outlet_pressure will be known and used instead
        if self._fuel_pump_outlet_pressure:
            return self._fuel_pump_outlet_pressure
        else:
            return self.combustion_chamber_pressure * self._fuel_pump_pressure_factor_first_guess

    @property
    def oxidizer_initial_flow_state(self):
        return FlowState(propellant_name=self.oxidizer_name,
                         temperature=self.oxidizer_initial_temperature,
                         pressure=self.oxidizer_initial_pressure,
                         mass_flow=self.main_oxidizer_flow,
                         type='oxidizer', )

    @property
    def fuel_initial_flow_state(self):
        return FlowState(propellant_name=self.fuel_name,
                         temperature=self.fuel_initial_temperature,
                         pressure=self.fuel_initial_pressure,
                         mass_flow=self.main_fuel_flow,
                         type='fuel', )

    @property
    def oxidizer(self):
        return Propellant(initial_flow_state=self.oxidizer_initial_flow_state,
                          burn_time=self.burn_time,
                          margin_factor=self.propellant_margin_factor,
                          propellant_density=self.oxidizer_density)

    @property
    def fuel(self):
        return Propellant(initial_flow_state=self.fuel_initial_flow_state,
                          burn_time=self.burn_time,
                          margin_factor=self.propellant_margin_factor,
                          propellant_density=self.fuel_density)

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
        return PressurantTank(structure_material=self.pressurant_tank_material,
                              safety_factor=self.pressurant_tank_safety_factor,
                              pressurant=self.pressurant)

    @property
    def oxidizer_tank(self):
        return Tank(inlet_flow_state=self.oxidizer_initial_flow_state,
                    propellant_volume=self.oxidizer.volume,
                    max_acceleration=self.max_acceleration,
                    ullage_factor=self.ullage_volume_factor,
                    pressurant_tank_volume=self.pressurant_tank.volume,
                    structure_material=self.oxidizer_tank_material,
                    safety_factor=self.tanks_structural_factor,
                    propellant_density=self.oxidizer_density, )

    @property
    def fuel_tank(self):
        return Tank(inlet_flow_state=self.fuel_initial_flow_state,
                    propellant_volume=self.fuel.volume,
                    max_acceleration=self.max_acceleration,
                    ullage_factor=self.ullage_volume_factor,
                    pressurant_tank_volume=None,
                    structure_material=self.fuel_tank_material,
                    safety_factor=self.tanks_structural_factor,
                    propellant_density=self.fuel_density, )

    @property
    def oxidizer_pump(self):
        return Pump(inlet_flow_state=self.oxidizer_tank.outlet_flow_state,
                    expected_outlet_pressure=self.oxidizer_pump_outlet_pressure,
                    efficiency=self.oxidizer_pump_efficiency,
                    specific_power=self.oxidizer_pump_specific_power,
                    propellant_density=self.oxidizer_density, )

    @property
    def fuel_pump(self):
        return Pump(inlet_flow_state=self.fuel_tank.outlet_flow_state,
                    expected_outlet_pressure=self.fuel_pump_outlet_pressure,
                    efficiency=self.fuel_pump_efficiency,
                    specific_power=self.fuel_pump_specific_power,
                    propellant_density=self.fuel_density, )

    @property
    def injector_inlet_flow_states(self):
        return self.cooling_channel_section.outlet_flow_state, self.oxidizer_pump.outlet_flow_state

    @property
    def injector(self):
        return Injector(inlet_flow_states=self.injector_inlet_flow_states,
                        combustion_chamber_pressure=self.combustion_chamber_pressure,
                        combustion_chamber_area=self.combustion_chamber.area,
                        structure_material=self.injector_material,
                        safety_factor=self.injector_safety_factor,
                        pressure_drop_factor=self.injector_pressure_drop_factor)

    @cached_property
    def nozzle(self):
        kwargs = {'throat_area': self.throat_area,
                  'expansion_ratio': self.expansion_ratio,
                  'area_ratio_chamber_throat': self.area_ratio_chamber_throat,
                  'conv_half_angle': self.convergent_half_angle,
                  'conv_throat_bend_ratio': self.convergent_throat_bend_ratio,
                  'conv_chamber_bend_ratio': self.convergent_chamber_bend_ratio,
                  'div_throat_half_angle': self.divergent_throat_half_angle,
                  'structure_material': self.nozzle_material,
                  'chamber_pressure': self.combustion_chamber_pressure,
                  'safety_factor': self.nozzle_safety_factor,}
        if self.nozzle_type == 'conical':
            return ConicalNozzle(**kwargs)
        elif self.nozzle_type == 'bell':
            return BellNozzle(**kwargs, div_exit_half_angle=self.divergent_exit_half_angle)
        else:
            raise ValueError('No proper nozzle_type was provided, choose from [bell, conical]')

    @cached_property
    def combustion_chamber(self):
        return CombustionChamber(throat_area=self.throat_area,
                                 combustion_chamber_pressure=self.combustion_chamber_pressure,
                                 area_ratio_chamber_throat=self.area_ratio_chamber_throat,
                                 characteristic_length=self.chamber_characteristic_length,
                                 convergent_volume_estimate=self.nozzle.conv_volume_estimate,
                                 safety_factor=self.combustion_chamber_safety_factor,
                                 structure_material=self.combustion_chamber_material,)

    @cached_property
    def thrust_chamber(self):
        return ThrustChamber(nozzle=self.nozzle,
                             chamber=self.combustion_chamber,
                             heat_capacity_ratio=self.cc_hot_gas_heat_capacity_ratio)

    @property
    def max_distance_from_throat_heat_transfer_section(self):
        if self.expansion_ratio_end_cooling:
            if self.verbose and self.distance_from_throat_end_cooling:
                warnings.warn(
                    'Expansion_ratio_end_cooling is given, distance_from_throat_end_cooling is ignored, but also provided')
            return self.thrust_chamber.get_distance_for_divergent_expansion_ratio(self.expansion_ratio_end_cooling)
        else:
            if self.distance_from_throat_end_cooling is None:
                if self.expansion_ratio > 20:
                    warnings.warn('No end of cooling provided, limited to expansion ratio of 20')
                    self._expansion_ratio_end_cooling = 20
                    return self.thrust_chamber.get_distance_for_divergent_expansion_ratio(20)
                else:
                    warnings.warn('No end of cooling provided, assumed to be end of nozzle')
                    self._expansion_ratio_end_cooling = self.expansion_ratio
                    return self.thrust_chamber.max_distance_from_throat
            else:
                r_end = self.thrust_chamber.get_radius(self.distance_from_throat_end_cooling)
                self._expansion_ratio_end_cooling = r_end**2 * constants.pi / self.throat_area
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
                'hot_gas_emissivity': self.hot_gas_emissivity,
                'heat_capacity_ratio': self.cc_hot_gas_heat_capacity_ratio,
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
    def maximum_coolant_outlet_temperature(self):
        return self.maximum_wall_temperature

    def calc_minimum_required_coolant_mass_flow(self):
        """Determine minimum coolant flow required to keep outlet temp below maximum."""
        ccs_flow_state = self.cooling_inlet_flow_state
        h_max = CoolProp.PropsSI('H',
                                 'T', self.maximum_coolant_outlet_temperature,
                                 'P', ccs_flow_state.pressure,
                                 ccs_flow_state.coolprop_name)
        h_in = ccs_flow_state.mass_specific_enthalpy
        delta_h_max = h_max - h_in
        q_tot = self.total_heat_transfer
        self.minimum_required_coolant_mass_flow = q_tot / abs(delta_h_max)

        if self.minimum_required_coolant_mass_flow > self.main_fuel_flow:
            raise ValueError(
                'The minimum required coolant flow is larger than the main fuel flow: cooling is not possible')

    @property
    def cooling_inlet_flow_state(self):
        return self.fuel_pump.outlet_flow_state

    @property
    def cooling_channel_section(self):
        return CoolingChannelSection(inlet_flow_state=self.cooling_inlet_flow_state,
                                     total_heat_transfer=self.total_heat_transfer,
                                     maximum_outlet_temperature=self.maximum_coolant_outlet_temperature,
                                     combustion_chamber_pressure=self.combustion_chamber_pressure,
                                     pressure_drop=self.cooling_section_pressure_drop,
                                     verbose=self.verbose,
                                     _is_temp_calc_needed=self._is_temp_calc_needed,
                                     pressure_drop_factor=self.cooling_pressure_drop_factor,)

    @property
    def expansion_ratio_end(self):
        if self.expansion_ratio_end_cooling:
            return self.expansion_ratio_end_cooling
        else:
            return self._expansion_ratio_end_cooling

    @property
    def fuel_pumps_power_required(self):
        return self.fuel_pump.power_required / self.shaft_mechanical_efficiency

    @property
    def oxidizer_pumps_power_required(self):
        return self.oxidizer_pump.power_required / self.shaft_mechanical_efficiency

    @property
    def pumps_power_required(self):
        return self.fuel_pumps_power_required + self.oxidizer_pumps_power_required

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
        return self.props_mass + self.tanks_mass + self.pumps_mass + self.pressurant.mass

    @property
    def final_mass(self):
        return self.mass - self.props_mass

    @property
    def dry_mass(self):
        return self.pumps_mass + self.combustion_chamber.mass + self.nozzle.mass

    @property
    def mass_ratio(self):
        return self.final_mass / self.mass

    def get_payload(self, delta_v: float) -> float:
        e_dv = exp(delta_v / (self.overall_specific_impulse * constants.g))
        m0 = self.mass
        mf = self.final_mass
        return (m0 - mf * e_dv) / (e_dv - 1)

    @property
    def ideal_delta_v(self):
        return self.overall_specific_impulse * log(1 / self.mass_ratio) * constants.g

    @property
    def gravity_delta_v(self, vertical_fraction: float = 0.2):
        return self.ideal_delta_v - constants.g * self.burn_time * vertical_fraction

    @property
    def overall_specific_impulse(self):
        return self.thrust / self.total_mass_flow / constants.g

    @cached_property
    def vacuum_thrust_coefficient(self):
        return self.ideal_thrust_coefficient + self.exit_pressure / self.combustion_chamber_pressure * self.expansion_ratio

    @cached_property
    def sea_level_thrust_coefficient(self):
        sea_level_pressure = 101325
        return self.vacuum_thrust_coefficient - sea_level_pressure / self.combustion_chamber_pressure * self.expansion_ratio

    @cached_property
    def chamber_ideal_specific_impulse(self):
        return self.ideal_thrust_coefficient * self.characteristic_velocity / constants.g
    
    @cached_property
    def chamber_vacuum_specific_impulse(self):
        return self.vacuum_thrust_coefficient * self.characteristic_velocity / constants.g

    @cached_property
    def chamber_sea_level_specific_impulse(self):
        return self.sea_level_thrust_coefficient * self.characteristic_velocity / constants.g

    @property
    def payload_delta_v(self):
        payload = 10  # [kg]
        return self.overall_specific_impulse * log(
            self.mass + payload / (self.mass - self.props_mass + payload)) * constants.g

    @property
    def payload_mass_ratio(self):
        if self.mass_u is None:
            ValueError('Payload mass not given, impossible to calculate mass ratio with payload mass')
        return (self.mass + self.mass_u - self.props_mass) / (self.mass + self.mass_u)

    @property
    def chamber_thrust(self):
        """For consistency in naming with OpenCycles"""
        return self.thrust
