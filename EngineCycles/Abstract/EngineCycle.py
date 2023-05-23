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
from EngineComponents.Base.Nozzle import Nozzle
from EngineComponents.Base.Pressurant import Pressurant, PressurantTank
from EngineComponents.Base.Propellant import Propellant
from EngineComponents.Base.Tank import Tank
from EngineComponents.Base.ThrustChamber import ThrustChamber
from EngineComponents.Base.Pump import Pump
from EngineComponents.Base.Merger import Merger
from EngineComponents.Abstract.FlowState import FlowState, ManualFlowState
from EngineComponents.Abstract.Material import Material
from EngineFunctions.BaseFunctions import format_fancy_name
from EngineFunctions.CEAFunctions import get_cea_dict, get_cea_chamber_dict
from EngineFunctions.IRTFunctions import get_expansion_ratio_from_p_ratio, \
    get_pressure_ratio_fsolve, get_throat_area, get_thrust_coefficient_from_ideal
from EngineFunctions.AssumeValueFunctions import get_characteristic_length, get_initial_propellant_temperature, \
    get_prandtl_number_estimate, get_turbulent_recovery_factor, get_specific_impulse_quality_factor, \
    get_propellant_mixture, get_mass_mixture_ratio
from EngineFunctions.EmpiricalRelations import get_chamber_throat_area_ratio_estimate


@dataclass
class EngineCycle:
    thrust: float  # [N]
    burn_time: float  # [s]
    combustion_chamber_pressure: float  # [Pa]
    oxidizer_name: str
    fuel_name: str
    max_acceleration: float  # [m/s]
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
    pressurant_name: str
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
    maximum_wall_temperature: float  # [K]
    thrust_chamber_wall_emissivity: float  # [-]
    hot_gas_emissivity: float  # [-]
    shaft_mechanical_efficiency: float  # [-]

    is_frozen: bool = True
    # Values that override other inputs (one of them is required)
    expansion_ratio: Optional[float] = None  # [-]
    pressure_ratio: Optional[float] = None  # [-]
    exit_pressure_forced: Optional[float] = None  # [Pa]

    # Values that can be estimated or are not necessarily required
    mass_mixture_ratio: Optional[float] = None  # [-]
    specific_impulse_quality_factor: Optional[float] = None  # [-]
    oxidizer_initial_temperature: Optional[float] = None  # [k]
    fuel_initial_temperature: Optional[float] = None  # [K]
    recovery_factor: Optional[float] = None  # [-]
    area_ratio_chamber_throat: Optional[float] = None  # [-]
    chamber_characteristic_length: Optional[float] = None  # [m]
    expansion_ratio_start_cooling: Optional[float] = None
    expansion_ratio_end_cooling: Optional[float] = None  # [-]
    distance_from_throat_end_cooling: Optional[float] = None  # [m]
    distance_from_throat_start_cooling: Optional[float] = None  # [m]
    cooling_section_pressure_drop: Optional[float] = None  # [Pa]
    injector_pressure_drop: Optional[float] = None  # [Pa]
    ambient_pressure: Optional[float] = None  # [Pa]
    fuel_pump_outlet_pressure_forced: Optional[float] = None  # [Pa]
    oxidizer_pump_outlet_pressure_forced: Optional[float] = None  # [Pa]

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
    heat_flow_rate: float = field(init=False, repr=False, default=0)
    minimum_required_coolant_mass_flow: float = field(init=False, repr=False, default=0)
    _fuel_pump_outlet_pressure: float = field(init=False, repr=False, default=None)
    _oxidizer_pump_outlet_pressure: float = field(init=False, repr=False, default=None)
    _expansion_ratio_end_cooling: float = field(init=False, repr=False, default=None)
    _cea_frozen: float = field(init=False, repr=False, default=None)
    _cea_frozenAtThroat: float = field(init=False, repr=False, default=None)
    verbose: bool = True
    iterate: bool = True

    def __post_init__(self):
        self.initialize_cea()
        self.set_initial_values()
        if self._ignore_cooling:
            # self._is_temp_calc_needed = False
            warnings.warn('!!_ignore_cooling flag has been set to True!!')
            self.set_pump_outlet_pressures()
        else:
            self.set_heat_transfer()
            self.set_pump_outlet_pressures()
            self.calc_minimum_required_coolant_mass_flow()
        if self.iterate:
            self.print_verbose_start_iteration_message()
            # self.print_verbose_iteration_message()

            self.iterate_flow()

            self.print_verbose_end_iteration_message()

    def initialize_cea(self):
        # Setting of internal variables
        self._cea_frozen, self._cea_frozenAtThroat = (1, 1) if self.is_frozen else (0, 0)
        if self.mass_mixture_ratio is None:
            self.mass_mixture_ratio = get_mass_mixture_ratio(
                propellant_mix=self.propellant_mix
            )
        self.resolve_expansion_choice()
        self.set_cea()
        if self.expansion_ratio is None:
            self.expansion_ratio = get_expansion_ratio_from_p_ratio(self.pressure_ratio,
                                                                    self.cc_hot_gas_heat_capacity_ratio)

    def resolve_expansion_choice(self):
        if self.exit_pressure_forced is not None:
            if self.pressure_ratio is not None or self.expansion_ratio is not None:
                raise ValueError(
                    'Please provide one and only one of exit_pressure_forced, expansion_ratio, pressure_ratio.'
                )
            self.pressure_ratio = self.combustion_chamber_pressure / self.exit_pressure_forced
        elif not ((self.pressure_ratio is None) ^ (self.expansion_ratio is None)):
            raise ValueError(
                'Neither or both the pressure_ratio and expansion_ratio are given. Provide one and only one.'
            )
        elif self.pressure_ratio is None:
            self.cc_hot_gas_heat_capacity_ratio = self.get_heat_capacity_ratio()
            self.pressure_ratio = get_pressure_ratio_fsolve(self.expansion_ratio,
                                                            self.cc_hot_gas_heat_capacity_ratio)

    def set_initial_values(self):
        """Set missing input values."""
        # Order is important
        if self.specific_impulse_quality_factor is None:
            self.specific_impulse_quality_factor = get_specific_impulse_quality_factor(
                propellant_mix=self.propellant_mix
            )
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
                propellant_mix=self.propellant_mix
            )
        if self.area_ratio_chamber_throat is None:
            self.area_ratio_chamber_throat = get_chamber_throat_area_ratio_estimate(
                throat_area=self.throat_area
            )
        if self.cc_hot_gas_prandtl_number is None:
            self.cc_hot_gas_prandtl_number = get_prandtl_number_estimate(
                heat_capacity_ratio=self.cc_hot_gas_heat_capacity_ratio
            )
        if self.recovery_factor is None:
            self.recovery_factor = get_turbulent_recovery_factor(
                prandtl_number=self.cc_hot_gas_prandtl_number
            )
        if self.cooling_section_pressure_drop is None:
            self.cooling_section_pressure_drop = self.combustion_chamber_pressure * self.cooling_pressure_drop_factor
        if self.injector_pressure_drop is None:
            self.injector_pressure_drop = self.combustion_chamber_pressure * self.injector_pressure_drop_factor

    def set_heat_transfer(self):
        heat_transfer = self.heat_transfer_section
        self.heat_flow_rate = heat_transfer.total_heat_transfer
        self.heat_flux_func = heat_transfer.heat_transfer_func

    def iterate_flow(self):
        raise NotImplementedError

    @property
    def verbose_iteration_name(self):
        raise NotImplementedError

    @property
    def verbose_iteration_actual(self):
        raise NotImplementedError

    @property
    def verbose_iteration_required(self):
        raise NotImplementedError

    def print_verbose_start_iteration_message(self):
        if self.verbose:
            print(f'Start {self.verbose_iteration_name} Iteration')

    def print_verbose_iteration_message(self):
        if self.verbose:
            print(f'Actual:   {self.verbose_iteration_actual:.5f} kg/s\t'
                  f'Required: {self.verbose_iteration_required:.5f} kg/s')

    def print_verbose_end_iteration_message(self):
        if self.verbose:
            print(f'{self.verbose_iteration_name} Set\n')

    @cached_property
    def propellant_mix(self):
        return get_propellant_mixture(fuel_name=self.fuel_name, oxidizer_name=self.oxidizer_name)

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
        # Get the heat_capacity_ratio only, to be able to estimate a pressure ratio, which is required for setting all
        # other CEA values
        kwargs = {key: value for key, value in self.cea_kwargs.items() if key not in ('eps', 'PcOvPe')}
        return get_cea_chamber_dict(**kwargs)['y_cc']

    def set_pump_outlet_pressures(self):
        Merger._warn_pressure = False
        self.calc_pump_outlet_pressures()
        Merger._warn_pressure = True

    def calc_pump_outlet_pressures(self):
        self._fuel_pump_outlet_pressure = self.fuel_pump_expected_pressure if self.fuel_pump_outlet_pressure_forced is None else self.fuel_pump_outlet_pressure_forced
        self._oxidizer_pump_outlet_pressure = self.oxidizer_pump_expected_pressure if self.oxidizer_pump_outlet_pressure_forced is None else self.oxidizer_pump_outlet_pressure_forced

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
    def thrust_coefficient(self):
        return get_thrust_coefficient_from_ideal(ideal_thrust_coefficient=self.ideal_thrust_coefficient,
                                                 chamber_pressure=self.combustion_chamber_pressure,
                                                 exit_pressure=self.exit_pressure,
                                                 expansion_ratio=self.expansion_ratio,
                                                 ambient_pressure=self.ambient_pressure, )

    @cached_property
    def chamber_equivalent_velocity(self):
        return self.thrust_coefficient * self.characteristic_velocity

    @cached_property
    def chamber_specific_impulse(self):
        return self.specific_impulse_quality_factor * self.chamber_equivalent_velocity / constants.g

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
        return self.chamber_mass_flow / (self.mass_mixture_ratio + 1)

    @property
    def chamber_oxidizer_flow(self):
        return (self.mass_mixture_ratio * self.chamber_mass_flow) / (self.mass_mixture_ratio + 1)

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
    def oxidizer_main_flow_state(self):
        return FlowState(propellant_name=self.oxidizer_name,
                         temperature=self.oxidizer_initial_temperature,
                         pressure=self.oxidizer_initial_pressure,
                         mass_flow=self.main_oxidizer_flow,
                         type='oxidizer', )

    @property
    def fuel_main_flow_state(self):
        return FlowState(propellant_name=self.fuel_name,
                         temperature=self.fuel_initial_temperature,
                         pressure=self.fuel_initial_pressure,
                         mass_flow=self.main_fuel_flow,
                         type='fuel', )

    @property
    def oxidizer(self):
        return Propellant(main_flow_state=self.oxidizer_main_flow_state,
                          burn_time=self.burn_time,
                          margin_factor=self.propellant_margin_factor)

    @property
    def fuel(self):
        return Propellant(main_flow_state=self.fuel_main_flow_state,
                          burn_time=self.burn_time,
                          margin_factor=self.propellant_margin_factor)

    @property
    def pressurant_initial_state(self):
        return FlowState(propellant_name=self.pressurant_name,
                         temperature=self.pressurant_initial_temperature,
                         pressure=self.pressurant_initial_pressure,
                         mass_flow=None,
                         type='pressurant')

    @property
    def pressurant(self):
        return Pressurant(oxidizer_volume=self.oxidizer.volume,
                          fuel_volume=self.fuel.volume,
                          fuel_tank_initial_pressure=self.fuel_initial_pressure,
                          oxidizer_tank_initial_pressure=self.oxidizer_initial_pressure,
                          margin_factor=self.pressurant_margin_factor,
                          initial_fluid_state=self.pressurant_initial_state,
                          final_pressure=self.pressurant_final_pressure,
                          propellant_tanks_ullage_factor=self.ullage_volume_factor)

    @property
    def pressurant_tank(self):
        return PressurantTank(structure_material=self.pressurant_tank_material,
                              safety_factor=self.pressurant_tank_safety_factor,
                              pressurant=self.pressurant)

    @property
    def oxidizer_tank(self):
        return Tank(inlet_flow_state=self.oxidizer_main_flow_state,
                    propellant_volume=self.oxidizer.volume,
                    max_acceleration=self.max_acceleration,
                    ullage_factor=self.ullage_volume_factor,
                    pressurant_tank_volume=self.pressurant_tank.volume,
                    structure_material=self.oxidizer_tank_material,
                    safety_factor=self.tanks_structural_factor, )

    @property
    def fuel_tank(self):
        return Tank(inlet_flow_state=self.fuel_main_flow_state,
                    propellant_volume=self.fuel.volume,
                    max_acceleration=self.max_acceleration,
                    ullage_factor=self.ullage_volume_factor,
                    pressurant_tank_volume=None,
                    structure_material=self.fuel_tank_material,
                    safety_factor=self.tanks_structural_factor)

    @property
    def oxidizer_pump(self):
        return Pump(inlet_flow_state=self.oxidizer_tank.outlet_flow_state,
                    expected_outlet_pressure=self.oxidizer_pump_outlet_pressure,
                    efficiency=self.oxidizer_pump_efficiency,
                    specific_power=self.oxidizer_pump_specific_power, )

    @property
    def fuel_pump(self):
        return Pump(inlet_flow_state=self.fuel_tank.outlet_flow_state,
                    expected_outlet_pressure=self.fuel_pump_outlet_pressure,
                    efficiency=self.fuel_pump_efficiency,
                    specific_power=self.fuel_pump_specific_power)

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
                        pressure_drop=self.injector_pressure_drop, )

    @cached_property
    def nozzle(self):
        return Nozzle(throat_area=self.throat_area,
                      expansion_ratio=self.expansion_ratio,
                      area_ratio_chamber_throat=self.area_ratio_chamber_throat,
                      conv_half_angle=self.convergent_half_angle,
                      conv_throat_bend_ratio=self.convergent_throat_bend_ratio,
                      conv_chamber_bend_ratio=self.convergent_chamber_bend_ratio,
                      div_throat_half_angle=self.divergent_throat_half_angle,
                      structure_material=self.nozzle_material,
                      chamber_pressure=self.combustion_chamber_pressure,
                      safety_factor=self.nozzle_safety_factor, )

    @cached_property
    def combustion_chamber(self):
        return CombustionChamber(throat_area=self.throat_area,
                                 combustion_chamber_pressure=self.combustion_chamber_pressure,
                                 area_ratio_chamber_throat=self.area_ratio_chamber_throat,
                                 characteristic_length=self.chamber_characteristic_length,
                                 convergent_volume_estimate=self.nozzle.conv_volume_estimate,
                                 safety_factor=self.combustion_chamber_safety_factor,
                                 structure_material=self.combustion_chamber_material, )

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
                self._expansion_ratio_end_cooling = r_end ** 2 * constants.pi / self.throat_area
                return self.distance_from_throat_end_cooling

    @property
    def min_distance_from_throat_heat_transfer_section(self):
        return self.distance_from_throat_start_cooling

    @property
    def chamber_flow_state(self):
        return ManualFlowState(propellant_name='ChamberGas',
                               temperature=self.combustion_temperature,
                               pressure=self.combustion_chamber_pressure,
                               mass_flow=self.chamber_mass_flow,
                               type='combusted',
                               _molar_mass=self.cc_hot_gas_molar_mass,
                               _specific_heat_capacity=self.cc_hot_gas_specific_heat_capacity,
                               _heat_capacity_ratio=self.cc_hot_gas_heat_capacity_ratio,
                               _dynamic_viscosity=self.cc_hot_gas_dynamic_viscosity,
                               _prandtl_number=self.cc_hot_gas_prandtl_number, )

    @property
    def convective_heat_transfer_args(self):
        return {'combustion_chamber_flow_state': self.chamber_flow_state,
                'maximum_wall_temperature': self.maximum_wall_temperature,
                'hot_gas_emissivity': self.hot_gas_emissivity,
                'thrust_chamber': self.thrust_chamber,
                'recovery_factor': self.recovery_factor,
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
            theoretical_total_convective_heat_transfer=self.theoretical_convective_heat_transfer.total_convective_heat_transfer,
        )

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
        q_tot = self.heat_flow_rate
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
                                     heat_flow_rate=self.heat_flow_rate,
                                     maximum_outlet_temperature=self.maximum_coolant_outlet_temperature,
                                     pressure_drop=self.cooling_section_pressure_drop,
                                     _is_temp_calc_needed=self._is_temp_calc_needed, )

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
    def fuel_pump_power_efficiency(self):
        return self.fuel_pumps_power_required / self.fuel_pump.mass_flow

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
    def mass_kwak(self):
        return self.props_mass + self.tanks_mass + self.pumps_mass + self.pressurant.mass

    @property
    def mass_ratio_kwak(self):
        return 1 / self.inverse_mass_ratio_kwak

    @property
    def inverse_mass_ratio_kwak(self):
        return 1 - self.props_mass / self.mass_kwak

    @property
    def ideal_delta_v_kwak(self):
        return self.overall_specific_impulse * log(self.mass_ratio_kwak) * constants.g

    @property
    def total_thrust_chamber_mass(self):
        return self.thrust_chamber.mass + self.injector.mass + self.cooling_channel_section.mass

    @property
    def chamber_propellant_mass(self):
        return self.chamber_mass_flow * self.burn_time * self.propellant_margin_factor

    @property
    def feed_system_mass(self):
        return self.pumps_mass

    @property
    def cc_prop_group_mass(self):
        return self.chamber_propellant_mass + self.tanks_mass + self.pressurant.mass

    @property
    def energy_source_mass(self):
        raise NotImplementedError

    @property
    def energy_source_ratio(self):
        return self.energy_source_mass / self.burn_time

    @property
    def cc_prop_group_ratio(self):
        return self.cc_prop_group_mass / self.burn_time

    @property
    def power_mass(self):
        return self.feed_system_mass

    @property
    def anti_power_mass(self):
        return self.tanks_plus_pressurant + self.chamber_propellant_mass

    @property
    def combo(self):
        return self.power_mass + self.anti_power_mass

    @property
    def power_ratio(self):
        return self.power_mass / self.burn_time

    @property
    def anti_power_ratio(self):
        return self.anti_power_mass / self.burn_time

    @property
    def initial_mass_ratio(self):
        return self.initial_mass / self.burn_time

    @property
    def feed_system_ratio(self):
        return self.feed_system_mass / self.burn_time

    @property
    def tanks_plus_pressurant(self):
        return self.tanks_mass + self.pressurant.mass

    @property
    def engine_dry_mass(self):
        return self.feed_system_mass + self.total_thrust_chamber_mass

    @property
    def dry_mass(self):
        return self.engine_dry_mass + self.tanks_mass

    @property
    def initial_mass(self):
        return self.final_mass + self.props_mass

    @property
    def final_mass(self):
        return self.dry_mass + self.pressurant.mass

    @property
    def inverse_mass_ratio(self):
        return self.final_mass / self.initial_mass

    @property
    def mass_ratio(self):
        return self.initial_mass / self.final_mass

    def get_payload(self, delta_v: float) -> float:
        e_dv = exp(delta_v / (self.overall_specific_impulse * constants.g))
        m0 = self.initial_mass
        mf = self.final_mass
        return (m0 - mf * e_dv) / (e_dv - 1)

    @property
    def perct(self):
        return (self.total_thrust_chamber_mass + self.feed_system_mass) / self.initial_mass

    @property
    def ideal_delta_v(self):
        return self.overall_specific_impulse * log(1 / self.inverse_mass_ratio) * constants.g

    @property
    def change_in_velocity(self):
        return self.ideal_delta_v

    @property
    def change_in_velocity_km(self):
        return self.change_in_velocity * 1e-3

    def get_payload_delta_v(self, payload_mass):
        mass_ratio = (self.final_mass + payload_mass) / (self.initial_mass + payload_mass)
        return self.overall_specific_impulse * log(1 / mass_ratio) * constants.g

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

    def adjusted_mass_ratio(self, payload=0, factor_0=1.1459, factor_f=1.9188, ):
        m_0_adj = self.initial_mass * factor_0 + payload
        m_f_adj = self.final_mass * factor_f + payload
        return m_0_adj / m_f_adj

    def adjusted_mass_ratio2(self, payload=0, factor_0=0.04316, factor_f=0.1027, ):
        m_0_adj = self.initial_mass * factor_0
        m_f_adj = self.initial_mass * factor_f
        return (self.initial_mass + m_0_adj + m_f_adj + payload) / (self.final_mass + m_f_adj + payload)

    @property
    def payload_mass_ratio(self):
        if self.mass_u is None:
            ValueError('Payload mass not given, impossible to calculate mass ratio with payload mass')
        return (self.mass + self.mass_u - self.props_mass) / (self.mass + self.mass_u)

    @property
    def chamber_thrust(self):
        """For consistency in naming with OpenCycles"""
        return self.thrust

    @property
    def components_list(self):
        return [
            'Fuel',
            'Oxidizer',
            'Pressurant',
            'Fuel Tank',
            'Oxidizer Tank',
            'Pressurant Tank',
            'Fuel Pump',
            'Oxidizer Pump',
            'Injector',
            'Combustion Chamber',
            'Nozzle',
            'Cooling Channel Section',
        ]

    @property
    def aggregate_masses(self):
        return {
            'Initial': self.initial_mass,
            'Final': self.final_mass,
            'Feed System': self.feed_system_mass,
            'Propellant': self.props_mass,
            'Chamber Propellant': self.chamber_propellant_mass,
            'Tanks': self.tanks_mass,
            'Pumps': self.pumps_mass,
            'Thrust Chamber': self.total_thrust_chamber_mass,
            'Engine Dry': self.engine_dry_mass,
            'Dry': self.dry_mass,
            'Power': self.power_mass,

        }

    @property
    def components_masses(self):
        return {
            name: getattr(self, format_fancy_name(name)).mass
            for name in self.components_list
        }

    @property
    def combined_info(self):
        return self.components_masses | self.aggregate_masses | {
            'Specific Impulse [s]': self.overall_specific_impulse,
            'Mass Ratio [-]': self.mass_ratio,
            'Velocity Change [m/s]': self.change_in_velocity,
        }

    def print_masses(self, decimals: int = 2):
        print('\nComponent Masses [kg]:')
        for key, value in self.components_masses.items():
            print(f'{key:<25}: {value:>10.{decimals}f}')
        print('\nAggregation Masses [kg]:')
        for key, value in self.aggregate_masses.items():
            print(f'{key:<25}: {value:>10.{decimals}f}')
