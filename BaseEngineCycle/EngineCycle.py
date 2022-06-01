from dataclasses import dataclass
from functools import cached_property
from math import sqrt, log
from typing import Optional

from scipy import constants as constants

from BaseEngineCycle.BaseFunctions import get_propellant_mix_name
from BaseEngineCycle.CombustionChamber import CombustionChamber, Injector
from BaseEngineCycle.Cooling import Coolant, CoolingChannels, HeatExchanger
from BaseEngineCycle.Nozzle import BellNozzle
from BaseEngineCycle.Pressurant import Pressurant, PressurantTank
from BaseEngineCycle.PropellantTank import Propellant, Tank
from BaseEngineCycle.ThrustChamber import ThrustChamber
from BaseEngineCycle.TurboPump import Pump
from cea import get_cea_values_dict
from irt import get_kerckhove, get_expansion_ratio


@dataclass
class EngineCycle:
    thrust: float  # [N]
    burn_time: float  # [s]
    combustion_chamber_pressure: float  # [Pa]
    exit_pressure: float  # [Pa]
    max_acceleration: float  # [m/s]
    mass_mixture_ratio: float  # [-]
    fuel_initial_pressure: float  # [Pa]
    fuel_pump_pressure_factor: float  # [Pa]
    fuel_pump_specific_power: float  # [W]
    fuel_pump_efficiency: float  # [-]
    fuel_density: float  # [kg/m3]
    oxidizer_initial_pressure: float  # [Pa]
    oxidizer_pump_pressure_factor: float  # [Pa]
    oxidizer_pump_specific_power: float  # [W]
    oxidizer_pump_efficiency: float  # [-]
    oxidizer_density: float  # [kg/m3]
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
    kwak_fix_cycle_type: str
    oxidizer_name: str
    fuel_name: str
    is_frozen: bool

    # Values that override other inputs
    area_ratio_forced: Optional[float] = None

    # Values that can be estimated
    recovery_factor: Optional[float] = None  # [-]
    chamber_throat_area_ratio: Optional[float] = None  # [-]

    # Values that can be estimated by CEA
    characteristic_velocity: Optional[float] = None  # [m/s]
    thrust_coefficient: Optional[float] = None  # [-]
    combustion_temperature: Optional[float] = None  # [K]
    cc_hot_gas_molar_mass: Optional[float] = None  # [kg/mol]
    cc_hot_gas_heat_capacity_ratio: Optional[float] = None  # [-]
    cc_hot_gas_dynamic_viscosity: Optional[float] = None  # [Pa*s]
    cc_hot_gas_prandtl_number: Optional[float] = None  # [-]
    cc_hot_gas_specific_heat_capacity: Optional[float] = None  # [J/(kg*K)]
    kwak_fix: bool = False

    def __post_init__(self):
        assert self.kwak_fix_cycle_type in ['gg', 'ep']
        self.set_cea()

    @cached_property
    def cea_dict(self):
        return {'characteristic_velocity': 'c_star',
                'thrust_coefficient': 'C_F',
                'combustion_temperature': 'T_C',
                'cc_hot_gas_molar_mass': 'mm_cc',
                'cc_hot_gas_heat_capacity_ratio': 'y_cc',
                'cc_hot_gas_dynamic_viscosity': 'mu_cc',
                'cc_hot_gas_prandtl_number': 'pr_cc',
                'cc_hot_gas_specific_heat_capacity': 'cp_cc'}

    def get_cea(self):
        return get_cea_values_dict(
            chamber_pressure=self.combustion_chamber_pressure, mixture_ratio=self.mass_mixture_ratio,
            exit_pressure=self.exit_pressure, fuel_name=self.fuel_name, ox_name=self.oxidizer_name,
            isfrozen=self.is_frozen, area_ratio=self.area_ratio_forced
        )

    def set_cea(self):
        cea_values = self.get_cea()
        # Checking if value is given, if not: assign value found by CEA. Despite this check for each attribute, it is
        # recommended to either provide none of the CEA properties or all of them
        cea_attributes = self.cea_dict.keys()
        for attribute in cea_attributes:
            cea_name = self.cea_dict[attribute]
            if getattr(self, attribute) is None:
                setattr(self, attribute, cea_values[cea_name])

    def reiterate(self):
        cea_values = self.get_cea()
        for attribute in self.cea_dict.keys():
            cea_name = self.cea_dict[attribute]
            setattr(self, attribute, cea_values[cea_name])

    @property
    def cstar_cf(self):
        return self.characteristic_velocity, self.thrust_coefficient

    @property
    def propellant_mix_name(self):
        return get_propellant_mix_name(fuel_name=self.fuel_name, oxidizer_name=self.oxidizer_name)

    @property
    def mass_flow(self):
        return self.thrust / (self.characteristic_velocity * self.thrust_coefficient)

    @property
    def throat_area(self):
        return self.mass_flow * sqrt(constants.R / self.cc_hot_gas_molar_mass * self.combustion_temperature) / (
                get_kerckhove(self.cc_hot_gas_heat_capacity_ratio) * self.combustion_chamber_pressure)

    @property
    def area_ratio(self):
        return get_expansion_ratio(self.pressure_ratio, self.cc_hot_gas_heat_capacity_ratio)

    @property
    def exit_area(self):
        return self.throat_area * self.area_ratio

    @property
    def pressure_ratio(self):
        return self.combustion_chamber_pressure / self.exit_pressure

    @property
    def simple_specific_impulse(self):
        return self.thrust / self.mass_flow / g

    @property
    def fuel_flow(self):
        return 1 / (self.mass_mixture_ratio + 1) * self.mass_flow

    @property
    def oxidizer_flow(self):
        return self.mass_mixture_ratio / (self.mass_mixture_ratio + 1) * self.mass_flow

    @property
    def delta_p_oxidizer_pump(self):
        return self.combustion_chamber_pressure * self.oxidizer_pump_pressure_factor - self.oxidizer_initial_pressure

    @property
    def delta_p_fuel_pump(self):
        return self.combustion_chamber_pressure * self.fuel_pump_pressure_factor - self.fuel_initial_pressure

    @property
    def oxidizer(self):
        return Propellant(mass_flow=self.oxidizer_flow,
                          burn_time=self.burn_time,
                          density=self.oxidizer_density,
                          type='oxidizer',
                          margin_factor=self.propellant_margin_factor)

    @property
    def fuel(self):
        return Propellant(mass_flow=self.fuel_flow,
                          burn_time=self.burn_time,
                          density=self.fuel_density,
                          type='fuel',
                          margin_factor=self.propellant_margin_factor)

    @property
    def pressurant(self):
        return Pressurant(fuel=self.fuel,
                          oxidizer=self.oxidizer,
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
        return Tank(max_acceleration=self.max_acceleration,
                    ullage_factor=self.ullage_volume_factor,
                    propellant=self.oxidizer,
                    pressurant_tank_volume=self.pressurant_tank.volume,
                    initial_pressure=self.oxidizer_initial_pressure,
                    material_density=self.tanks_material_density,
                    yield_strength=self.tanks_yield_strength,
                    safety_factor=self.tanks_structural_factor,
                    kwak_fix_cycle_type=self.kwak_fix_cycle_type,
                    kwak_fix=self.kwak_fix)

    @property
    def fuel_tank(self):
        return Tank(max_acceleration=self.max_acceleration,
                    ullage_factor=self.ullage_volume_factor,
                    propellant=self.fuel,
                    pressurant_tank_volume=None,
                    initial_pressure=self.fuel_initial_pressure,
                    material_density=self.tanks_material_density,
                    yield_strength=self.tanks_yield_strength,
                    safety_factor=self.tanks_structural_factor,
                    kwak_fix_cycle_type=self.kwak_fix_cycle_type,
                    kwak_fix=self.kwak_fix)

    @property
    def oxidizer_pump(self):
        return Pump(propellant=self.oxidizer,
                    pressure_increase=self.delta_p_oxidizer_pump,
                    efficiency=self.oxidizer_pump_efficiency,
                    specific_power=self.oxidizer_pump_specific_power,
                    mass_flow=self.oxidizer_flow)

    @property
    def fuel_pump(self):
        return Pump(propellant=self.fuel,
                    pressure_increase=self.delta_p_fuel_pump,
                    efficiency=self.fuel_pump_efficiency,
                    specific_power=self.fuel_pump_specific_power,
                    mass_flow=self.fuel_flow)

    @property
    def combustion_chamber(self):
        return CombustionChamber(throat_area=self.throat_area,
                                 combustion_chamber_pressure=self.combustion_chamber_pressure,
                                 propellant_mix=self.propellant_mix_name,
                                 area_ratio_chamber_throat=self.chamber_throat_area_ratio,
                                 safety_factor=self.combustion_chamber_safety_factor,
                                 yield_strength=self.combustion_chamber_yield_strength,
                                 material_density=self.combustion_chamber_material_density)

    @property
    def injector(self):
        return Injector(combustion_chamber_pressure=self.combustion_chamber_pressure,
                        combustion_chamber_area=self.combustion_chamber.area,
                        material_density=self.injector_material_density,
                        safety_factor=self.injector_safety_factor,
                        yield_strength=self.injector_yield_strength,
                        propellant_is_gas=self.injector_propellant_is_gas)

    @property
    def nozzle(self):
        return BellNozzle(throat_area=self.throat_area,
                          area_ratio=self.area_ratio,
                          chamber_radius=self.combustion_chamber.radius,
                          conv_half_angle=self.convergent_half_angle,
                          conv_throat_bend_ratio=self.convergent_throat_bend_ratio,
                          conv_chamber_bend_ratio=self.convergent_chamber_bend_ratio,
                          div_throat_half_angle=self.divergent_throat_half_angle,
                          div_exit_half_angle=self.divergent_exit_half_angle)

    @property
    def thrust_chamber(self):
        return ThrustChamber(nozzle=self.nozzle,
                             chamber=self.combustion_chamber,
                             injector=self.injector,
                             heat_capacity_ratio=self.cc_hot_gas_heat_capacity_ratio)

    @property
    def coolant(self):
        return Coolant(heat_capacity_liquid=1,
                       heat_capacity_gas=1,
                       heat_of_vaporization=1,
                       molar_mass=1,
                       boiling_temperature_1_bar=1)

    @property
    def cooling_channels(self):
        return CoolingChannels(coolant=self.coolant,
                               total_heat_transfer=self.heat_exchanger.total_heat_transfer,
                               outlet_pressure=self.combustion_chamber_pressure,
                               inlet_temperature=1,
                               mass_flow=1)

    @property
    def heat_exchanger(self):
        return HeatExchanger(combustion_temperature=self.combustion_temperature,
                             combustion_chamber_pressure=self.combustion_chamber_pressure,
                             dynamic_viscosity=self.cc_hot_gas_dynamic_viscosity,
                             specific_heat_capacity=self.cc_hot_gas_specific_heat_capacity,
                             mass_flow=self.mass_flow,
                             maximum_wall_temperature=self.maximum_wall_temperature,
                             thrust_chamber_wall_emissivity=self.thrust_chamber_wall_emissivity,
                             hot_gas_emissivity=self.hot_gas_emissivity,
                             heat_capacity_ratio=self.cc_hot_gas_heat_capacity_ratio,
                             convective_coefficient_mode=self.convective_coefficient_mode,
                             thrust_chamber=self.thrust_chamber,
                             recovery_factor=self.recovery_factor,
                             prandtl_number=self.cc_hot_gas_prandtl_number)

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
    def payload_delta_v(self):
        payload = 10  # [kg]
        return self.simple_specific_impulse * log(self.mass + payload / (self.mass - self.props_mass + payload)) * g

    @property
    def payload_mass_ratio(self):
        if self.mass_u is None:
            ValueError('Payload mass not given, impossible to calculate mass ratio with payload mass')
        return (self.mass + self.mass_u - self.props_mass) / (self.mass + self.mass_u)
