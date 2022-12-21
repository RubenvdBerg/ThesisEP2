import warnings
from dataclasses import dataclass, field
from math import log
from typing import Optional
from EngineCycles.BaseEngineCycle.FlowComponent import FlowComponent
from EngineCycles.BaseEngineCycle.FlowState import FlowState


@dataclass
class ElectricComponent:
    electric_energy_efficiency: float  # -
    specific_power: float  # W/kg
    output_power: float  # W

    @property
    def input_power(self):
        return self.output_power / self.electric_energy_efficiency

    @property
    def mass(self):
        return self.output_power / self.specific_power


@dataclass
class ElectricMotor(ElectricComponent):
    electric_heat_loss_factor: float 
    oxidizer_pump_inlet_flow_state: FlowState
    oxidizer_leakage_factor: float
    magnet_temp_limit: float

    @property
    def power_heat_loss(self):
        return self.output_power * .015

    def calc_cooling(self):
        cp_ox = self.oxidizer_pump_inlet_flow_state.specific_heat_capacity
        m_ox_leak = self.oxidizer_pump_inlet_flow_state.mass_flow * self.oxidizer_leakage_factor
        deltaT = (self.magnet_temp_limit - 50) - self.oxidizer_pump_inlet_flow_state.temperature
        required_mass_flow = self.power_heat_loss / (cp_ox * deltaT)
        if required_mass_flow > m_ox_leak:
            warnings.warn('The expected cooling due to oxidizer leakage through the electric motor is too low. The motor magnets will be too hot and possibly demagnetize.')
            print(f'-------------------------------------------------{required_mass_flow} required, {m_ox_leak} expected')


@dataclass
class Inverter(ElectricComponent):
    pass


@dataclass
class Battery(ElectricComponent):
    specific_energy: float  # J/kg
    battery_packing_factor: float  # -
    burn_time: float  # s

    # Overwrite ElectricComponent attribute
    electric_energy_efficiency: Optional[float] = field(init=False, default=None)

    def __post_init__(self):
        self.electric_energy_efficiency = self.eta_e

    @property
    def eta_e(self):
        return 0.093 * log(self.burn_time) + 0.3301

    @property
    def total_energy(self):
        return self.input_power * self.burn_time

    @property
    def energy_heat_loss(self):
        return self.total_energy * (1 - self.eta_e)

    @property
    def power_heat_loss(self):
        return self.input_power * (1 - self.eta_e)

    @property
    def coolant_flow_required(self):
        return (1 - self.eta_e) * self.total_energy / (self.fuel_specific_heat * self.coolant_allowable_temperature_change * self.burn_time)

    @property  # Overwrite ElectricComponent mass
    def mass(self):
        return self.battery_packing_factor * self.output_power * max(1 / self.specific_power, self.burn_time / (self.specific_energy * self.eta_e))


@dataclass
class BatteryCooler(FlowComponent):
    """Component that adjusts its outlet flow to be equal to the coolant flow required. Used to iterate until flows
    through pump and battery are matching in the EP Cycle"""
    power_heat_loss: float = 0  # [W]
    outlet_pressure_required: float = 0  # [Pa]
    coolant_allowable_temperature_change: float = 0  # [K]
    coolant_specific_heat_capacity: Optional[float] = None  # [J/(kgK)]

    def __post_init__(self):
        if self.coolant_specific_heat_capacity is None:
            self.coolant_specific_heat_capacity = self.inlet_flow_state.specific_heat_capacity

    @property
    def mass_flow_change(self):
        return self.coolant_flow_required - self.inlet_mass_flow

    @property
    def pressure_change(self):
        return self.outlet_pressure_required - self.inlet_pressure

    @property
    def coolant_flow_required(self):
        return self.power_heat_loss / (self.coolant_specific_heat_capacity * self.coolant_allowable_temperature_change)
