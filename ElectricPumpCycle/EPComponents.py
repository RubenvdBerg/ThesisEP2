from dataclasses import dataclass, field
from math import log
from typing import Optional
from BaseEngineCycle.FlowComponent import FlowComponent


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
    pass


@dataclass
class Inverter(ElectricComponent):
    pass


@dataclass
class Battery(ElectricComponent):
    specific_energy: float  # J/kg
    battery_packing_factor: float  # -
    burn_time: float  # s
    fuel_specific_heat: float  # J/(kg*K)
    coolant_allowable_temperature_change: float  # K
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
    def heat_loss(self):
        return self.total_energy * (1 - self.eta_e)

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
    coolant_flow_required: float = 0  # [kg/s]
    outlet_pressure_required: float = 0  # [Pa]

    @property
    def mass_flow_change(self):
        return self.coolant_flow_required - self.inlet_mass_flow

    @property
    def pressure_change(self):
        return self.outlet_pressure_required - self.inlet_pressure
