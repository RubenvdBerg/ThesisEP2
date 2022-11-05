from dataclasses import dataclass
from typing import Optional
from EngineCycles.Functions.BaseFunctions import only_one_none
from EngineCycles.BaseEngineCycle.FlowComponent import FlowComponent
from functools import cached_property


@dataclass
class BaseTurbine:
    """Finds the missing value based on the basic turbine power equation
    (eq. 10-18, Rocket Propulsion Elements 9th Ed., Sutton)
    """
    pump_power_required: float = 0  # [W]
    efficiency: float = 0  # [-]
    specific_heat_capacity: float = 0  # [J/kg]
    heat_capacity_ratio: float = 0  # [-]

    # Two out of these three must be given
    pressure_ratio: Optional[float] = None  # [-]
    mass_flow: Optional[float] = None  # [kg/s]
    inlet_temperature: Optional[float] = None  # [K]

    def __post_init__(self):
        if not only_one_none(self.pressure_ratio, self.mass_flow, self.inlet_temperature):
            raise ValueError('Exactly two out of pressure_ratio, mass_flow, and inlet_temperature must be provided')

    @property
    def mass_flow_required(self):
        if self.mass_flow is None:
            return (self.pump_power_required
                    / (self.efficiency * self.specific_heat_capacity * self.inlet_temperature
                       * (1 - self.pressure_ratio ** ((1 - self.heat_capacity_ratio) / self.heat_capacity_ratio))))
        else:
            raise ValueError(self.error_message('mass_flow'))

    @property
    def inlet_temperature_required(self):
        if self.inlet_temperature is None:
            return (self.pump_power_required
                    / (self.efficiency * self.specific_heat_capacity * self.mass_flow
                       * (1 - self.pressure_ratio ** ((1 - self.heat_capacity_ratio) / self.heat_capacity_ratio))))
        else:
            raise ValueError(self.error_message('inlet_temperature'))

    @property
    def pressure_ratio_required(self):
        if self.inlet_temperature is None:
            return (self.pump_power_required
                    / (self.efficiency * self.specific_heat_capacity * self.mass_flow
                       * (1 - self.pressure_ratio ** ((1 - self.heat_capacity_ratio) / self.heat_capacity_ratio))))
        else:
            raise ValueError(self.error_message('pressure_ratio'))

    def error_message(self, attribute_name):
        return f'Turbine {attribute_name} is already given [{getattr(self, attribute_name):}], also estimating ' \
               f'{attribute_name}_required makes no sense'


@dataclass
class Turbine(FlowComponent):
    pump_power_required: float = 0  # [W]
    efficiency: float = 0  # [-]
    pressure_ratio: float = 0  # [-]
    gas_heat_capacity_ratio: Optional[float] = None  # [-]
    gas_specific_heat_capacity: Optional[float] = None  # [J/kg]

    @cached_property
    def drive_gas_heat_capacity_ratio(self):
        if self.gas_heat_capacity_ratio is None:
            return self.inlet_flow_state.heat_capacity_ratio
        else:
            return self.gas_heat_capacity_ratio

    @cached_property
    def drive_gas_specific_heat_capacity(self):
        if self.gas_specific_heat_capacity is None:
            return self.inlet_flow_state.specific_heat_capacity
        else:
            return self.gas_specific_heat_capacity

    @property
    def mass_flow_required(self):
        return BaseTurbine(
            pump_power_required=self.pump_power_required,
            efficiency=self.efficiency,
            specific_heat_capacity=self.drive_gas_specific_heat_capacity,
            heat_capacity_ratio=self.drive_gas_heat_capacity_ratio,
            pressure_ratio=self.pressure_ratio,
            inlet_temperature=self.inlet_temperature
        ).mass_flow_required

    @property
    def temperature_ratio(self):
        return self.pressure_ratio ** ((self.drive_gas_heat_capacity_ratio - 1) / self.drive_gas_heat_capacity_ratio)

    @property
    def temperature_change(self):
        return self.inlet_temperature / self.temperature_ratio - self.inlet_temperature

    @property
    def pressure_change(self):
        return self.inlet_pressure / self.pressure_ratio - self.inlet_pressure