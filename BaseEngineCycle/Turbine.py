from dataclasses import dataclass
from typing import Optional

from BaseEngineCycle.BaseFunctions import only_one_none


@dataclass
class Turbine:
    pump_power_required: float  # [W]
    efficiency: float  # [-]
    specific_heat_capacity: float  # [J/(kg*K)]
    heat_capacity_ratio: float  # [-]
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
            raise ValueError(self.error_message('inlet_temperature'))

    def error_message(self, attribute_name):
        return f'Turbine {attribute_name} is already given [{getattr(self, attribute_name):}], also getting ' \
               f'{attribute_name}_required makes no sense'