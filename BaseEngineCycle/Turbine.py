from dataclasses import dataclass
from typing import Optional

from BaseEngineCycle.BaseFunctions import only_one_none
from CoolProp.CoolProp import PropsSI

@dataclass
class Turbine:
    pump_power_required: float  # [W]
    efficiency: float  # [-]

    # Can be estimated if propellant state is given
    specific_heat_capacity: Optional[float] = None  # [J/(kg*K)]
    heat_capacity_ratio: Optional[float] = None  # [-]

    # Two out of these three must be given
    pressure_ratio: Optional[float] = None  # [-]
    mass_flow: Optional[float] = None  # [kg/s]
    inlet_temperature: Optional[float] = None  # [K]

    # Required for turbine drive gas properties estimation
    inlet_pressure: Optional[float] = None  # [Pa]
    coolprop_name: Optional[str] = None  # CoolPropName of propellant driving the turbine

    def __post_init__(self):
        if not only_one_none(self.pressure_ratio, self.mass_flow, self.inlet_temperature):
            raise ValueError('Exactly two out of pressure_ratio, mass_flow, and inlet_temperature must be provided')

        # Drive gas properties estimation
        coolprop_state_args = ('T', self.inlet_temperature, 'P', self.inlet_pressure, self.coolprop_name)
        if self.specific_heat_capacity is None:
            self.specific_heat_capacity = PropsSI('C', *coolprop_state_args)
        if self.heat_capacity_ratio is None:
            specific_heat_capacity_const_volume = PropsSI('O', *coolprop_state_args)
            self.heat_capacity_ratio = self.specific_heat_capacity / specific_heat_capacity_const_volume

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
        
    @property
    def outlet_temperature(self):
        y = self.heat_capacity_ratio
        return self.inlet_temperature * self.pressure_ratio**((y-1)/y)

    @property
    def outlet_pressure(self):
        return self.pressure_ratio * self.inlet_pressure

    def error_message(self, attribute_name):
        return f'Turbine {attribute_name} is already given [{getattr(self, attribute_name):}], also estimating ' \
               f'{attribute_name}_required makes no sense'