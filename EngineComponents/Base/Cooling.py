from dataclasses import dataclass
import CoolProp.CoolProp as CoolProp
from EngineComponents.Abstract.FlowComponent import FlowComponent


@dataclass
class CoolingChannelSection(FlowComponent):
    heat_flow_rate: float = 0  # [W]
    maximum_outlet_temperature: float = 0  # [K]
    pressure_drop: float = 0  # [Pa]
    _is_temp_calc_needed: bool = True

    @property
    def pressure_change(self):
        return -self.pressure_drop

    @property
    def increase_mass_specific_enthalpy(self):
        return self.heat_flow_rate / self.inlet_flow_state.mass_flow

    @property
    def outlet_mass_specific_enthalpy(self):
        return self.inlet_flow_state.mass_specific_enthalpy + self.increase_mass_specific_enthalpy

    @property
    def calculated_outlet_temperature(self):
        try:
            return CoolProp.PropsSI('T',
                                    'H', self.outlet_mass_specific_enthalpy,
                                    'P', self.outlet_pressure,
                                    self.inlet_flow_state.coolprop_name)
        except ValueError:
            raise ValueError('The coolant reached an unacceptable state (see previous error in the stack)')

    @property
    def temperature_change(self):
        if self._is_temp_calc_needed:
            return self.calculated_outlet_temperature - self.inlet_temperature
        else:
            return self.maximum_outlet_temperature - self.inlet_temperature

    @property
    def mass(self):
        return 0
