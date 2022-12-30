from dataclasses import dataclass
from typing import Optional
import CoolProp.CoolProp as CoolProp
from EngineComponents.Abstract.FlowComponent import FlowComponent


@dataclass
class CoolingChannelSection(FlowComponent):
    total_heat_transfer: float = 0  # [W]
    maximum_outlet_temperature: float = 0  # [K]
    pressure_drop: Optional[float] = None  # [Pa]
    combustion_chamber_pressure: Optional[float] = None  # [Pa]
    pressure_drop_factor: Optional[float] = None  # [-]
    verbose: bool = True
    _is_temp_calc_needed: bool = True

    def __post_init__(self):
        self.coolprop_name = self.inlet_flow_state.coolprop_name
        CoolProp.set_reference_state(self.coolprop_name, 'NBP')

    @property
    def pressure_change(self):
        if self.pressure_drop is None and self.combustion_chamber_pressure is None:
            raise ValueError(
                'One of [pressure_drop] and [pressure_drop_factor](to estimate the pressure_drop) must be given'
            )
        if self.pressure_drop is None:
            # Humble 1995 p.209 suggest pressure drop to be 10% - 20% of chamber pressure
            return -self.combustion_chamber_pressure * self.pressure_drop_factor
        else:
            return -self.pressure_drop

    @property
    def increase_mass_specific_enthalpy(self):
        return self.total_heat_transfer / self.inlet_flow_state.mass_flow

    @property
    def outlet_mass_specific_enthalpy(self):
        return self.inlet_flow_state.mass_specific_enthalpy + self.increase_mass_specific_enthalpy

    @property
    def ideal_outlet_temperature(self):
        try:
            return CoolProp.PropsSI('T',
                                    'H', self.outlet_mass_specific_enthalpy,
                                    'P', self.inlet_pressure,
                                    self.coolprop_name)
        except ValueError:
            raise ValueError('The coolant reached an unacceptable state (see previous error in the stack)')

    @property
    def temperature_change(self):
        if self._is_temp_calc_needed:
            return self.ideal_outlet_temperature - self.inlet_temperature
        else:
            return self.maximum_outlet_temperature - self.inlet_temperature