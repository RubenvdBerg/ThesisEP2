from dataclasses import dataclass, field
from typing import Optional
import CoolProp.CoolProp as CoolProp
from BaseEngineCycle.FlowComponent import FlowComponent


@dataclass
class CoolingChannelSection(FlowComponent):
    total_heat_transfer: float  # [W]
    pressure_drop: Optional[float] = None  # [Pa]
    combustion_chamber_pressure: Optional[float] = None  # [Pa]
    _pressure_drop_ratio: float = field(init=False, default=.15)  # [-]
    verbose: bool = True

    def __post_init__(self):
        self.coolprop_name = self.inlet_flow_state.coolprop_name
        CoolProp.set_reference_state(self.coolprop_name, 'NBP')

    @property
    def pressure_change(self):
        if self.pressure_drop is None and self.combustion_chamber_pressure is None:
            raise ValueError(
                'One of [pressure_drop] and [combustion_chamber_pressure](to estimate the pressure_drop) must be given'
            )
        if self.pressure_drop is None:
            # Humble 1995 p.209 suggest pressure drop to be 10% - 20% of chamber pressure
            return -self.combustion_chamber_pressure * self._pressure_drop_ratio
        else:
            return -self.pressure_drop

    @property
    def default_inlet_temperature(self):
        _temp_in_dict = {'Hydrogen': 20.25, 'Oxygen': 90.15, 'n-Dodecane': 293.15, 'Methane': 111.66}
        return _temp_in_dict[self.coolprop_name]

    @property
    def increase_mass_specific_enthalpy(self):
        return self.total_heat_transfer / self.inlet_flow_state.mass_flow

    @property
    def outlet_mass_specific_enthalpy(self):
        return self.inlet_flow_state.mass_specific_enthalpy + self.increase_mass_specific_enthalpy

    @property
    def ideal_outlet_temperature(self):
        return CoolProp.PropsSI('T', 'H', self.outlet_mass_specific_enthalpy, 'P', self.outlet_pressure,
                                self.coolprop_name)

    @property
    def temperature_change(self):
        return self.ideal_outlet_temperature - self.inlet_temperature
