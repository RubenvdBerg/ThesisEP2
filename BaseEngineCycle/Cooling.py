from dataclasses import dataclass, field
from typing import Optional
import CoolProp.CoolProp as CoolProp


@dataclass
class CoolingChannelSection:
    coolprop_name: str
    total_heat_transfer: float  # [W]
    outlet_pressure: float  # [Pa]
    mass_flow: float  # [kg/s]
    pressure_drop: Optional[float] = None  # [Pa]
    inlet_temperature: Optional[float] = None  # [K]
    _pressure_drop_ratio: float = field(init=False, default=.15)  # [-]
    _outlet_temp_estimate: float = field(init=False, default=200)  # [K]
    verbose: bool = True

    def __post_init__(self):
        CoolProp.set_reference_state(self.coolprop_name, 'NBP')
        # Set optional variables
        if self.inlet_temperature is None:
            self.inlet_temperature = self.default_inlet_temperature
        if self.pressure_drop is None:
            # Humble 1995 p.209 suggest pressure drop to be 10% - 20% of chamber/outlet pressure
            self.pressure_drop = self.outlet_pressure * self._pressure_drop_ratio

    @property
    def default_inlet_temperature(self):
        _temp_in_dict = {'Hydrogen': 20.25, 'Oxygen': 90.15, 'n-Dodecane': 293.15, 'Methane': 111.66}
        return _temp_in_dict[self.coolprop_name]

    @property
    def inlet_pressure(self):
        return self.outlet_pressure + self.pressure_drop

    @property
    def inlet_mass_specific_enthalpy(self):
        return CoolProp.PropsSI('H', 'T', self.inlet_temperature, 'P', self.inlet_pressure, self.coolprop_name)

    @property
    def increase_mass_specific_enthalpy(self):
        return self.total_heat_transfer / self.mass_flow

    @property
    def outlet_mass_specific_enthalpy(self):
        return self.inlet_mass_specific_enthalpy + self.increase_mass_specific_enthalpy

    @property
    def outlet_temperature(self):
        return CoolProp.PropsSI('T', 'H', self.outlet_mass_specific_enthalpy, 'P', self.outlet_pressure,
                                self.coolprop_name)
