from dataclasses import dataclass
from typing import Optional
from EngineComponents.Abstract.FlowComponent import FlowComponent
import CoolProp.CoolProp as CoolProp


@dataclass
class Pump(FlowComponent):
    expected_outlet_pressure: float = 0  # [Pa]
    efficiency: float = 0  # [-]
    specific_power: float = 0  # [W/kg]

    def __post_init__(self):
        if self.inlet_flow_state.phase != 'liquid':
            raise ValueError('Pump inlet flow is not liquid')

    @property
    def propellant_density(self):
        return self.inlet_flow_state.density

    @property
    def volumetric_flow_rate(self):
        return self.inlet_mass_flow / self.propellant_density

    @property
    def power_required(self):
        if self.pressure_change < 0:
            raise ValueError(
                'Negative pressure change required over pump, increase combustion pressure or decrease tank pressure')
        return self.volumetric_flow_rate * self.pressure_change / self.efficiency

    @property
    def mass(self):
        return self.power_required / self.specific_power

    @property
    def pressure_change(self):
        return self.expected_outlet_pressure - self.inlet_pressure

    @property
    def enthalpy_change(self):
        return self.pressure_change / (self.propellant_density * self.efficiency)

    @property
    def outlet_enthalpy(self):
        return self.inlet_flow_state.mass_specific_enthalpy + self.enthalpy_change

    @property
    def estimated_outlet_temperature(self):
        return CoolProp.PropsSI('T',
                                'H', self.outlet_enthalpy,
                                'P', self.expected_outlet_pressure,
                                self.inlet_flow_state.coolprop_name)

    @property
    def temperature_change(self):
        return self.estimated_outlet_temperature - self.inlet_temperature