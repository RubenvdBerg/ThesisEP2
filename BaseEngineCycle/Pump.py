from dataclasses import dataclass

from BaseEngineCycle.Propellant import Propellant
from BaseEngineCycle.FlowComponent import FlowComponent


@dataclass
class Pump(FlowComponent):
    pressure_increase: float  # [Pa]
    efficiency: float  # [-]
    specific_power: float  # [W/kg]

    @property
    def volumetric_flow_rate(self):
        return self.mass_flow / self.inlet_flow_state.density

    @property
    def power_required(self):
        return self.volumetric_flow_rate * self.pressure_increase / self.efficiency

    @property
    def mass(self):
        return self.power_required / self.specific_power

    @property
    def pressure_change(self):
        return self.pressure_increase

# @dataclass
# class Pump:
#     propellant: Propellant
#     mass_flow: float  # [kg/s]
#     pressure_increase: float  # [Pa]
#     efficiency: float  # [-]
#     specific_power: float  # [W/kg]
#     inlet_pressure: float  # [Pa]
#
#     @property
#     def volumetric_flow_rate(self):
#         return self.mass_flow / self.propellant.density
#
#     @property
#     def power_required(self):
#         return self.volumetric_flow_rate * self.pressure_increase / self.efficiency
#
#     @property
#     def mass(self):
#         return self.power_required / self.specific_power
#
#     @property
#     def outlet_pressure(self):
#         return self.inlet_pressure + self.pressure_increase
