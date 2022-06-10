from dataclasses import dataclass

from BaseEngineCycle.Propellant import Propellant


@dataclass
class Pump:
    propellant: Propellant
    mass_flow: float  # [kg/s]
    pressure_increase: float  # [Pa]
    efficiency: float  # [-]
    specific_power: float  # [W/kg]

    @property
    def volumetric_flow_rate(self):
        return self.mass_flow / self.propellant.density

    @property
    def power_required(self):
        return self.volumetric_flow_rate * self.pressure_increase / self.efficiency

    @property
    def mass(self):
        return self.power_required / self.specific_power


