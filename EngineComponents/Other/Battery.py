from dataclasses import dataclass, field
from math import log
from typing import Optional
from EngineComponents.Abstract.ElectricalComponent import ElectricalComponent


@dataclass
class Battery(ElectricalComponent):
    specific_energy: float  # J/kg
    battery_packing_factor: float  # -
    burn_time: float  # s

    # Overwrite ElectricComponent attribute
    electric_energy_efficiency: Optional[float] = field(init=False, default=None)

    def __post_init__(self):
        self.electric_energy_efficiency = self.eta_e

    @property
    def eta_e(self):
        return min(0.093 * log(self.burn_time) + 0.3301, .99)

    @property
    def total_energy(self):
        return self.input_power * self.burn_time

    @property
    def energy_heat_loss(self):
        return self.total_energy * (1 - self.eta_e)

    @property
    def power_heat_loss(self):
        return self.input_power * (1 - self.eta_e)

    @property
    def optimal_time(self):
        return self.specific_energy / self.specific_power

    @property
    def power_mass(self):
        return self.battery_packing_factor * self.output_power / self.specific_power

    @property
    def energy_mass(self):
        return self.battery_packing_factor * self.output_power * self.burn_time / (self.specific_energy * self.eta_e)

    @property  # Overwrite ElectricComponent mass
    def mass(self):
        return max(self.power_mass, self.energy_mass)