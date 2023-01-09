from dataclasses import dataclass


@dataclass
class ElectricalComponent:
    electric_energy_efficiency: float  # -
    specific_power: float  # W/kg
    output_power: float  # W

    @property
    def input_power(self):
        return self.output_power / self.electric_energy_efficiency

    @property
    def mass(self):
        return self.output_power / self.specific_power