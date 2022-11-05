from dataclasses import dataclass
from EngineCycles.ElectricPumpCycle.EPComponents import Battery


@dataclass
class KwakBattery(Battery):
    electric_motor_efficiency: float  # [-]
    inverter_efficiency: float  # [-]

    @property
    def total_energy(self):
        return self.output_power * self.burn_time

    # Overwrite Super _mass_
    @property
    def mass(self):
        return (self.battery_packing_factor * self.output_power / (self.electric_motor_efficiency * self.inverter_efficiency)
                * max(1 / self.specific_power, self.burn_time / (self.specific_energy * self.eta_e)))
