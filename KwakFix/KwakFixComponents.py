from dataclasses import dataclass
from EngineComponents.Other.Battery import Battery
from EngineComponents.Base.Tank import Tank
from EngineComponents.Base.Pump import Pump
from EngineComponents.Base.Propellant import Propellant


@dataclass
class KwakBattery(Battery):
    electric_motor_efficiency: float  # [-]
    inverter_efficiency: float  # [-]

    @property
    def total_energy(self):
        return self.output_power * self.burn_time

    @property
    def power_heat_loss(self):
        return self.output_power * (1 - self.eta_e)

    @property
    def mass(self):
        return (self.battery_packing_factor * self.output_power / (
                    self.electric_motor_efficiency * self.inverter_efficiency)
                * max(1 / self.specific_power, self.burn_time / (self.specific_energy * self.eta_e)))


@dataclass
class KwakTank(Tank):
    manual_propellant_density: float = 0

    @property
    def propellant_density(self):
        return self.manual_propellant_density


@dataclass
class KwakPump(Pump):
    manual_propellant_density: float = 0

    @property
    def propellant_density(self):
        return self.manual_propellant_density


@dataclass
class KwakPropellant(Propellant):
    manual_propellant_density: float = 0

    @property
    def density(self):
        return self.manual_propellant_density
