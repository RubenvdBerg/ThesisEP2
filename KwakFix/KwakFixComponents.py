import warnings
from dataclasses import dataclass
from math import pi
from typing import Optional

from BaseEngineCycle.Propellant import Propellant
from BaseEngineCycle.Structure import PressureStructure
from BaseEngineCycle.Tank import Tank
from ElectricPumpCycle.EPComponents import Battery


@dataclass
class KwakTank(Tank):
    _kwak_fix_cycle_type: str = ''

    @property  # Override Tank property
    def initial_head(self):
        if self._kwak_fix_cycle_type == 'ep':
            if self.inlet_flow_state.type == 'oxidizer':
                return 1.91839449096392
            elif self.inlet_flow_state.type == 'fuel':
                return 1.65560478870526
        elif self._kwak_fix_cycle_type == 'gg':
            if self.inlet_flow_state.type == 'oxidizer':
                return 1.92539045861846
            elif self.inlet_flow_state.type == 'fuel':
                return 1.71033897378923
        else:
            return super().initial_head


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
        return (self.battery_packing_factor * self.output_power
                / (self.electric_motor_efficiency * self.inverter_efficiency)
                * max(1 / self.specific_power, self.burn_time / (self.specific_energy * self.eta_e)))
