from BaseEngineCycle.EngineCycle import EngineCycle
from GasGeneratorCycle.GGCycle import GasGeneratorCycle
from ElectricPumpCycle.EPCycle import ElectricPumpCycle
from dataclasses import dataclass
from KwakFix.KwakFixComponents import KwakTank, KwakBattery
from typing import Literal


@dataclass
class KwakFixEngineCycle(EngineCycle):
    _kwak_fix_cycle_type: Literal['ep', 'gg', 'base'] = 'base'

    @property
    def oxidizer_tank(self):
        return KwakTank(material_density=self.tanks_material_density,
                        safety_factor=self.tanks_structural_factor,
                        yield_strength=self.tanks_yield_strength,
                        propellant=self.oxidizer,
                        max_acceleration=self.max_acceleration,
                        ullage_factor=self.ullage_volume_factor,
                        initial_pressure=self.oxidizer_initial_pressure,
                        pressurant_tank_volume=self.pressurant_tank.volume,
                        _kwak_fix_cycle_type=self._kwak_fix_cycle_type)

    @property
    def fuel_tank(self):
        return KwakTank(max_acceleration=self.max_acceleration,
                        ullage_factor=self.ullage_volume_factor,
                        propellant=self.fuel,
                        pressurant_tank_volume=None,
                        initial_pressure=self.fuel_initial_pressure,
                        material_density=self.tanks_material_density,
                        yield_strength=self.tanks_yield_strength,
                        safety_factor=self.tanks_structural_factor,
                        _kwak_fix_cycle_type=self._kwak_fix_cycle_type)


@dataclass
class KwakFixGasGeneratorCycle(GasGeneratorCycle, KwakFixEngineCycle):
    _kwak_fix_cycle_type: str = 'gg'

    def __post_init__(self):
        super().__post_init__()
        self.exhaust_thrust_contribution = 0.002137040978335770

    @property
    def gg_propellant_mass(self):
        return self.gg_mass_flow * self.burn_time


@dataclass
class KwakFixElectricPumpCycle(ElectricPumpCycle, KwakFixEngineCycle):
    _kwak_fix_cycle_type: str = 'ep'

    @property
    def battery(self):
        return KwakBattery(specific_power=self.battery_specific_power,
                           specific_energy=self.battery_specific_energy,
                           battery_packing_factor=self.battery_structural_factor,
                           output_power=self.electric_motor.output_power,
                           burn_time=self.burn_time,
                           fuel_specific_heat=self.fuel_specific_heat,
                           coolant_allowable_temperature_change=self.coolant_allowable_temperature_change,
                           inverter_efficiency=self.inverter.electric_energy_efficiency,
                           electric_motor_efficiency=self.electric_motor.electric_energy_efficiency)
