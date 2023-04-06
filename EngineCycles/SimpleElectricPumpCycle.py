import warnings
from dataclasses import dataclass, replace, field
from typing import Optional
from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineComponents.Base.Pump import Pump
from EngineComponents.Base.Tank import Tank
from EngineComponents.Other.ElectricMotor import SimpleElectricMotor
from EngineComponents.Other.Battery import Battery
from EngineComponents.Other.Inverter import Inverter


@dataclass
class SimpleElectricPumpCycle(EngineCycle):
    electric_motor_specific_power: float = 0  # W/kg
    inverter_specific_power: float = 0  # W/kg
    battery_specific_power: float = 0  # W/kg
    battery_specific_energy: float = 0  # J/kg
    electric_motor_efficiency: float = 0  # -
    inverter_efficiency: float = 0  # -
    battery_structural_factor: float = 0  # -
    iterate: bool = False

    @property
    def fuel_pump(self):
        return Pump(inlet_flow_state=self.fuel_tank.outlet_flow_state,
                    expected_outlet_pressure=self.fuel_pump_outlet_pressure,
                    efficiency=self.fuel_pump_efficiency,
                    specific_power=self.fuel_pump_specific_power, )

    @property
    def electric_motor(self):
        return SimpleElectricMotor(specific_power=self.electric_motor_specific_power,
                                   electric_energy_efficiency=self.electric_motor_efficiency,
                                   output_power=self.pumps_power_required, )

    @property
    def inverter(self):
        return Inverter(specific_power=self.inverter_specific_power,
                        electric_energy_efficiency=self.inverter_efficiency,
                        output_power=self.electric_motor.input_power)

    @property
    def battery(self):
        return Battery(specific_power=self.battery_specific_power,
                       specific_energy=self.battery_specific_energy,
                       battery_packing_factor=self.battery_structural_factor,
                       output_power=self.inverter.input_power,
                       burn_time=self.burn_time, )

    @property
    def cooling_inlet_flow_state(self):
        return self.fuel_pump.outlet_flow_state

    @property
    def oxidizer_tank(self):
        return Tank(inlet_flow_state=self.oxidizer_main_flow_state,
                    propellant_volume=self.oxidizer.volume,
                    max_acceleration=self.max_acceleration,
                    ullage_factor=self.ullage_volume_factor,
                    pressurant_tank_volume=None,
                    structure_material=self.oxidizer_tank_material,
                    safety_factor=self.tanks_structural_factor, )

    @property
    def tanks_mass(self):
        return self.oxidizer_tank.mass + self.fuel_tank.mass

    @property
    def feed_system_mass(self):
        return super().feed_system_mass + self.electric_motor.mass + self.inverter.mass

    @property
    def dry_mass(self):
        return super().dry_mass + self.battery.mass

    @property
    def mass_kwak(self):
        return super().mass_kwak + self.battery.mass + self.inverter.mass + self.electric_motor.mass - self.pressurant.mass
