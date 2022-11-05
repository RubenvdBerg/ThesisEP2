from EngineCycles.BaseEngineCycle.EngineCycle import EngineCycle
from EngineCycles.GasGeneratorCycle.GGCycle import GasGeneratorCycle
from EngineCycles.ElectricPumpCycle.EPCycle import ElectricPumpCycle
from EngineCycles.BaseOpenCycle import OpenCycle
from dataclasses import dataclass
from scipy import constants
from KwakFix.KwakFixComponents import KwakBattery
from typing import Literal


@dataclass
class KwakFixGasGeneratorCycle(GasGeneratorCycle, EngineCycle):
    _iteration_done: bool = False

    def __post_init__(self):
        super().__post_init__()

    def iterate_mass_flow(self):
        m_tu = self.turbine.mass_flow_required
        m_gg_o = m_tu * self.gg_mass_mixture_ratio / (1 + self.gg_mass_mixture_ratio)
        m_gg_f = m_tu * 1 / (1 + self.gg_mass_mixture_ratio)
        self._iterative_turbine_mass_flow += m_gg_f + m_gg_o
        m_tu = self.turbine.mass_flow_required
        m_gg_o = m_tu * self.gg_mass_mixture_ratio / (1 + self.gg_mass_mixture_ratio)
        m_gg_f = m_tu * 1 / (1 + self.gg_mass_mixture_ratio)
        self.mo = self.main_oxidizer_flow - m_gg_o
        self.mf = self.main_fuel_flow - m_gg_f
        self._iterative_turbine_mass_flow = m_tu
        self._iteration_done = True

    @property
    def turbine_mass_flow_initial_guess(self):
        return 0

    @property
    def total_mass_flow(self):
        return self.main_fuel_flow + self.main_oxidizer_flow

    @property
    def overall_specific_impulse(self):
        return self.chamber_mass_flow * self.ideal_thrust_coefficient * self.characteristic_velocity / self.total_mass_flow / constants.g

    @property
    def chamber_mass_flow(self):
        if self._iteration_done:
            return self.mf + self.mo
        else:
            return self.base_mass_flow

    @property
    def chamber_fuel_flow(self):
        if self._iteration_done:
            return self.mf
        else:
            return super().chamber_fuel_flow

    @property
    def chamber_oxidizer_flow(self):
        if self._iteration_done:
            return self.mo
        else:
            return super().chamber_oxidizer_flow

    @property
    def gg_propellant_mass(self):
        return self.gg_mass_flow * self.burn_time


@dataclass
class KwakFixElectricPumpCycle(ElectricPumpCycle, EngineCycle):
    pass
    # @property
    # def battery(self):
    #     return KwakBattery(specific_power=self.battery_specific_power,
    #                        specific_energy=self.battery_specific_energy,
    #                        battery_packing_factor=self.battery_structural_factor,
    #                        output_power=self.pump_power_required,
    #                        burn_time=self.burn_time,
    #                        fuel_specific_heat=self.fuel_specific_heat,
    #                        coolant_allowable_temperature_change=self.battery_coolant_temperature_change,
    #                        inverter_efficiency=self.inverter.electric_energy_efficiency,
    #                        electric_motor_efficiency=self.electric_motor.electric_energy_efficiency)
