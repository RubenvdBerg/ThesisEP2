from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from dataclasses import dataclass
from scipy import constants
from KwakFix.KwakFixComponents import KwakBattery


@dataclass
class KwakEngineCycle(EngineCycle):

    def set_pump_outlet_pressures(self):
        Merger._warn_pressure = False
        Merger._warn_pressure = True

    @property
    def oxidizer_pump_outlet_pressure(self):
        return self.combustion_chamber_pressure * self._oxidizer_pump_pressure_factor_first_guess

    @property
    def fuel_pump_outlet_pressure(self):
        return self.combustion_chamber_pressure * self._fuel_pump_pressure_factor_first_guess


@dataclass
class KwakFixGasGeneratorCycle(GasGeneratorCycle, EngineCycle):
    _iteration_done: bool = False

    def set_initial_values(self):
        super().set_initial_values()
        if self.gg_gas_molar_mass:
            r = constants.gas_constant / self.gg_gas_molar_mass
            self.gg_gas_density = self.gg_pressure / (r * self.turbine_maximum_temperature)

    def iterate_mass_flow(self):
        m_tu = self.turbine.mass_flow_required
        m_gg_o = m_tu * self.gg_mass_mixture_ratio / (1 + self.gg_mass_mixture_ratio)
        m_gg_f = m_tu * 1 / (1 + self.gg_mass_mixture_ratio)
        self._iterative_turbine_mass_flow += m_gg_f + m_gg_o
        m_tu = self.turbine.mass_flow_required
        m_gg_o = m_tu * self.gg_mass_mixture_ratio / (1 + self.gg_mass_mixture_ratio)
        m_gg_f = m_tu * 1 / (1 + self.gg_mass_mixture_ratio)
        # self._iterative_turbine_mass_flow = m_gg_f + m_gg_o
        # m_tu = self.turbine.mass_flow_required
        # m_gg_o = m_tu * self.gg_mass_mixture_ratio / (1 + self.gg_mass_mixture_ratio)
        # m_gg_f = m_tu * 1 / (1 + self.gg_mass_mixture_ratio)
        self.mo = self.main_oxidizer_flow - m_gg_o
        self.mf = self.main_fuel_flow - m_gg_f
        self._iterative_turbine_mass_flow = m_tu
        self._iteration_done = True


    @property
    def turbine_mass_flow_initial_guess(self):
        return 0

    @property
    def overall_specific_impulse(self):
        """Calculate specific impulse without accounting for turbine exhaust thrust contribution."""
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

    @property
    def battery(self):
        return KwakBattery(specific_power=self.battery_specific_power,
                           specific_energy=self.battery_specific_energy,
                           battery_packing_factor=self.battery_structural_factor,
                           output_power=self.pumps_power_required,
                           burn_time=self.burn_time,
                           inverter_efficiency=self.inverter.electric_energy_efficiency,
                           electric_motor_efficiency=self.electric_motor.electric_energy_efficiency)
