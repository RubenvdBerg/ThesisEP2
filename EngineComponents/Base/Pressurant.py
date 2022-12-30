from dataclasses import dataclass
from scipy import constants as constants
from EngineComponents.Abstract.PressureComponent import PressureComponent


@dataclass
class Pressurant:
    fuel_volume: float  # [m3]
    oxidizer_volume: float  # [m3]
    propellant_tanks_ullage_factor: float  # [Pa]
    fuel_tank_initial_pressure: float  # [Pa]
    oxidizer_tank_initial_pressure: float  # [Pa]
    margin_factor: float  # [-]
    initial_pressure: float  # [Pa]
    final_pressure: float  # [Pa]
    heat_capacity_ratio: float  # [-]
    molar_mass: float  # [kg/mol]!!
    initial_temperature: float  # [K]

    @property
    def specific_gas_constant(self):
        return round(constants.R / self.molar_mass)

    @property
    def mass(self):
        fact_m = self.margin_factor
        fact_u = self.propellant_tanks_ullage_factor
        y = self.heat_capacity_ratio
        R_sp = self.specific_gas_constant
        T_0 = self.initial_temperature
        otp = self.oxidizer_tank_initial_pressure
        ftp = self.fuel_tank_initial_pressure
        ov = self.oxidizer_volume
        fv = self.fuel_volume
        p0 = self.initial_pressure
        p1 = self.final_pressure
        return fact_m * fact_u * y / (R_sp * T_0) * (otp * ov + ftp * fv) / (1 - (p1 / p0))


@dataclass
class PressurantTank(PressureComponent):
    pressurant: Pressurant

    @property
    def initial_pressure(self):
        return self.pressurant.initial_pressure

    @property
    def final_pressure(self):
        return self.pressurant.final_pressure

    @property
    def volume(self):
        return (self.pressurant.mass * self.pressurant.specific_gas_constant
                * self.pressurant.initial_temperature / self.pressurant.initial_pressure)

    @property
    def maximum_expected_operating_pressure(self):
        return self.initial_pressure