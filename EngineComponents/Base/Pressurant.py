from dataclasses import dataclass, field
from EngineComponents.Abstract.PressureComponent import NewPressureComponent
from EngineComponents.Abstract.FlowState import FlowState


@dataclass
class Pressurant:
    fuel_volume: float  # [m3]
    oxidizer_volume: float  # [m3]
    propellant_tanks_ullage_factor: float  # [Pa]
    fuel_tank_initial_pressure: float  # [Pa]
    oxidizer_tank_initial_pressure: float  # [Pa]
    margin_factor: float  # [-]
    final_pressure: float  # [Pa]
    initial_fluid_state: FlowState

    @property
    def mass(self):
        fact_m = self.margin_factor
        fact_u = self.propellant_tanks_ullage_factor
        y = self.initial_fluid_state.heat_capacity_ratio
        R_sp = self.initial_fluid_state.specific_gas_constant
        T_0 = self.initial_fluid_state.temperature
        otp = self.oxidizer_tank_initial_pressure
        ftp = self.fuel_tank_initial_pressure
        ov = self.oxidizer_volume
        fv = self.fuel_volume
        p0 = self.initial_fluid_state.pressure
        p1 = self.final_pressure
        return fact_m * fact_u * y / (R_sp * T_0) * (otp * ov + ftp * fv) / (1 - (p1 / p0))

    @property
    def volume(self):
        return (self.mass * self.initial_fluid_state.specific_gas_constant
                * self.initial_fluid_state.temperature / self.initial_fluid_state.pressure)


@dataclass
class PressurantTank(NewPressureComponent):
    pressurant: Pressurant
    geometry: str = field(init=False, default='sphere')

    @property
    def volume(self):
        return self.pressurant.volume

    @property
    def max_pressure(self):
        return self.pressurant.initial_fluid_state.pressure
