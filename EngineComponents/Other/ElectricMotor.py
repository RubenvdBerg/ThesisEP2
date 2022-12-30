import warnings
from dataclasses import dataclass

from EngineComponents.Abstract.ElectricComponent import ElectricComponent
from EngineComponents.Abstract.FlowState import FlowState


@dataclass
class ElectricMotor(ElectricComponent):
    electric_heat_loss_factor: float
    oxidizer_pump_inlet_flow_state: FlowState
    oxidizer_leakage_factor: float
    magnet_temp_limit: float

    @property
    def power_heat_loss(self):
        return self.output_power * .015

    def calc_cooling(self):
        cp_ox = self.oxidizer_pump_inlet_flow_state.specific_heat_capacity
        m_ox_leak = self.oxidizer_pump_inlet_flow_state.mass_flow * self.oxidizer_leakage_factor
        deltaT = (self.magnet_temp_limit - 50) - self.oxidizer_pump_inlet_flow_state.temperature
        required_mass_flow = self.power_heat_loss / (cp_ox * deltaT)
        if required_mass_flow > m_ox_leak:
            warnings.warn('The expected cooling due to oxidizer leakage through the electric motor is too low. The motor magnets will be too hot and possibly demagnetize.')
            print(f'-------------------------------------------------{required_mass_flow} required, {m_ox_leak} expected')