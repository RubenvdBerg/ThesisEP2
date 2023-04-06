import warnings
from dataclasses import dataclass

from EngineComponents.Abstract.ElectricalComponent import ElectricalComponent
from EngineComponents.Abstract.FlowState import FlowState

@dataclass
class SimpleElectricMotor(ElectricalComponent):
    pass

@dataclass
class ElectricMotor(ElectricalComponent):
    oxidizer_pump_inlet_flow_state: FlowState
    oxidizer_leakage_factor: float
    magnet_temp_limit: float
    electric_heat_loss_factor: float

    @property
    def power_heat_loss(self):
        return self.output_power * self.electric_heat_loss_factor

    def calc_cooling(self):
        cp_ox = self.oxidizer_pump_inlet_flow_state.specific_heat_capacity
        m_ox_leak = self.oxidizer_pump_inlet_flow_state.mass_flow * self.oxidizer_leakage_factor
        temp_margin = 50
        deltaT = (self.magnet_temp_limit - temp_margin) - self.oxidizer_pump_inlet_flow_state.temperature
        required_mass_flow = self.power_heat_loss / (cp_ox * deltaT)
        if required_mass_flow > m_ox_leak:
            warnings.warn('The expected cooling due to assumed oxidizer leakage through the electric motor is too low. '
                          'The motor magnets will be too hot and will possibly demagnetize.')
