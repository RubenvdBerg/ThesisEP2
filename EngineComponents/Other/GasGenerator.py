from EngineComponents.Abstract.PressureComponent import PressureComponent
from EngineComponents.Abstract.FlowState import ManualFlowState, DefaultFlowState, FlowState
from dataclasses import dataclass
from typing import Optional


@dataclass
class GasGenerator(PressureComponent):
    """Handles mass estimation of the Gas Generator and merging of inlet flows.

    mass_mixture_ratio is required for clarity, so it is available as attribute for other components"""
    stay_time: float = 0  # [s]
    pressure: float = 0  # [Pa]
    combustion_temperature: float = 0  # [K]
    oxidizer_inlet_flow_state: FlowState = DefaultFlowState()
    fuel_inlet_flow_state: FlowState = DefaultFlowState()
    mass_mixture_ratio: float = 0  # [-]

    # Gas properties
    specific_heat_capacity: Optional[float] = None  # [J/(kgK)]
    heat_capacity_ratio: Optional[float] = None  # [-]
    molar_mass: Optional[float] = None  # [kg/mol]
    gas_density: Optional[float] = None  # [kg/m3]

    @property
    def outlet_mass_flow(self) -> float:
        return sum((self.oxidizer_inlet_flow_state.mass_flow, self.fuel_inlet_flow_state.mass_flow))

    @property
    def outlet_flow_state(self) -> FlowState:
        return ManualFlowState(propellant_name='ExhaustGas',
                               temperature=self.combustion_temperature,
                               pressure=self.pressure,
                               mass_flow=self.outlet_mass_flow,
                               type='combusted',
                               _specific_heat_capacity=self.specific_heat_capacity,
                               _heat_capacity_ratio=self.heat_capacity_ratio,
                               _molar_mass=self.molar_mass,
                               _density=self.gas_density,)

    # Mass Calculation Properties
    @property
    def maximum_expected_operating_pressure(self):
        return self.pressure

    @property
    def volume(self):
        return self.stay_time * self.outlet_mass_flow / self.gas_density
