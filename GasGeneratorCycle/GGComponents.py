from BaseEngineCycle.Structure import Structure
from BaseEngineCycle.FlowComponent import BaseFlowComponent
from BaseEngineCycle.FlowState import FlowState, DefaultFlowState
from dataclasses import dataclass


@dataclass
class GasGenerator(Structure):
    """Handles mass estimation of the Gas Generator and merging of inlet flows. mass_mixture_ratio is required for
    clarity, so it is available as attribute for other components, as well as inheritance

    !!Does not estimate temperature and pressure changes within GasGenerator, but simply assumes the outlet flow to be
    equal to the given inlet conditions of the turbine!!
    """
    gas_constant: float = 0  # [J/(kg*K)]
    stay_time: float = 0  # [s]
    turbine_temp_limit: float = 0  # [K]
    turbine_inlet_pressure: float = 0  # [Pa]
    mass_mixture_ratio: float = 0  # [-]
    oxidizer_inlet_flow_state: FlowState = DefaultFlowState()
    fuel_inlet_flow_state: FlowState = DefaultFlowState()

    @property
    def outlet_temperature(self):
        """Gas generator operated at limit of turbine, temperature wise"""
        return self.turbine_temp_limit

    @property
    def outlet_pressure(self):
        return self.turbine_inlet_pressure

    @property
    def outlet_mass_flow(self) -> float:
        return sum((self.oxidizer_inlet_flow_state.mass_flow, self.fuel_inlet_flow_state.mass_flow))

    @property
    def outlet_flow_state(self) -> FlowState:
        return FlowState(propellant_name='ExhaustGas',
                         temperature=self.outlet_temperature,
                         pressure=self.outlet_pressure,
                         mass_flow=self.outlet_mass_flow,
                         type='burnt')

    @property
    def mass(self):  # Overwrite Structure.mass()
        return (self.safety_factor * 3 / 2 * self.material_density / self.yield_strength
                * self.stay_time * self.inlet_pressure / self.gas_density * self.mass_flow)

    @property
    def gas_density(self):
        """Estimate gas generator exhaust gas density with ideal gas law"""
        return self.inlet_pressure / (self.gas_constant * self.turbine_temp_limit)


@dataclass
class ComplexGasGenerator(GasGenerator):
    pass