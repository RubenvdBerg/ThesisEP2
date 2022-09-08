from BaseEngineCycle.EngineCycle import EngineCycle
from BaseOpenCycle.Turbine import Turbine
from BaseOpenCycle.SecondaryExhaust import SecondaryExhaust
from BaseEngineCycle.FlowState import DefaultFlowState, FlowState
from dataclasses import dataclass, field, replace
from scipy.constants import g
from typing import Optional
import warnings
from functools import cached_property


# Abstract class that has the attributes GasGenerator and OpenExpander share
@dataclass
class OpenEngineCycle(EngineCycle):
    turbine_efficiency: float = 0  # [-]
    turbopump_specific_power: float = 0  # [W/kg]
    exhaust_thrust_contribution: float = 0  # [-]
    exhaust_expansion_ratio: float = 0  # [-]
    turbine_pressure_ratio: Optional[float] = None  # [-]
    turbine_inlet_pressure: Optional[float] = None  # [Pa]
    turbine_outlet_pressure: Optional[float] = None  # [Pa]

    turbine_gas_specific_heat_capacity: Optional[float] = None  # [J/kg]
    turbine_gas_heat_capacity_ratio: Optional[float] = None  # [-]
    turbine_gas_molar_mass: Optional[float] = None

    # Iteration attribute, not required at init
    _iterative_turbine_mass_flow: float = field(init=False, repr=False)  # [kg/s]

    # Override from EngineCycle, are assigned later, not required at init
    fuel_pump_specific_power: float = field(init=False, repr=False, default=None)
    oxidizer_pump_specific_power: float = field(init=False, repr=False, default=None)

    def __post_init__(self):
        super().__post_init__()
        # Combination of turbine, oxidizer-pump and fuel-pump seen as a single turbopump with single specific power
        self.fuel_pump_specific_power = self.oxidizer_pump_specific_power = self.turbopump_specific_power
        self.resolve_turbine_pressure_ratio_choice()
        self._iterative_turbine_mass_flow = self.turbine_mass_flow_initial_guess
        if self.iterate:
            self.iterate_mass_flow()
        self.check_exhaust_thrust_contribution()

    def turbine_flow_error_larger_than_accuracy(self):
        error = abs(self.turbine.mass_flow_required - self._iterative_turbine_mass_flow)
        margin = self.turbine.mass_flow_required * self.iteration_accuracy
        return error > margin

    def iterate_mass_flow(self):
        if self.verbose:
            print('Iterate Turbine Mass Flow')
        while self.turbine_flow_error_larger_than_accuracy():
            if self.verbose:
                print(f'Actual:  {self._iterative_turbine_mass_flow:.5f} kg/s')
                print(f'Required:{self.turbine.mass_flow_required:.5f} kg/s\n')
            self._iterative_turbine_mass_flow = self.turbine.mass_flow_required
        if self.verbose:
            print(f'Turbine Mass Flow Set\n')

    def reiterate(self):
        if self.verbose:
            print('Start reiteration')
        self.update_cea()
        self.iterate_mass_flow()

    def resolve_turbine_pressure_ratio_choice(self):
        if not ((self.turbine_pressure_ratio is None) ^ (self.turbine_outlet_pressure is None)):
            raise ValueError('Both or neither turbine_pressure_ratio and turbine_outlet_pressure are provided. '
                             'Provide one and only one')
        elif self.turbine_pressure_ratio is None:
            self.turbine_pressure_ratio = self._turbine_inlet_pressure / self.turbine_outlet_pressure

    def check_exhaust_thrust_contribution(self):
        print(f'Guessed exhaust thrust contribution {self.exhaust_thrust_contribution * self.thrust * 1e-3:.3f} kN')
        print(f'Estimated exhaust thrust contribution {self.turbine_exhaust.thrust * 1e-3:.3f} kN')

    @cached_property
    def _turbine_inlet_pressure(self):
        if self.turbine_inlet_pressure is None:
            warnings.warn('No turbine inlet pressure provided, estimated to be equal to the pressure of the flow '
                          'before the main chamber injector')
            return self.combustion_chamber_pressure * (1 + self._injector_pressure_drop_factor)
        else:
            return self.turbine_inlet_pressure

    @property
    def turbine_mass_flow_initial_guess(self):
        return 0.0

    @property
    def turbine_base_inlet_flow_state(self):
        return FlowState(propellant_name=self.fuel_name,
                         temperature=0.0,
                         pressure=self._turbine_inlet_pressure,
                         mass_flow=None,
                         type='fuel', )

    @property
    def turbine_inlet_flow_state(self):
        return replace(self.turbine_base_inlet_flow_state,
                       mass_flow=self._iterative_turbine_mass_flow)

    @property
    def turbine(self):
        return Turbine(inlet_flow_state=self.turbine_inlet_flow_state,
                       pump_power_required=self.pump_power_required,
                       efficiency=self.turbine_efficiency,
                       pressure_ratio=self.turbine_pressure_ratio,
                       gas_heat_capacity_ratio=self.turbine_gas_heat_capacity_ratio,
                       gas_specific_heat_capacity=self.turbine_gas_specific_heat_capacity,)

    @property
    def turbine_exhaust(self):
        return SecondaryExhaust(inlet_flow_state=self.turbine.outlet_flow_state,
                                expansion_ratio=self.exhaust_expansion_ratio,
                                ambient_pressure=self.ambient_pressure,
                                gas_heat_capacity_ratio=self.turbine_gas_heat_capacity_ratio,
                                gas_molar_mass=self.turbine_gas_molar_mass,)

    @property
    def chamber_mass_flow(self):
        return (1 - self.exhaust_thrust_contribution) * self.base_mass_flow

    @property
    def total_mass_flow(self):
        return self.chamber_mass_flow + self._iterative_turbine_inlet_flow_state.mass_flow

    @property  # Override EngineCycle.simple_specific_impulse to use total mass flow
    def simple_specific_impulse(self):
        return self.thrust / self.total_mass_flow / g
