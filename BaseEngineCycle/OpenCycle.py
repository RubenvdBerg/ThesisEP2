from BaseEngineCycle.EngineCycle import EngineCycle
from BaseEngineCycle.Turbine import Turbine
from BaseEngineCycle.SecondaryExhaust import SecondaryExhaust
from dataclasses import dataclass, field
from scipy.constants import g
from typing import Optional
import warnings
from functools import cached_property


# Abstract class that has the attributes GasGenerator and OpenExpander share
@dataclass
class OpenEngineCycle(EngineCycle):
    turbine_gas_specific_heat_capacity: float = 0  # [J/(kg*K)]
    turbine_gas_heat_capacity_ratio: float = 0  # [-]

    turbine_efficiency: float = 0  # [-]
    turbopump_specific_power: float = 0  # [W/kg]
    exhaust_thrust_contribution: float = 0  # [-]
    exhaust_expansion_ratio: float = 0  # [-]
    turbine_pressure_ratio: Optional[float] = None  # [-]
    turbine_inlet_pressure: Optional[float] = None  # [Pa]
    turbine_outlet_pressure: Optional[float] = None  # [Pa]

    # Iteration attribute, not required at init
    turbine_mass_flow: float = field(init=False, default=0)  # [kg/s]

    # Override from EngineCycle, are assigned later, not required at init
    fuel_pump_specific_power: float = field(init=False, repr=False, default=None)
    oxidizer_pump_specific_power: float = field(init=False, repr=False, default=None)

    def __post_init__(self):
        super().__post_init__()
        # Total turbine, oxidizer-pump, fuel-pump-combinations seen as a single turbopump with single specific power
        self.fuel_pump_specific_power = self.oxidizer_pump_specific_power = self.turbopump_specific_power
        if self.turbine_inlet_pressure is None:
            self.turbine_inlet_pressure = self.combustion_chamber_pressure + self.injector.pressure_drop
            warnings.warn('No turbine inlet pressure provided, estimated to be equal to the pressure of the flow '
                          'before the main chamber injector')
        self.resolve_turbine_pressure_ratio_choice()
        self.turbine_mass_flow = self.turbine_mass_flow_initial_guess
        if self.iterate:
            self.iterate_mass_flow()
        self.check_exhaust_thrust_contribution()

    def iterate_mass_flow(self):
        if self.verbose:
            print('Iterate Turbine Mass Flow')
        while abs(self.turbine.mass_flow_required
                  - self.turbine_mass_flow) > self.turbine.mass_flow_required * self.iteration_accuracy:
            if self.verbose:
                print(f'Actual:  {self.turbine_mass_flow:.5f} kg/s')
                print(f'Required:{self.turbine.mass_flow_required:.5f} kg/s\n')
            self.turbine_mass_flow = self.turbine.mass_flow_required
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
            self.turbine_pressure_ratio = self.turbine_inlet_pressure / self.turbine_outlet_pressure

    def check_exhaust_thrust_contribution(self):
        print(f'Expected exhaust thrust contribution {self.exhaust_thrust_contribution*self.thrust*1e-3:.2f} kN')
        print(f'Estimated exhaust thrust contribution {self.turbine_exhaust.thrust*1e-3:.2f} kN')

    @property
    def turbine_mass_flow_initial_guess(self):
        return 0.0

    @property
    def turbine_inlet_temperature(self):
        return None

    @cached_property
    def turbine_gas_coolprop_name(self):
        return self.coolprop_name(self.fuel.name)

    @property
    def turbine(self):
        return Turbine(pump_power_required=self.pump_power_required,
                       inlet_temperature=self.turbine_inlet_temperature,
                       efficiency=self.turbine_efficiency,
                       specific_heat_capacity=self.turbine_gas_specific_heat_capacity,
                       heat_capacity_ratio=self.turbine_gas_heat_capacity_ratio,
                       pressure_ratio=self.turbine_pressure_ratio,
                       inlet_pressure=self.turbine_inlet_pressure,
                       coolprop_name=self.turbine_gas_coolprop_name)
    @property
    def turbine_exhaust(self):
        return SecondaryExhaust(inlet_pressure=self.turbine.outlet_pressure,
                                inlet_temperature=self.turbine.outlet_temperature,
                                expansion_ratio=self.exhaust_expansion_ratio,
                                mass_flow=self.turbine_mass_flow,
                                coolprop_name=self.turbine_gas_coolprop_name)

    @property
    def chamber_mass_flow(self):
        return (1 - self.exhaust_thrust_contribution) * self.base_mass_flow

    @property
    def total_mass_flow(self):
        return self.chamber_mass_flow + self.turbine_mass_flow

    @property  # Override EngineCycle.simple_specific_impulse to use total mass flow
    def simple_specific_impulse(self):
        return self.thrust / self.total_mass_flow / g

