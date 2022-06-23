from BaseEngineCycle.EngineCycle import EngineCycle
from BaseEngineCycle.Turbine import Turbine
from dataclasses import dataclass, field
from scipy.constants import g


# Abstract class that has the attributes GasGenerator and OpenExpander share
@dataclass
class OpenEngineCycle(EngineCycle):
    turbine_gas_specific_heat_capacity: float = 0  # [J/(kg*K)]
    turbine_gas_heat_capacity_ratio: float = 0  # [-]
    turbine_pressure_ratio: float = 0  # [-]
    turbine_efficiency: float = 0  # [-]
    turbopump_specific_power: float = 0  # [W/kg]
    exhaust_thrust_contribution: float = 0  # [-]

    # Iteration attribute, not required at init
    turbine_mass_flow: float = field(init=False, default=0)  # [kg/s]

    # Override from EngineCycle, are assigned later, not required at init
    fuel_pump_specific_power: float = field(init=False, repr=False, default=None)
    oxidizer_pump_specific_power: float = field(init=False, repr=False, default=None)

    def __post_init__(self):
        super().__post_init__()
        # Total turbine, oxidizer-pump, fuel-pump-combinations seen as a single turbopump with single specific power
        self.fuel_pump_specific_power = self.oxidizer_pump_specific_power = self.turbopump_specific_power
        if self.iterate:
            self.iterate_mass_flow()

    def iterate_mass_flow(self):
        if self.verbose:
            print('Iterate Turbine Mass Flow')
        while self.turbine_mass_flow * 1.001 < self.turbine.mass_flow_required:
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

    @property
    def turbine_inlet_temperature(self):
        return None

    @property
    def turbine(self):
        return Turbine(pump_power_required=self.pump_power_required,
                       inlet_temperature=self.turbine_inlet_temperature,
                       efficiency=self.turbine_efficiency,
                       specific_heat_capacity=self.turbine_gas_specific_heat_capacity,
                       heat_capacity_ratio=self.turbine_gas_heat_capacity_ratio,
                       pressure_ratio=self.turbine_pressure_ratio)

    @property
    def chamber_mass_flow(self):
        return (1 - self.exhaust_thrust_contribution) * self.mass_flow

    @property
    def chamber_fuel_flow(self):
        return 1 / (self.mass_mixture_ratio + 1) * self.chamber_mass_flow

    @property
    def chamber_oxidizer_flow(self):
        return self.mass_mixture_ratio / (self.mass_mixture_ratio + 1) * self.mass_flow

    @property
    def total_mass_flow(self):
        return self.chamber_mass_flow + self.turbine_mass_flow

    @property  # Override EngineCycle.simple_specific_impulse to use total mass flow
    def simple_specific_impulse(self):
        return self.thrust / self.total_mass_flow / g
