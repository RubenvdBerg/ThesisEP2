from EngineCycles.BaseEngineCycle.EngineCycle import EngineCycle
from EngineCycles.BaseOpenCycle.Turbine import Turbine
from EngineCycles.BaseOpenCycle.SecondaryExhaust import SecondaryExhaust
from EngineCycles.BaseEngineCycle.FlowState import FlowState
from dataclasses import dataclass, field
from scipy.constants import g
from typing import Optional
import warnings


# Abstract class that has the attributes GasGenerator and OpenExpander share
@dataclass
class OpenEngineCycle(EngineCycle):
    turbine_efficiency: float = 0  # [-]
    turbopump_specific_power: float = 0  # [W/kg]
    exhaust_thrust_contribution: float = 0  # [-]
    exhaust_expansion_ratio: float = 0  # [-]
    turbine_maximum_temperature: float = 0  # [K]
    turbine_pressure_ratio: Optional[float] = None  # [-]
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
        if self.verbose and self.iterate:
            self.flow_check()
            self.check_exhaust_thrust_contribution()
            self.do_pressure_check()

    def flow_check(self, turbine_compare_flow_state: Optional[FlowState] = None):
        if turbine_compare_flow_state is None:
            turbine_compare_flow_state = self.default_turbine_flow_check_state
        turbine_inlet_state = self.turbine.inlet_flow_state
        # Line below ensures that warning after if statement is given multiple times if encountered multiple times, instead of only once
        warnings.simplefilter('always', UserWarning)
        if not turbine_compare_flow_state.almost_equal(turbine_inlet_state):
            if self.verbose:
                warnings.warn(
                    'FlowStates at the Turbine Inlet and the upstream Outlet FlowState are not equal after iteration')
                print('Param.         : \t     XX-out \t      TU-in')
                for key, gg_item in turbine_compare_flow_state.print_pretty_dict.items():
                    tu_item = turbine_inlet_state.print_pretty_dict[key]
                    print(f'{key: <15}: \t {gg_item: >10} \t {tu_item: >10}')

    @property
    def default_turbine_flow_check_state(self) -> FlowState:
        raise NotImplementedError

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
        if self.verbose:
            warnings.warn(f'\nExhaust thrust contribution initial value (used in program) different from final estimate:\n'
                          f'Used         : {self.exhaust_thrust_contribution * self.thrust * 1e-3:.3f} kN\n'
                          f'Estimate     : {self.turbine_exhaust.thrust * 1e-3:.3f} kN')

    @property
    def turbine_mass_flow(self):
        return self._iterative_turbine_mass_flow

    @property
    def turbine_mass_flow_initial_guess(self):
        return 0.0

    @property
    def turbine_inlet_temperature(self):
        return 0.0

    @property
    def turbine_inlet_pressure(self):
        return self.cooling_channel_section.outlet_pressure

    @property
    def turbine_inlet_flow_state(self):
        return FlowState(propellant_name=self.fuel_name,
                         temperature=self.turbine_inlet_temperature,
                         pressure=self.turbine_inlet_pressure,
                         mass_flow=self._iterative_turbine_mass_flow,
                         type='fuel',)

    @property
    def turbine(self):
        return Turbine(inlet_flow_state=self.turbine_inlet_flow_state,
                       pump_power_required=self.pump_power_required,
                       efficiency=self.turbine_efficiency,
                       pressure_ratio=self.turbine_pressure_ratio,
                       gas_heat_capacity_ratio=self.turbine_gas_heat_capacity_ratio,
                       gas_specific_heat_capacity=self.turbine_gas_specific_heat_capacity, )

    @property
    def turbine_exhaust(self):
        return SecondaryExhaust(inlet_flow_state=self.turbine.outlet_flow_state,
                                expansion_ratio=self.exhaust_expansion_ratio,
                                ambient_pressure=self.ambient_pressure,
                                gas_heat_capacity_ratio=self.turbine_gas_heat_capacity_ratio,
                                gas_molar_mass=self.turbine_gas_molar_mass, )

    @property
    def chamber_mass_flow(self):
        return (1 - self.exhaust_thrust_contribution) * self.base_mass_flow

    @property
    def total_mass_flow(self):
        return self.chamber_mass_flow + self.turbine_mass_flow

    @property
    def overall_specific_impulse(self):
        return self.thrust / self.total_mass_flow / g

    @property
    def secondary_specific_impulse(self):
        return self.turbine_exhaust.specific_impulse
