from EngineCycles.BaseEngineCycle.EngineCycle import EngineCycle
from EngineCycles.BaseOpenCycle.Turbine import Turbine
from EngineCycles.BaseOpenCycle.SecondaryExhaust import SecondaryExhaust
from EngineCycles.BaseEngineCycle.FlowState import FlowState
from EngineCycles.BaseEngineCycle.Splitter import Splitter
from dataclasses import dataclass, field, replace
from scipy.constants import g
from math import isclose
from typing import Optional
import warnings
from functools import cached_property


# Abstract class that has the attributes GasGenerator and OpenExpander share, the turbine_mass_flow property needs to be
# used somewhere to increase the main flows, such that this iterative process works.
# OpenEngineCycle does NOT work on its own!
@dataclass
class OpenEngineCycle(EngineCycle):
    turbine_efficiency: float = 0  # [-]
    turbopump_specific_power: float = 0  # [W/kg]
    turbine_maximum_temperature: float = 0  # [K]
    secondary_specific_impulse_correction_factor: float = 1
    exhaust_expansion_ratio: Optional[float] = None  # [-]
    exhaust_exit_pressure: Optional[float] = None  # [Pa]
    turbine_pressure_ratio: Optional[float] = None  # [-]
    turbine_outlet_pressure_forced: Optional[float] = None  # [Pa]

    # Iteration attribute, not required at init
    _iterative_turbine_mass_flow: float = field(init=False, repr=False, default=0.00001)  # [kg/s]
    _exhaust_thrust_contribution: float = field(init=False, repr=False, default=.01)  # [-]

    # Override from EngineCycle, are assigned later, not required at init
    fuel_pump_specific_power: float = field(init=False, repr=False, default=None)
    oxidizer_pump_specific_power: float = field(init=False, repr=False, default=None)

    def __post_init__(self):
        super().__post_init__()
        if self.iterate:
            self.iterate_mass_flow()

    def iterate_mass_flow(self):
        if self.verbose:
            print('Iterate Turbine Mass Flow')
        self._iterative_turbine_mass_flow = self.turbine_mass_flow_initial_guess
        while self.turbine_flow_error_larger_than_accuracy():
            if self.verbose:
                print(f'Actual:  {self._iterative_turbine_mass_flow:.5f} kg/s')
                print(f'Required:{self.turbine.mass_flow_required:.5f} kg/s\n')
            self._iterative_turbine_mass_flow = self.turbine.mass_flow_required
            self._exhaust_thrust_contribution = self.secondary_exhaust.thrust / self.thrust

        if self.verbose:
            print(f'Turbine Mass Flow Set\n')

    @cached_property
    def turbine_mass_flow_initial_guess(self):
        return 0.0

    def set_initial_values(self):
        super().set_initial_values()
        # Combination of turbine, oxidizer-pump and fuel-pump seen as a single turbopump with single specific power
        self.fuel_pump_specific_power = self.oxidizer_pump_specific_power = self.turbopump_specific_power

    def turbine_flow_error_larger_than_accuracy(self):
        error = abs(self.turbine.mass_flow_required - self.turbine_mass_flow)
        margin = self.turbine.mass_flow_required * self.iteration_accuracy
        return error > margin

    @property
    def chamber_thrust(self):
        return (1 - self._exhaust_thrust_contribution) * self.thrust

    @property
    def turbine_mass_flow(self):
        return self._iterative_turbine_mass_flow

    @property
    def turbine_inlet_flow_state(self):
        return self.cooling_channel_section.outlet_flow_state

    @property
    def turbine(self):
        return Turbine(inlet_flow_state=self.turbine_inlet_flow_state,
                       power_required=self.pumps_power_required,
                       efficiency=self.turbine_efficiency,
                       pressure_ratio=self.turbine_pressure_ratio,
                       outlet_pressure_forced=self.turbine_outlet_pressure_forced, )

    @property
    def secondary_exhaust(self):
        return SecondaryExhaust(inlet_flow_state=self.turbine.outlet_flow_state,
                                expansion_ratio=self.exhaust_expansion_ratio,
                                exit_pressure=self.exhaust_exit_pressure,
                                ambient_pressure=self.ambient_pressure,
                                specific_impulse_correction_factor=self.secondary_specific_impulse_correction_factor, )

    @property
    def secondary_specific_impulse(self):
        return self.secondary_exhaust.specific_impulse


@dataclass
class OpenEngineCycle_DoubleTurbine(EngineCycle):
    turbine_maximum_temperature: float = 0  # [K]
    turbopump_specific_power: float = 0
    oxidizer_turbine_efficiency: float = 0  # [-]
    fuel_turbine_efficiency: float = 0  # [-]
    oxidizer_shaft_mechanical_efficiency: float = 1  # [-]
    fuel_shaft_mechanical_efficiency: float = 1  # [-]

    oxidizer_secondary_specific_impulse_correction_factor: float = 1
    fuel_secondary_specific_impulse_correction_factor: float = 1

    oxidizer_exhaust_expansion_ratio: Optional[float] = None  # [-]
    oxidizer_exhaust_exit_pressure_forced: Optional[float] = None  # [Pa]
    oxidizer_turbine_pressure_ratio: Optional[float] = None  # [-]
    oxidizer_turbine_outlet_pressure_forced: Optional[float] = None  # [Pa]
    fuel_exhaust_exit_pressure_forced: Optional[float] = None  # [Pa]
    fuel_exhaust_expansion_ratio: Optional[float] = None  # [-]
    fuel_turbine_pressure_ratio: Optional[float] = None  # [-]
    fuel_turbine_outlet_pressure_forced: Optional[float] = None  # [Pa]

    # Iteration attribute, not required at init
    _iterative_oxidizer_turbine_mass_flow: float = field(init=False, repr=False, default=0.00001)  # [kg/s]
    _iterative_fuel_turbine_mass_flow: float = field(init=False, repr=False, default=0.00001)  # [kg/s]
    _iterative_turbine_mass_flow: float = field(init=False, repr=False, default=0.0000001)  # [kg/s]
    _exhaust_thrust_contribution: float = field(init=False, repr=False, default=.01)  # [-]

    # Override from EngineCycle, are assigned later, not required at init
    fuel_pump_specific_power: float = field(init=False, repr=False, default=None)
    oxidizer_pump_specific_power: float = field(init=False, repr=False, default=None)

    def __post_init__(self):
        super().__post_init__()
        if self.iterate:
            self.iterate_mass_flow()

    def iterate_mass_flow(self):
        if self.verbose:
            print('Iterate Turbine Mass Flow')
        self._iterative_oxidizer_turbine_mass_flow = self.oxidizer_turbine_mass_flow_initial_guess
        self._iterative_fuel_turbine_mass_flow = self.fuel_turbine_mass_flow_initial_guess
        while self.turbine_flow_error_larger_than_accuracy():
            req_ox = self.oxidizer_turbine.mass_flow_required
            req_fu = self.fuel_turbine.mass_flow_required
            if self.verbose:
                print(f'Actual:  {self.turbine_mass_flow:.5f} kg/s')
                print(f'Required:{req_ox + req_fu:.5f} kg/s\n')
            self._iterative_oxidizer_turbine_mass_flow = req_ox
            self._iterative_fuel_turbine_mass_flow = req_fu
            self._exhaust_thrust_contribution = (self.oxidizer_secondary_exhaust.thrust
                                                 + self.fuel_secondary_exhaust.thrust) / self.thrust

        if self.verbose:
            print(f'Turbine Mass Flow Set\n')

    @cached_property
    def oxidizer_turbine_mass_flow_initial_guess(self):
        return 0.0

    @cached_property
    def fuel_turbine_mass_flow_initial_guess(self):
        return 0.0

    def set_initial_values(self):
        super().set_initial_values()
        self.fuel_pump_specific_power = self.oxidizer_pump_specific_power = self.turbopump_specific_power

    def turbine_flow_error_larger_than_accuracy(self):
        required = self.oxidizer_turbine.mass_flow_required + self.fuel_turbine.mass_flow_required
        error = abs(required - self.turbine_mass_flow)
        margin = required * self.iteration_accuracy
        return error > margin

    @property
    def chamber_thrust(self):
        return (1 - self._exhaust_thrust_contribution) * self.thrust

    @property
    def turbine_mass_flow(self):
        return self._iterative_fuel_turbine_mass_flow + self._iterative_oxidizer_turbine_mass_flow

    @property
    def turbine_inlet_flow_state(self):
        return self.cooling_channel_section.outlet_flow_state

    @property
    def turbine_splitter(self):
        return Splitter(inlet_flow_state=self.turbine_inlet_flow_state,
                        required_outlet_mass_flows=(self._iterative_oxidizer_turbine_mass_flow,),
                        outlet_flow_names=('oxidizer_turbine', 'fuel_turbine'))

    @property
    def fuel_pumps_power_required(self):
        return self.fuel_pump.power_required / self.shaft_mechanical_efficiency

    @property
    def fuel_turbine(self):
        return Turbine(inlet_flow_state=self.turbine_splitter.outlet_flow_states['fuel_turbine'],
                       power_required=self.fuel_pumps_power_required,
                       efficiency=self.fuel_turbine_efficiency,
                       pressure_ratio=self.fuel_turbine_pressure_ratio,
                       outlet_pressure_forced=self.fuel_turbine_outlet_pressure_forced, )

    @property
    def fuel_secondary_exhaust(self):
        return SecondaryExhaust(
            inlet_flow_state=self.fuel_turbine.outlet_flow_state,
            expansion_ratio=self.fuel_exhaust_expansion_ratio,
            exit_pressure=self.fuel_exhaust_exit_pressure_forced,
            ambient_pressure=self.ambient_pressure,
            specific_impulse_correction_factor=self.fuel_secondary_specific_impulse_correction_factor,
        )

    @property
    def oxidizer_pumps_power_required(self):
        return self.oxidizer_pump.power_required / self.shaft_mechanical_efficiency

    @property
    def oxidizer_turbine(self):
        return Turbine(inlet_flow_state=self.turbine_splitter.outlet_flow_states['oxidizer_turbine'],
                       power_required=self.oxidizer_pumps_power_required,
                       efficiency=self.oxidizer_turbine_efficiency,
                       pressure_ratio=self.oxidizer_turbine_pressure_ratio,
                       outlet_pressure_forced=self.oxidizer_turbine_outlet_pressure_forced, )

    @property
    def oxidizer_secondary_exhaust(self):
        return SecondaryExhaust(
            inlet_flow_state=self.oxidizer_turbine.outlet_flow_state,
            expansion_ratio=self.oxidizer_exhaust_expansion_ratio,
            exit_pressure=self.oxidizer_exhaust_exit_pressure_forced,
            ambient_pressure=self.ambient_pressure,
            specific_impulse_correction_factor=self.oxidizer_secondary_specific_impulse_correction_factor,
        )
