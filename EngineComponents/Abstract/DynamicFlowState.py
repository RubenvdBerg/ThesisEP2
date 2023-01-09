from dataclasses import dataclass, field
from functools import cached_property
from typing import Optional
from CoolProp.CoolProp import PropsSI
from EngineComponents.Abstract.FlowState import FlowState


@dataclass
class DynamicFlowState(FlowState):
    """
    Adds flow speed property to FlowState object: Differentiate between static and total conditions

    .pressure and .temperature attributes are replaced by the total properties of the flow
    """
    pressure: float = field(init=False, repr=False, default=None)
    temperature: float = field(init=False, repr=False, default=None)
    total_temperature: float
    total_pressure: float
    given_flow_speed: float = field(repr=False)

    _iteration_accuracy: float = 1e-3
    _max_iterations: float = 5
    _static_temperature: float = field(init=False, repr=False, default=None)
    _static_pressure: float = field(init=False, repr=False, default=None)
    verbose: bool = False
    _static_temperature_initial_guess: Optional[float] = None
    _static_pressure_initial_guess: Optional[float] = None

    def __post_init__(self):
        # Total state as initial guess for static state unless initial guess specified
        self._static_temperature = self.total_temperature if self._static_temperature_initial_guess is None else self._static_temperature_initial_guess
        self._static_pressure = self.total_pressure if self._static_pressure_initial_guess is None else self._static_pressure_initial_guess
        self.iterate()

    def iterate(self):
        """
        Density is dependent on the static pressure and static pressure is calculated by subtracting the
        dynamic pressure (which requires the density) from the total pressure. Thus; iteration is required.
        Same for cp, density, and static temperature.
        """
        if self.verbose:
            print('DynamicFlowStateIteration:')
        iterations = 0
        while (self.error_too_large(self._static_pressure, self.static_pressure)
               or self.error_too_large(self._static_temperature, self.static_temperature)):
            if self.verbose:
                print(f'Cur.:{self._static_temperature:.5e} K, {self._static_pressure:.6e} Pa\n'
                      f'Exp.:{self.static_temperature:.6e} K, {self.static_pressure:.6e} Pa\n'
                      f'Mach:{self.mach:.3f}')
            iterations += 1
            if iterations > self._max_iterations:
                break
            self._static_temperature = self.static_temperature
            self._static_pressure = self.static_pressure
            if self.mach > 1:
                print(f'Total state: {self.total_pressure * 1e-5:.1f} bar, {self.total_temperature:.1f} K')
                print(f'Static state: {self._static_pressure * 1e-5:.1f} bar, {self._static_temperature:.1f} K')
                raise ValueError(
                    'CoolingFlowState flow_speed is higher than Mach 1, decrease the mass flow or increase '
                    'the total flow area to decrease the flow speed')

    @property
    def state_inputs(self):
        """Make flow properties dependent on iteration variables, otherwise requesting static temp/pressure is an endless loop"""
        return 'T', self._static_temperature, 'P', self._static_pressure, self.coolprop_name

    def error_too_large(self, current: float, expected: float):
        error = abs((current - expected) / expected)
        return error > self._iteration_accuracy

    @cached_property
    def flow_speed(self):
        """Reroute flow_speed, which is required, so it can be overridden in child class CoolantFlowState"""
        return self.given_flow_speed

    @property
    def mach(self):
        return self.flow_speed / self.speed_of_sound

    @property
    def dynamic_temp(self):
        return .5 * self.flow_speed ** 2 / self.specific_heat_capacity

    @property
    def dynamic_pressure(self):
        return .5 * self.density * self.flow_speed ** 2

    @property
    def static_temperature(self):
        return self.total_temperature - self.dynamic_temp

    @property
    def static_pressure(self):
        return self.total_pressure - self.dynamic_pressure

    def get_reynolds(self, linear_dimension: float, flow_speed: float = None):
        if flow_speed is None:
            flow_speed = self.flow_speed
        return super().get_reynolds(flow_speed=flow_speed, linear_dimension=linear_dimension)

    @property
    def total_state_inputs(self):
        return 'T', self.total_temperature, 'P', self.total_pressure, self.coolprop_name

    @property
    def total_mass_specific_enthalpy(self):
        return PropsSI('H', *self.total_state_inputs)


@dataclass
class CoolantFlowState(DynamicFlowState):
    """Same as DynamicFlowState, but internally calculates flow speed from mass flux instead"""
    total_flow_area: float = 0
    given_flow_speed: float = field(init=False, repr=False, default=None)

    @property
    def flow_speed(self):
        return self.mass_flow / (self.density * self.total_flow_area)