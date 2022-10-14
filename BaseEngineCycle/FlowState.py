from dataclasses import dataclass, field
from CoolProp.CoolProp import PropsSI
from functools import cached_property
from typing import Literal, Optional


@dataclass
class FlowState:
    """
    Keep track of the state of the (propellant) flow and access associated properties easily as attributes

    Flow is assumed to have no flow speed for simplicity: static = total conditions, see DynamicFlowState if flow speed
    is needed
    """
    propellant_name: str
    temperature: float  # [K]
    pressure: float  # [Pa]
    mass_flow: Optional[float]  # [kg/s]
    type: Literal['oxidizer', 'fuel', 'combusted']

    @property
    def print_pretty_dict(self):
        from collections import defaultdict
        fstrings = defaultdict(lambda: '', {'temperature': '.0f', 'pressure': '.3e', 'mass_flow': '.3e'})
        return {key: f'{item:{fstrings[key]}}' for key, item in vars(self).items()}

    @property
    def coolprop_name(self):
        p_name = self.propellant_name.upper()

        def matches_any(patterns: list[str, ...]):
            return any(pattern in p_name for pattern in patterns)

        if matches_any(['RP', 'ROCKETPROPELLANT']):
            return 'n-Dodecane'
        elif matches_any(['H2', 'HYDROGEN']):
            return 'Hydrogen'
        elif matches_any(['O2', 'OXYGEN', 'LOX']):
            return 'Oxygen'
        elif matches_any(['CH4', 'METHANE']):
            return 'Methane'
        else:
            raise ValueError('No matching coolprop_name was recognized for propellant_name')

    @cached_property
    def molar_mass(self):
        return PropsSI('MOLAR_MASS', self.coolprop_name)

    @property
    def state_inputs(self):
        return 'T', self.temperature, 'P', self.pressure, self.coolprop_name

    @property
    def specific_heat_capacity(self):
        return PropsSI('CPMASS', *self.state_inputs)

    @property
    def specific_heat_capacity_const_volume(self):
        return PropsSI('CVMASS', *self.state_inputs)

    @property
    def heat_capacity_ratio(self):
        return self.specific_heat_capacity / self.specific_heat_capacity_const_volume

    @property
    def density(self):
        if self.coolprop_name == 'n-Dodecane':
            # Known correction RP-1 generally 3-4% heavier than dodecane
            return PropsSI('DMASS', *self.state_inputs) * 1.04
        return PropsSI('DMASS', *self.state_inputs)

    @property
    def mass_specific_enthalpy(self):
        return PropsSI('HMASS', *self.state_inputs)

    @property
    def prandtl_number(self):
        return PropsSI('PRANDTL', *self.state_inputs)

    @property
    def conductivity(self):
        return PropsSI('L', *self.state_inputs)

    @property
    def dynamic_viscosity(self):
        return PropsSI('V', *self.state_inputs)

    @property
    def speed_of_sound(self):
        return PropsSI('A', *self.state_inputs)

    def propssi(self, string_input: str):
        return PropsSI(string_input, *self.state_inputs)

    def get_reynolds(self, linear_dimension: float, flow_speed: float):
        return self.density * flow_speed * linear_dimension / self.dynamic_viscosity

    def almost_equal(self, other: 'FlowState', margin: float = 1e-8) -> bool:
        """Checks if flowstates have the same fields but leaves some margin for floating point errors"""
        equals_list = []
        for own_val, other_val in zip(self.__dict__.values(), other.__dict__.values()):
            if type(own_val) == float:
                equals_list.append( abs(own_val - other_val) / own_val < margin)
            else:
                equals_list.append(own_val == other_val)
        return all(equals_list)


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
    _flow_speed: float = field(repr=False)

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

                self._static_temperature = self.static_temperature
                self._static_pressure = self.static_pressure
                if self.mach > 1:
                    raise ValueError('CoolingFlowState flow_speed is higher than Mach 1, increase the amount of '
                                     'channels or total flow area to decrease the flow speed')
            iterations += 1
            if iterations > self._max_iterations:
                break

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
        return self._flow_speed

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
    _flow_speed: float = field(init=False, repr=False, default=None)

    @property
    def flow_speed(self):
        return self.mass_flow / (self.density * self.total_flow_area)



@dataclass
class DefaultFlowState(FlowState):
    """When interface is derived from the default value provided, this will make linters shut up"""
    propellant_name: str = 'Default'
    temperature: float = 0  # [K]
    pressure: float = 0  # [Pa]
    mass_flow: float = 0  # [kg/s]
    type: str = 'Default'

    def coolprop_name(self):
        raise ValueError('Called a function or class that requires a FlowState without defining it')


@dataclass
class ManualFlowState(FlowState):
    """
    Override the FlowState with manually added flow properties, instead of properties derived from temp. and prressure

    Be careful, since changing the temperature and/or pressure obviously will not lead to an update of the given flow
    properties!
    """
    _molar_mass: Optional[float] = None
    _specific_heat_capacity: Optional[float] = None
    _specific_heat_capacity_const_volume: Optional[float] = None
    _heat_capacity_ratio: Optional[float] = None
    _density: Optional[float] = None
    _mass_specific_enthalpy: Optional[float] = None
    _prandtl_number: Optional[float] = None
    _conductivty: Optional[float] = None
    _dynamic_viscosity: Optional[float] = None
    _speed_of_sound: Optional[float] = None

    @cached_property
    def molar_mass(self):
        return self._molar_mass

    @cached_property
    def specific_heat_capacity(self):
        return self._specific_heat_capacity

    @cached_property
    def specific_heat_capacity_const_volume(self):
        return self._specific_heat_capacity_const_volume

    @property
    def heat_capacity_ratio(self):
        if self._heat_capacity_ratio is None:
            try:
                return self.specific_heat_capacity / self.specific_heat_capacity_const_volume
            except TypeError:
                return None
        else:
            return self._heat_capacity_ratio

    @cached_property
    def density(self):
        return self._density

    @cached_property
    def mass_specific_enthalpy(self):
        return self._mass_specific_enthalpy

    @cached_property
    def prandtl_number(self):
        return self._prandtl_number

    @cached_property
    def conductivity(self):
        return self._conductivty

    @cached_property
    def dynamic_viscosity(self):
        return self._dynamic_viscosity

    @cached_property
    def speed_of_sound(self):
        return self._speed_of_sound
