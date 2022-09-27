from dataclasses import dataclass
from CoolProp.CoolProp import PropsSI
from functools import cached_property
from typing import Literal, Optional


@dataclass
class FlowState:
    propellant_name: str
    temperature: float  # [K]
    pressure: float  # [Pa]
    mass_flow: Optional[float]  # [kg/s]
    type: Literal['oxidizer', 'fuel', 'combusted']

    @property
    def print_pretty_dict(self):
        from collections import defaultdict
        fstrings = defaultdict(lambda: '', {'temperature': '.0f', 'pressure': '.3e', 'mass_flow': '.3e'})
        return {key: f'{item:{fstrings[key]}}' for key, item in self.__dict__.items()}

    @cached_property
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

    def get_reynolds(self, flow_speed: float, linear_dimension: float):
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
class DefaultFlowState(FlowState):
    propellant_name: str = 'Default'
    temperature: float = 0  # [K]
    pressure: float = 0  # [Pa]
    mass_flow: float = 0  # [kg/s]
    type: str = 'Default'

    def coolprop_name(self):
        raise ValueError('Called a function or class that requires a FlowState without defining it')


@dataclass
class ManualFlowState(FlowState):
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