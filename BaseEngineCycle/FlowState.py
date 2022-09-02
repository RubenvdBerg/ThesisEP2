from dataclasses import dataclass
from CoolProp.CoolProp import PropsSI
from functools import cached_property
from typing import Literal


@dataclass
class FlowState:
    propellant_name: str
    temperature: float  # [K]
    pressure: float  # [Pa]
    mass_flow: float  # [kg/s]
    type: Literal['oxidizer', 'fuel']

    @cached_property
    def coolprop_name(self):
        p_name = self.propellant_name.upper()
        if 'RP' in p_name:
            return 'n-Dodecane'
        elif 'H2' in p_name:
            return 'Hydrogen'
        elif 'O2' in p_name or 'OX' in p_name:
            return 'Oxygen'
        elif 'CH4' in p_name:
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

    def propssi(self, string_input: str):
        return PropsSI(string_input, *self.state_inputs)


@dataclass
class DefaultFlowState(FlowState):
    propellant_name: str = 'Default'
    temperature: float = 0  # [K]
    pressure: float = 0  # [Pa]
    mass_flow: float = 0  # [kg/s]
    type: str = 'Default'

    def coolprop_name(self):
        raise ValueError('Called a function or class that requires a FlowState without defining it')
