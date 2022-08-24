from dataclasses import dataclass
from CoolProp.CoolProp import PropsSI


@dataclass
class FlowState:
    propellant_name: str
    temperature: float  # [K]
    pressure: float  # [Pa]
    mass_flow: float  # [kg/s]

    @cached_property
    def coolprop_name(self):
        p_name = self.propellant_name.upper()
        if 'RP' in p_name:
            return 'n-Dodecane'
        if 'H2' in p_name:
            return 'Hydrogen'
        if 'O2' in p_name or 'OX' in p_name:
            return 'Oxygen'
        if 'CH4' in p_name:
            return 'Methane'

    @property
    def molar_mass(self):
        return PropsSI('MOLAR_MASS',self.coolprop_name)

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
        return PropsSI('DMASS', *self.state_inputs)

    @property
    def mass_specific_enthalpy(self):
        return PropsSI('HMASS', *self.state_inputs)

    def propssi(self, string_input: str):
        return PropsSI(string_input, *self.state_inputs)