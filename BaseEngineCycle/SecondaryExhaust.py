from dataclasses import dataclass, field
from CoolProp.CoolProp import PropsSI
from irt import get_thrust_coefficient, get_pressure_ratio_fsolve, get_characteristic_velocity, get_throat_area
from scipy.constants import g, pi
from typing import Optional


@dataclass
class SecondaryExhaust:
    coolprop_name: str
    inlet_pressure: float  # [Pa]
    inlet_temperature: float  # [K]
    expansion_ratio: float  # [-]
    mass_flow: float  # [kg/s]
    ambient_pressure: Optional[float] = None  # [Pa]
    _ideal_expansion: bool = field(default=False, init=False)

    def __post_init__(self):
        if self.ambient_pressure is None:
            self._ideal_expansion = True

    @property
    def inlet_state_input(self):
        return ('T', self.inlet_temperature, 'P', self.inlet_pressure, self.coolprop_name)

    @property
    def specific_heat_capacity(self):
        return PropsSI('CPMASS', *self.inlet_state_input)

    @property
    def specific_heat_capacity_const_volume(self):
        return PropsSI('CVMASS', *self.inlet_state_input)

    @property
    def heat_capacity_ratio(self):
        return self.specific_heat_capacity / self.specific_heat_capacity_const_volume

    @property
    def molar_mass(self):
        return PropsSI('MOLAR_MASS', self.coolprop_name)

    @property
    def pressure_ratio(self):
        return get_pressure_ratio_fsolve(expansion_ratio=self.expansion_ratio,
                                         heat_capacity_ratio=self.heat_capacity_ratio)

    @property
    def characteristic_velocity(self):
        return get_characteristic_velocity(molar_mass=self.molar_mass,
                                           chamber_temperature=self.inlet_temperature,
                                           heat_capacity_ratio=self.heat_capacity_ratio)

    @property
    def thrust_coefficient(self):
        return get_thrust_coefficient(pressure_ratio=self.pressure_ratio,
                                      heat_capacity_ratio=self.heat_capacity_ratio,
                                      expansion_ratio=self.expansion_ratio,
                                      chamber_pressure=self.inlet_pressure,
                                      ambient_pressure=self.ambient_pressure,
                                      ideal_expansion=self._ideal_expansion)

    @property
    def equivalent_velocity(self):
        return self.characteristic_velocity * self.thrust_coefficient

    @property
    def specific_impulse(self):
        return self.equivalent_velocity / g

    @property
    def exit_pressure(self):
        return self.inlet_pressure / self.pressure_ratio

    @property
    def throat_area(self):
        return get_throat_area(molar_mass=self.molar_mass,
                               heat_capacity_ratio=self.heat_capacity_ratio,
                               chamber_temperature=self.inlet_temperature,
                               mass_flow=self.mass_flow,
                               chamber_pressure=self.inlet_pressure)

    @property
    def exit_area(self):
        return self.expansion_ratio * self.throat_area

    @property
    def thrust(self):
        if self._ideal_expansion:
            return self.equivalent_velocity * self.mass_flow
        else:
            return (self.equivalent_velocity * self.mass_flow
                    + ((self.exit_pressure - self.ambient_pressure) / self.mass_flow) * self.exit_area)


if __name__ == '__main__':
    test_exhaust = SecondaryExhaust(inlet_pressure=2e5,
                                    inlet_temperature=470.5,
                                    coolprop_name='Methane',
                                    mass_flow=.87,
                                    expansion_ratio=20)
    print(test_exhaust.specific_impulse)
