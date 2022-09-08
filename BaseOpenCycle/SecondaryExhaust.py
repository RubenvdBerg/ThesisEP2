from dataclasses import dataclass, field
from CoolProp.CoolProp import PropsSI
from irt import get_thrust_coefficient, get_pressure_ratio_fsolve, get_characteristic_velocity, get_throat_area
from scipy.constants import g, pi, gas_constant
from typing import Optional
from functools import cached_property
from BaseEngineCycle.FlowComponent import FlowComponent


@dataclass
class SecondaryExhaust(FlowComponent):
    expansion_ratio: float = 0  # [-]
    gas_heat_capacity_ratio: Optional[float] = None  # [-]
    gas_molar_mass: Optional[float] = None  # [kg/mol]
    ambient_pressure: Optional[float] = None  # [Pa]
    _ideal_expansion: bool = field(default=False, init=False)

    def __post_init__(self):
        if self.ambient_pressure is None:
            self._ideal_expansion = True

    @cached_property
    def _gas_heat_capacity_ratio(self):
        if self.gas_heat_capacity_ratio is None:
            return self.inlet_flow_state.heat_capacity_ratio
        else:
            return self.gas_heat_capacity_ratio

    @cached_property
    def gas_molar_mass(self):
        if self.gas_molar_mass is None:
            return self.inlet_flow_state.molar_mass
        else:
            return self.gas_molar_mass

    @property
    def pressure_ratio(self):
        return get_pressure_ratio_fsolve(expansion_ratio=self.expansion_ratio,
                                         heat_capacity_ratio=self._gas_heat_capacity_ratio)

    @property
    def characteristic_velocity(self):
        return get_characteristic_velocity(molar_mass=self.gas_molar_mass,
                                           chamber_temperature=self.inlet_temperature,
                                           heat_capacity_ratio=self._gas_heat_capacity_ratio)

    @property
    def thrust_coefficient(self):
        return get_thrust_coefficient(pressure_ratio=self.pressure_ratio,
                                      heat_capacity_ratio=self._gas_heat_capacity_ratio,
                                      expansion_ratio=self.expansion_ratio,
                                      chamber_pressure=self.inlet_pressure,
                                      ambient_pressure=self.ambient_pressure,
                                      ideal_expansion=self._ideal_expansion)

    @property
    def temperature_ratio(self):
        return self.pressure_ratio**((self._gas_heat_capacity_ratio - 1) / self._gas_heat_capacity_ratio)

    @property
    def equivalent_velocity(self):
        return self.characteristic_velocity * self.thrust_coefficient

    @property
    def specific_impulse(self):
        return self.equivalent_velocity / g

    @property
    def pressure_change(self):
        return self.inlet_pressure / self.pressure_ratio - self.inlet_pressure

    @property
    def temperature_change(self):
        return self.inlet_temperature / self.temperature_ratio - self.inlet_temperature

    @property
    def throat_area(self):
        return get_throat_area(molar_mass=self.inlet_flow_state.molar_mass,
                               heat_capacity_ratio=self._gas_heat_capacity_ratio,
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
    from BaseEngineCycle.FlowState import FlowState
    f1 = FlowState(propellant_name='Methane', temperature=470.5, pressure=2e5, mass_flow=.87, type='fuel')
    test_exhaust = SecondaryExhaust(inlet_flow_state=f1,
                                    expansion_ratio=20)
    print(test_exhaust.specific_impulse)
