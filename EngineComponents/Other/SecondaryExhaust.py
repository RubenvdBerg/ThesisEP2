from dataclasses import dataclass, field
from EngineFunctions.IRTFunctions import get_thrust_coefficient, get_pressure_ratio_fsolve, \
    get_characteristic_velocity, get_throat_area, get_expansion_ratio_from_p_ratio, is_choked
from scipy.constants import g
from typing import Optional
from EngineComponents.Abstract.FlowComponent import FlowComponent
import warnings


@dataclass
class SecondaryExhaust(FlowComponent):
    specific_impulse_correction_factor: float = 1
    exit_pressure: Optional[float] = None
    expansion_ratio: Optional[float] = None  # [-]

    ambient_pressure: Optional[float] = None  # [Pa]
    pressure_ratio:Optional[float] = field(init=False, repr=False)  # [-]
    flow_heat_capacity_ratio: Optional[float] = field(init=False, repr=False)  # [-]
    flow_molar_mass: Optional[float] = field(init=False, repr=False)  # [kg/mol]

    def __post_init__(self):
        self.set_flow_properties()
        self.resolve_expansion_choice()
        self.check_choked_flow()
        # self.pressure_check()

    def set_flow_properties(self):
        """Some magic to give the assumed replacement of RP1, i.e. Dodecane, similar limits."""
        try:
            self.flow_heat_capacity_ratio = self.inlet_flow_state.heat_capacity_ratio
        except ValueError:
            if 'RP' in self.inlet_flow_state.propellant_name and 230 < self.inlet_temperature < 264:
                # Set to simple constant
                self.flow_heat_capacity_ratio = 1.25
            else:
                raise
        self.flow_molar_mass = self.inlet_flow_state.molar_mass

    def resolve_expansion_choice(self):
        if not ((self.expansion_ratio is None) ^ (self.exit_pressure is None)):
            raise ValueError('Neither or both the secondary expansion_ratio and secondary exit_pressure are given. Provide one and only one')
        elif self.expansion_ratio is None:
            self.pressure_ratio = self.inlet_pressure / self.exit_pressure
            self.expansion_ratio = get_expansion_ratio_from_p_ratio(pressure_ratio=self.pressure_ratio,
                                                                    heat_capacity_ratio=self.flow_heat_capacity_ratio)
        else:
            self.pressure_ratio = get_pressure_ratio_fsolve(expansion_ratio=self.expansion_ratio,
                                                            heat_capacity_ratio=self.flow_heat_capacity_ratio)

    def check_choked_flow(self):
        if not is_choked(self.pressure_ratio, self.flow_heat_capacity_ratio):
            warnings.warn('Flow in the secondary exhaust is not choked. Ideal rocket theory used to estimate the specific impulse is invalid!')

    # def pressure_check(self):
    #     if self.ambient_pressure is not None:
    #         if self.ambient_pressure > self.inlet_pressure:
    #             raise ValueError('Turbine outlet pressure lower than ambient pressure. Decrease the turbine pressure ratio or directly increase the turbine outlet pressure.')
    #         elif self.ambient_pressure > self.outlet_pressure:
    #             raise ValueError('SecondaryExhaust exit pressure is lower than the ambient pressure. Decrease the expansion ratio or directly increase the exit pressure.')

    @property
    def characteristic_velocity(self):
        return get_characteristic_velocity(molar_mass=self.flow_molar_mass,
                                           chamber_temperature=self.inlet_temperature,
                                           heat_capacity_ratio=self.flow_heat_capacity_ratio)

    @property
    def thrust_coefficient(self):
        return get_thrust_coefficient(pressure_ratio=self.pressure_ratio,
                                      heat_capacity_ratio=self.flow_heat_capacity_ratio,
                                      expansion_ratio=self.expansion_ratio,
                                      chamber_pressure=self.inlet_pressure,
                                      ambient_pressure=self.ambient_pressure, )

    @property
    def temperature_ratio(self):
        return self.pressure_ratio**((self.flow_heat_capacity_ratio - 1) / self.flow_heat_capacity_ratio)

    @property
    def equivalent_velocity(self):
        return self.characteristic_velocity * self.thrust_coefficient * self.specific_impulse_correction_factor

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
                               heat_capacity_ratio=self.flow_heat_capacity_ratio,
                               chamber_temperature=self.inlet_temperature,
                               mass_flow=self.mass_flow,
                               chamber_pressure=self.inlet_pressure)

    @property
    def exit_area(self):
        return self.expansion_ratio * self.throat_area

    @property
    def thrust(self):
        """Reminder: pressure term already accounted for in equivalent_velocity through pressure term"""
        return self.equivalent_velocity * self.mass_flow


if __name__ == '__main__':
    from EngineComponents.Abstract.FlowState import FlowState
    f1 = FlowState(propellant_name='Methane', temperature=470.5, pressure=2e5, mass_flow=.87, type='fuel')
    test_exhaust = SecondaryExhaust(inlet_flow_state=f1,
                                    expansion_ratio=20)
    print(test_exhaust.specific_impulse)
