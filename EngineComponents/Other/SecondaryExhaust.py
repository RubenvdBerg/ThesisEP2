from dataclasses import dataclass, field
from EngineFunctions.IRTFunctions import get_thrust_coefficient, get_pressure_ratio_fsolve, \
    get_characteristic_velocity, get_throat_area, get_expansion_ratio_from_p_ratio, is_choked
from scipy.constants import g
from math import radians, pi, sqrt, tan
from typing import Optional
from functools import cached_property
from EngineComponents.Abstract.FlowComponent import FlowComponent
from EngineComponents.Abstract.StructuralComponent import StructuralComponent
import warnings


@dataclass
class SecondaryExhaust(FlowComponent, StructuralComponent):
    specific_impulse_quality_factor: float = 1
    exit_pressure: Optional[float] = None
    expansion_ratio: Optional[float] = None  # [-]

    ambient_pressure: Optional[float] = None  # [Pa]
    pressure_ratio: Optional[float] = field(init=False, repr=False)  # [-]
    flow_heat_capacity_ratio: Optional[float] = field(init=False, repr=False)  # [-]
    flow_molar_mass: Optional[float] = field(init=False, repr=False) # [kg/mol]
    _divergence_angle: float = radians(15)

    def __post_init__(self):
        self.set_flow_properties()
        self.resolve_expansion_choice()
        self.check_choked_flow()

    def set_flow_properties(self):
        """Give the assumed replacement of RP1, i.e. Dodecane, similar temperature limits as RP1."""
        self.flow_molar_mass = self.inlet_flow_state.molar_mass
        try:
            self.flow_heat_capacity_ratio = self.inlet_flow_state.heat_capacity_ratio
        except ValueError:
            if 'RP' in self.inlet_flow_state.propellant_name and 230 < self.inlet_temperature < 264:
                # Set to simple constant
                self.flow_heat_capacity_ratio = 1.25
            else:
                raise

    def resolve_expansion_choice(self):
        if not ((self.expansion_ratio is None) ^ (self.exit_pressure is None)):
            raise ValueError('Neither or both the secondary expansion_ratio and secondary exit_pressure are given. '
                             'Provide one and only one')
        elif self.expansion_ratio is None:
            self.pressure_ratio = self.inlet_pressure / self.exit_pressure
            self.expansion_ratio = get_expansion_ratio_from_p_ratio(pressure_ratio=self.pressure_ratio,
                                                                    heat_capacity_ratio=self.flow_heat_capacity_ratio)
        else:
            self.pressure_ratio = get_pressure_ratio_fsolve(expansion_ratio=self.expansion_ratio,
                                                            heat_capacity_ratio=self.flow_heat_capacity_ratio)
            self.exit_pressure = self.inlet_pressure / self.pressure_ratio

    def check_choked_flow(self):
        if not is_choked(self.pressure_ratio, self.flow_heat_capacity_ratio):
            warnings.warn('Flow in a secondary exhaust is not choked.'
                          ' Ideal rocket theory, which is used to estimate the specific impulse, is invalid!')

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
                                      ambient_pressure=self.ambient_pressure,
                                      summerfield_criterion=False)

    @property
    def temperature_ratio(self):
        return self.pressure_ratio ** ((self.flow_heat_capacity_ratio - 1) / self.flow_heat_capacity_ratio)

    @property
    def equivalent_velocity(self):
        return self.characteristic_velocity * self.thrust_coefficient * self.specific_impulse_quality_factor

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
        return get_throat_area(molar_mass=self.flow_molar_mass,
                               heat_capacity_ratio=self.flow_heat_capacity_ratio,
                               chamber_temperature=self.inlet_temperature,
                               mass_flow=self.mass_flow,
                               chamber_pressure=self.inlet_pressure)

    @property
    def exit_area(self):
        return self.expansion_ratio * self.throat_area

    @property
    def throat_radius(self):
        return sqrt(self.throat_area / pi)

    @property
    def exit_radius(self):
        return sqrt(self.exit_area / pi)

    @property
    def length(self):
        return (self.exit_radius - self.throat_radius) / tan(self._divergence_angle)

    @cached_property
    def surface_area(self):
        # Curved surface area of frustrum with throat-/exit radius and divergence angle
        return (self.exit_area - self.throat_area) / tan(self._divergence_angle)

    @cached_property
    def wall_thickness(self):
        calculated_thickness = self.inlet_pressure * self.exit_radius / self.structure_material.yield_strength
        if calculated_thickness < self.structure_material.minimal_thickness:
            warnings.warn(
                'Calculated thickness of secondary exhaust is lower than minimum thickness of material. The latter is used.')
            return self.structure_material.minimal_thickness
        return calculated_thickness

    @cached_property
    def mass(self):
        return self.safety_factor * self.surface_area * self.wall_thickness * self.structure_material.density

    @property
    def thrust(self):
        """Reminder: pressure term already accounted for in equivalent_velocity through thrust coefficient"""
        return self.equivalent_velocity * self.mass_flow
