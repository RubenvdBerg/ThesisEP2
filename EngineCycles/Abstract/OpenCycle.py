from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineComponents.Abstract.Material import Material
from EngineComponents.Other.Turbine import Turbine
from EngineComponents.Other.SecondaryExhaust import SecondaryExhaust
from EngineComponents.Base.Splitter import Splitter
from dataclasses import dataclass, field
from typing import Optional
from functools import cached_property


# Abstract class that has the attributes GasGenerator and OpenExpander share, the turbine_mass_flow property needs to be
# used somewhere to increase the main flows, such that this iterative process works.
# OpenEngineCycle does NOT work on its own!
@dataclass
class OpenEngineCycle(EngineCycle):
    turbine_efficiency: float = 0  # [-]
    turbopump_specific_power: float = 0  # [W/kg]
    turbine_maximum_temperature: float = 0  # [K]
    exhaust_material: Material = None
    exhaust_safety_factor: float = 0  # [-]
    secondary_specific_impulse_quality_factor: Optional[float] = None
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

    def iterate_flow(self):
        self._iterative_turbine_mass_flow = self.turbine_mass_flow_initial_guess
        while self.turbine_flow_error_larger_than_accuracy():
            self.print_verbose_iteration_message()
            self._iterative_turbine_mass_flow = self.turbine.mass_flow_required
            self._exhaust_thrust_contribution = self.secondary_exhaust.thrust / self.thrust

    @property
    def verbose_iteration_name(self):
        return 'Turbine Mass Flow'

    @property
    def verbose_iteration_actual(self):
        return self._iterative_turbine_mass_flow

    @property
    def verbose_iteration_required(self):
        return self.turbine.mass_flow_required

    @cached_property
    def turbine_mass_flow_initial_guess(self):
        return 0.0

    def set_initial_values(self):
        super().set_initial_values()
        # Combination of turbine, oxidizer-pump and fuel-pump seen as a single turbopump with single specific power
        self.fuel_pump_specific_power = self.oxidizer_pump_specific_power = self.turbopump_specific_power
        if self.secondary_specific_impulse_quality_factor is None:
            self.secondary_specific_impulse_quality_factor = self.specific_impulse_quality_factor

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
                                specific_impulse_quality_factor=self.secondary_specific_impulse_quality_factor,
                                safety_factor=self.exhaust_safety_factor,
                                structure_material=self.exhaust_material)

    @property
    def secondary_specific_impulse(self):
        return self.secondary_exhaust.specific_impulse

    @property
    def feed_system_mass(self):
        return super().feed_system_mass + self.turbine.mass

    @property
    def engine_dry_mass(self):
        return super().engine_dry_mass + self.secondary_exhaust.mass

    @property
    def turbine_propellant_mass(self):
        return self.turbine_mass_flow * self.burn_time * self.propellant_margin_factor


@dataclass
class OpenEngineCycle_DoubleTurbine(EngineCycle):
    turbine_maximum_temperature: float = 0  # [K]
    turbopump_specific_power: float = 0
    oxidizer_turbine_efficiency: float = 0  # [-]
    fuel_turbine_efficiency: float = 0  # [-]
    exhaust_material: Material = None
    exhaust_safety_factor: float = 0  # [-]

    oxidizer_secondary_specific_impulse_quality_factor: float = 1
    fuel_secondary_specific_impulse_quality_factor: float = 1

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

    def iterate_flow(self):
        self._iterative_oxidizer_turbine_mass_flow = self.oxidizer_turbine_mass_flow_initial_guess
        self._iterative_fuel_turbine_mass_flow = self.fuel_turbine_mass_flow_initial_guess
        while self.turbine_flow_error_larger_than_accuracy():
            self.print_verbose_iteration_message()
            self._iterative_fuel_turbine_mass_flow = self.fuel_turbine.mass_flow_required
            self._iterative_oxidizer_turbine_mass_flow = self.oxidizer_turbine.mass_flow_required
            self._exhaust_thrust_contribution = self.exhaust_total_thrust / self.thrust

    @property
    def verbose_iteration_name(self):
        return 'Turbine Mass Flow'

    @property
    def verbose_iteration_actual(self):
        return self.turbine_mass_flow

    @property
    def verbose_iteration_required(self):
        return self.fuel_turbine.mass_flow_required + self.oxidizer_turbine.mass_flow_required

    @property
    def exhaust_total_thrust(self):
        return self.oxidizer_secondary_exhaust.thrust + self.fuel_secondary_exhaust.thrust

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
    def fuel_turbine_inlet_flow_state(self):
        return self.turbine_splitter.outlet_flow_states['fuel_turbine']

    @property
    def oxidizer_turbine_inlet_flow_state(self):
        return self.turbine_splitter.outlet_flow_states['oxidizer_turbine']

    @property
    def turbine_splitter(self):
        return Splitter(inlet_flow_state=self.turbine_inlet_flow_state,
                        required_outlet_mass_flows=(self._iterative_oxidizer_turbine_mass_flow,),
                        outlet_flow_names=('oxidizer_turbine', 'fuel_turbine'))

    @property
    def fuel_turbine(self):
        return Turbine(inlet_flow_state=self.fuel_turbine_inlet_flow_state,
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
            specific_impulse_quality_factor=self.fuel_secondary_specific_impulse_quality_factor,
            structure_material=self.exhaust_material,
            safety_factor=self.exhaust_safety_factor,
        )

    @property
    def oxidizer_turbine(self):
        return Turbine(inlet_flow_state=self.oxidizer_turbine_inlet_flow_state,
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
            specific_impulse_quality_factor=self.oxidizer_secondary_specific_impulse_quality_factor,
            structure_material=self.exhaust_material,
            safety_factor=self.exhaust_safety_factor,
        )

    @property
    def feed_system_mass(self):
        return super().feed_system_mass + self.fuel_turbine.mass + self.oxidizer_turbine.mass

    @property
    def engine_dry_mass(self):
        return super().engine_dry_mass + self.fuel_secondary_exhaust.mass + self.oxidizer_secondary_exhaust.mass

    @property
    def turbine_propellant_mass(self):
        return self.turbine_mass_flow * self.burn_time * self.propellant_margin_factor
