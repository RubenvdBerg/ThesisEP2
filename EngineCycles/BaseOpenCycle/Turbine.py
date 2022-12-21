from dataclasses import dataclass
from typing import Optional
from EngineCycles.Functions.BaseFunctions import only_one_none
from EngineCycles.BaseEngineCycle.FlowComponent import FlowComponent
from functools import cached_property


@dataclass
class BaseTurbine:
    """Finds the missing value based on the basic turbine power equation
    (eq. 10-18, Rocket Propulsion Elements 9th Ed., Sutton)
    """
    power_required: float = 0  # [W]
    efficiency: float = 0  # [-]
    specific_heat_capacity: float = 0  # [J/kg]
    heat_capacity_ratio: float = 0  # [-]

    # Two out of these three must be given
    pressure_ratio: Optional[float] = None  # [-]
    mass_flow: Optional[float] = None  # [kg/s]
    inlet_temperature: Optional[float] = None  # [K]

    def __post_init__(self):
        if not only_one_none(self.pressure_ratio, self.mass_flow, self.inlet_temperature):
            raise ValueError('Exactly two out of pressure_ratio, mass_flow, and inlet_temperature must be provided')

    @property
    def mass_flow_required(self):
        if self.mass_flow is None:
            return (self.power_required
                    / (self.efficiency * self.specific_heat_capacity * self.inlet_temperature
                       * (1 - self.pressure_ratio ** ((1 - self.heat_capacity_ratio) / self.heat_capacity_ratio))))
        else:
            raise ValueError(self.error_message('mass_flow'))

    @property
    def inlet_temperature_required(self):
        if self.inlet_temperature is None:
            return (self.power_required
                    / (self.efficiency * self.specific_heat_capacity * self.mass_flow
                       * (1 - self.pressure_ratio ** ((1 - self.heat_capacity_ratio) / self.heat_capacity_ratio))))
        else:
            raise ValueError(self.error_message('inlet_temperature'))

    @property
    def pressure_ratio_required(self):
        if self.inlet_temperature is None:
            return (self.power_required
                    / (self.efficiency * self.specific_heat_capacity * self.mass_flow
                       * (1 - self.pressure_ratio ** ((1 - self.heat_capacity_ratio) / self.heat_capacity_ratio))))
        else:
            raise ValueError(self.error_message('pressure_ratio'))

    def error_message(self, attribute_name):
        return f'Turbine {attribute_name} is already given [{getattr(self, attribute_name):}], also estimating ' \
               f'{attribute_name}_required makes no sense'


@dataclass
class Turbine(FlowComponent):
    power_required: float = 0  # [W]
    efficiency: float = 0  # [-]
    pressure_ratio: Optional[float] = None  # [-]
    outlet_pressure_forced: Optional[float] = None  # [Pa]
    # gas_heat_capacity_ratio: Optional[float] = None  # [-]
    # gas_specific_heat_capacity: Optional[float] = None  # [J/kg]

    def __post_init__(self):
        self.resolve_pressure_ratio_choice()

    def resolve_pressure_ratio_choice(self):
        if not ((self.pressure_ratio is None) ^ (self.outlet_pressure_forced is None)):
            raise ValueError('Both or neither the pressure_ratio and the outlet_pressure of the turbine are provided. '
                             'Provide one and only one')
        elif self.pressure_ratio is None:
            self.pressure_ratio = self.inlet_pressure / self.outlet_pressure_forced

    # @cached_property
    # def drive_gas_heat_capacity_ratio(self):
    #     if self.gas_heat_capacity_ratio is None:
    #         return self.inlet_flow_state.heat_capacity_ratio
    #     else:
    #         return self.gas_heat_capacity_ratio
    #
    # @cached_property
    # def drive_gas_specific_heat_capacity(self):
    #     if self.gas_specific_heat_capacity is None:
    #         return self.inlet_flow_state.specific_heat_capacity
    #     else:
    #         return self.gas_specific_heat_capacity

    @property
    def mass_flow_required(self):
        return BaseTurbine(
            power_required=self.power_required,
            efficiency=self.efficiency,
            specific_heat_capacity=self.inlet_flow_state.specific_heat_capacity,
            heat_capacity_ratio=self.inlet_flow_state.heat_capacity_ratio,
            pressure_ratio=self.pressure_ratio,
            inlet_temperature=self.inlet_temperature
        ).mass_flow_required

    @property
    def temperature_change(self):
        if self.mass_flow == 0:
            return 0
        dt = -1 * self.power_required / (self.mass_flow * self.inlet_flow_state.specific_heat_capacity)
        if (-1*dt) > self.inlet_temperature:
            raise ValueError('Turbine temperature decrease is larger than turbine inlet temperature')
        return dt

    @property
    def pressure_change(self):
        return self.inlet_pressure / self.pressure_ratio - self.inlet_pressure