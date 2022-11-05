from dataclasses import dataclass, field
from typing import Optional
import CoolProp.CoolProp as CoolProp
from EngineCycles.BaseEngineCycle.FlowComponent import FlowComponent


@dataclass
class CoolingChannelSection(FlowComponent):
    total_heat_transfer: float = 0  # [W]
    maximum_outlet_temperature: float = 0  # [K]
    pressure_drop: Optional[float] = None  # [Pa]
    combustion_chamber_pressure: Optional[float] = None  # [Pa]
    _pressure_drop_ratio: float = field(init=False, default=.15)  # [-]
    _wall_temp_factor: float = .9
    verbose: bool = True
    _instance_created: bool = False
    _is_temp_calc_needed: bool = True

    def __post_init__(self):
        self.coolprop_name = self.inlet_flow_state.coolprop_name
        CoolProp.set_reference_state(self.coolprop_name, 'NBP')

        # Ugly fix to prevent recursion, see EngineCycle.injector_inlet_flow_states for explanation
        CoolingChannelSection.__instance_created = True

    @property
    def pressure_change(self):
        if self.pressure_drop is None and self.combustion_chamber_pressure is None:
            raise ValueError(
                'One of [pressure_drop] and [combustion_chamber_pressure](to estimate the pressure_drop) must be given'
            )
        if self.pressure_drop is None:
            # Humble 1995 p.209 suggest pressure drop to be 10% - 20% of chamber pressure
            return -self.combustion_chamber_pressure * self._pressure_drop_ratio
        else:
            return -self.pressure_drop

    @property
    def increase_mass_specific_enthalpy(self):
        return self.total_heat_transfer / self.inlet_flow_state.mass_flow

    @property
    def outlet_mass_specific_enthalpy(self):
        return self.inlet_flow_state.mass_specific_enthalpy + self.increase_mass_specific_enthalpy

    @property
    def ideal_outlet_temperature(self):
        return CoolProp.PropsSI('T',
                                'H', self.outlet_mass_specific_enthalpy,
                                'P', self.inlet_pressure,
                                self.coolprop_name)

    @property
    def temperature_change(self):
        if self._is_temp_calc_needed:
            return self.ideal_outlet_temperature - self.inlet_temperature
        else:
            return 0

    @property
    def min_mass_flow(self):
        enthalpy_in = self.inlet_flow_state.mass_specific_enthalpy
        enthalpy_max = CoolProp.PropsSI('H',
                                        'T', self.maximum_outlet_temperature,
                                        'P', self.inlet_pressure,
                                        self.coolprop_name)
        delta_h_max = enthalpy_max - enthalpy_in
        q_tot = self.total_heat_transfer
        return abs(q_tot / delta_h_max)
