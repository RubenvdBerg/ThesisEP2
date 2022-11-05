from dataclasses import dataclass, field
from scipy.constants import gas_constant
from typing import Optional
from EngineCycles.BaseOpenCycle.OpenCycle import OpenEngineCycle
from EngineCycles.BaseEngineCycle.Splitter import Splitter
from EngineCycles.BaseEngineCycle.FlowState import FlowState
from EngineCycles.GasGeneratorCycle.GGComponents import GasGenerator


@dataclass
class GasGeneratorCycle(OpenEngineCycle):
    # TODO: dataclass inheritance is stupid see EP-class
    gg_gas_specific_gas_constant: float = 0  # [J/(kg*K)]
    gg_mass_mixture_ratio: float = 0  # [-]
    gg_stay_time: float = 0  # [s]
    gg_structural_factor: float = 0  # [-]
    gg_material_density: float = 0  # [kg/m3]
    gg_yield_strength: float = 0  # [Pa]
    gg_pressure: Optional[float] = None  # [Pa]
    _is_temp_calc_needed: bool = field(init=False, repr=False, default=False)

    # No longer required, but calculated from gg_specific_gas_constant
    turbine_gas_molar_mass: float = field(init=False)

    def __post_init__(self):
        self.turbine_gas_molar_mass = gas_constant / self.gg_gas_specific_gas_constant
        super().__post_init__()

    @property
    def default_turbine_flow_check_state(self) -> FlowState:
        return self.gas_generator.outlet_flow_state

    @property
    def turbine_mass_flow_initial_guess(self):
        """Initial guess based on verification engines. If no iteration is requested start at 0 to clearly show flows
        without any turbine requirements"""
        return .03 * self.base_mass_flow if self.iterate else 0

    @property
    def turbine_inlet_flow_state(self):
        """Turbine operates with gas generator exhaust at maximum allowable temperature."""
        return self.gas_generator.outlet_flow_state

    @property
    def post_oxidizer_pump_splitter(self):
        """Splits the flow into the required chamber oxidizer flow and 'extra' flow, which will be equal to the required
        gas generator oxidizer flow after iteration"""
        return Splitter(inlet_flow_state=self.oxidizer_pump.outlet_flow_state,
                        required_outlet_mass_flows=(self.chamber_oxidizer_flow,),
                        outlet_flow_names=('main', 'gg'))

    @property
    def post_fuel_pump_splitter(self):
        """Splits the flow into the required chamber fuel flow and 'extra' flow, which will be equal to the required gas
        generator fuel flow after iteration"""
        return Splitter(inlet_flow_state=self.fuel_pump.outlet_flow_state,
                        required_outlet_mass_flows=(self.chamber_fuel_flow,),
                        outlet_flow_names=('main', 'gg'))

    @property
    def gg_mass_flow(self):  # Must be equal
        return self._iterative_turbine_mass_flow

    @property
    def main_fuel_flow(self):  # Override EngineCycle flows
        return (self.chamber_fuel_flow
                + 1 / (self.gg_mass_mixture_ratio + 1) * self.gg_mass_flow)

    @property
    def main_oxidizer_flow(self):  # Override EngineCycle flows
        return (self.chamber_oxidizer_flow
                + self.gg_mass_mixture_ratio / (self.gg_mass_mixture_ratio + 1) * self.gg_mass_flow)

    @property
    def cooling_inlet_flow_state(self):
        """Adjusting default EngineCycle connection (from fuel pump to cooling) to account for splitter inbetween"""
        return self.post_fuel_pump_splitter.outlet_flow_state_main

    @property
    def injector_inlet_flow_states(self):
        return self.cooling_channel_section.outlet_flow_state, self.post_oxidizer_pump_splitter.outlet_flow_state_main

    @property
    def gas_generator(self):
        return GasGenerator(oxidizer_inlet_flow_state=self.post_oxidizer_pump_splitter.outlet_flow_state_gg,
                            fuel_inlet_flow_state=self.post_fuel_pump_splitter.outlet_flow_state_gg,
                            mass_mixture_ratio=self.gg_mass_mixture_ratio,
                            gas_constant=self.gg_gas_specific_gas_constant,
                            stay_time=self.gg_stay_time,
                            pressure=self.combustion_chamber_pressure if self.gg_pressure is None else self.gg_pressure,
                            turbine_temp_limit=self.turbine_maximum_temperature,
                            safety_factor=self.gg_structural_factor,
                            material_density=self.gg_material_density,
                            yield_strength=self.gg_yield_strength,
                            )

    @property
    def gg_propellant_mass(self):
        return self.gg_mass_flow * self.burn_time * self.propellant_margin_factor

    @property
    def cc_propellant_mass(self):
        return self.chamber_mass_flow * self.burn_time * self.propellant_margin_factor

    @property
    def feed_system_mass(self):
        return self.gas_generator.mass + self.pumps_mass

    @property
    def mass(self):
        return (self.cc_propellant_mass
                + self.gg_propellant_mass
                + self.feed_system_mass
                + self.tanks_mass
                + self.pressurant.mass)


if __name__ == '__main__':
    import arguments as args

    print(GasGeneratorCycle(**args.desgin_arguments, **args.base_arguments, **args.gg_arguments))
