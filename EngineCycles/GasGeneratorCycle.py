from dataclasses import dataclass, field
from typing import Optional
import warnings
from scipy import constants
from EngineCycles.Abstract.OpenCycle import OpenEngineCycle, OpenEngineCycle_DoubleTurbine
from EngineComponents.Base.Splitter import Splitter
from EngineComponents.Other.GasGenerator import GasGenerator
from EngineComponents.Abstract.Material import Material
from EngineComponents.Abstract.FlowState import FlowState, DefaultFlowState, ManualFlowState
from EngineFunctions.CEAFunctions import get_gas_generator_mmr, get_cea_dict_gg
from EngineFunctions.EmpiricalRelations import get_gas_generator_mmr_rp1

# Baseclass that can either inherit from single or double turbine OpenCycle (see next classes)
@dataclass
class GasGeneratorCycle_Mixin:
    # TODO: dataclass inheritance is stupid see EP-class
    gg_stay_time: float = 0  # [s]
    gg_structural_factor: float = 0  # [-]
    gg_material: Material = None
    gg_base_flow_state: Optional[FlowState] = None
    gg_is_frozen: Optional[bool] = None
    gg_pressure: Optional[float] = None  # [Pa]
    gg_mass_mixture_ratio: Optional[float] = None  # [-]
    _cea_gg_mmr: float = field(init=False, default=None)  # [-]

    def __post_init__(self):
        super().__post_init__()

    def set_initial_values(self):
        super().set_initial_values()
        if self.gg_is_frozen is None:
            self.gg_is_frozen = self.is_frozen
        if self.gg_pressure is None:
            self.gg_pressure = self.combustion_chamber_pressure
        if self.gg_mass_mixture_ratio is None:
            self.gg_mass_mixture_ratio = get_gas_generator_mmr(temperature_limit=self.turbine_maximum_temperature,
                                                               **self.cea_gg_kwargs)
            if self.gg_mass_mixture_ratio < 0.5 and 'RP' in self.fuel_name:
                self._cea_gg_mmr = self.gg_mass_mixture_ratio
                # CEA for fuel-rich RP1/LOX combustion is quite a ways off, this empirical relation is better
                self.gg_mass_mixture_ratio = get_gas_generator_mmr_rp1(temperature_limit=self.turbine_maximum_temperature)
        if self.gg_base_flow_state is None:
            # If an empirical mixture ratio is used (only happens if it's a better approximation than CEA)
            # this ratio is used in further calculation, EXCEPT for calculation of the gg_bas_flow_state below
            mixture_ratio = self.gg_mass_mixture_ratio if self._cea_gg_mmr is None else self._cea_gg_mmr
            cea_dict = get_cea_dict_gg(MR=mixture_ratio,
                                       frozen=1 if self.gg_is_frozen else 0,
                                       frozenAtThroat=1 if self.gg_is_frozen else 0,
                                       **self.cea_gg_kwargs)
            self.gg_base_flow_state = ManualFlowState(propellant_name='ExhaustGas',
                                                      temperature=cea_dict['T_C'],
                                                      pressure=self.gg_pressure,
                                                      mass_flow=None,
                                                      type='combusted',
                                                      _heat_capacity_ratio=cea_dict['y_cc'],
                                                      _specific_heat_capacity=cea_dict['cp_cc'],
                                                      _molar_mass=cea_dict['mm_cc'],
                                                      _density=cea_dict['rho_cc'])
        else:
            if self.gg_base_flow_state.pressure is None:
                self.gg_base_flow_state.pressure = self.gg_pressure

        self.check_gg_temp_and_pressure()

    @property
    def cea_gg_kwargs(self):
        return {'fuelName': self.fuel_name, 'oxName': self.oxidizer_name, 'Pc': self.gg_pressure, }

    @property
    def ideal_gas_densty(self):
        r = constants.gas_constant / self.gg_base_flow_state.molar_mass
        return self.gg_pressure / (r * self.gg_base_flow_state.temperature)

    def check_gg_temp_and_pressure(self):
        if self.gg_base_flow_state.temperature > self.turbine_maximum_temperature * 1.01:
            warnings.warn(
                f'The combustion temperature of the gas generator is higher than the maximum allowed turbine inlet '
                f'temperature [{self.turbine_maximum_temperature}]. The combustion temperature '
                f'[{self.gg_base_flow_state.temperature}] was either provided manually or calculated from a manually '
                f'provided mixture ratio [{self.gg_mass_mixture_ratio}]'
            )
        if self.gg_base_flow_state.pressure != self.gg_pressure:
            raise ValueError('Pressure provided through gg_base_flow_state must be the same as gg_pressure.')

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
        return self.turbine_mass_flow

    @property
    def gg_oxidizer_flow(self):
        return self.gg_mass_mixture_ratio / (self.gg_mass_mixture_ratio + 1) * self.gg_mass_flow

    @property
    def gg_fuel_flow(self):
        return 1 / (self.gg_mass_mixture_ratio + 1) * self.gg_mass_flow

    @property
    def main_fuel_flow(self):  # Override EngineCycle flows
        return self.chamber_fuel_flow + self.gg_fuel_flow

    @property
    def main_oxidizer_flow(self):  # Override EngineCycle flows
        return self.chamber_oxidizer_flow + self.gg_oxidizer_flow

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
                            stay_time=self.gg_stay_time,
                            base_flow_state=self.gg_base_flow_state,
                            outlet_temperature=self.turbine_maximum_temperature,
                            safety_factor=self.gg_structural_factor,
                            structure_material=self.gg_material,
                            )

    @property
    def oxidizer_pump_expected_pressure(self):
        p1 = self.combustion_chamber_pressure - self.injector.pressure_change
        p2 = self.gg_pressure
        return max(p1, p2)

    @property
    def gg_propellant_mass(self):
        return self.gg_mass_flow * self.burn_time * self.propellant_margin_factor

    @property
    def cc_propellant_mass(self):
        return self.chamber_mass_flow * self.burn_time * self.propellant_margin_factor

    @property
    def feed_system_mass(self):
        return super().feed_system_mass + self.gas_generator.mass

    @property
    def mass_kwak(self):
        return super().mass + self.gas_generator.mass

    @property
    def engine_dry_mass(self):
        return super().engine_dry_mass + self.gas_generator.mass


@dataclass
class GasGeneratorCycle(GasGeneratorCycle_Mixin, OpenEngineCycle):
    @property
    def turbine_mass_flow_initial_guess(self):
        """Initial guess based on verification engines. If no iteration is requested start at 0 to clearly show flows
        without any turbine requirements"""
        return .03 * self.base_mass_flow if self.iterate else 0


@dataclass
class GasGeneratorCycle_DoubleTurbine(GasGeneratorCycle_Mixin, OpenEngineCycle_DoubleTurbine):

    @property
    def fuel_turbine_mass_flow_initial_guess(self):
        return .015 * self.base_mass_flow if self.iterate else 0

    @property
    def oxidizer_turbine_mass_flow_initial_guess(self):
        return .015 * self.base_mass_flow if self.iterate else 0


@dataclass
class GasGeneratorCycle_DoubleTurbineSeries(GasGeneratorCycle_DoubleTurbine):
    """GasGeneratorCycle_DoubleTurbine assumes the turbines are in parallel and have separate exhausts. This
    configuration has turbines in series and a single exhaust. The fuel secondary exhaust and associated variables are
    removed."""

    # Variables no longer needed
    fuel_secondary_specific_impulse_quality_factor: float = field(init=False, repr=False)
    fuel_exhaust_expansion_ratio: Optional[float] = field(init=False, repr=False)
    fuel_exhaust_exit_pressure_forced: Optional[float] = field(init=False, repr=False)

    # Properties no longer needed
    @property
    def turbine_splitter(self):
        raise NotImplementedError

    @property
    def fuel_secondary_exhaust(self):
        raise NotImplementedError

    # Adjusted inlet flows
    @property
    def fuel_turbine_inlet_flow_state(self):
        return self.gas_generator.outlet_flow_state

    @property
    def oxidizer_turbine_inlet_flow_state(self):
        return self.fuel_turbine.outlet_flow_state

    # Adjusted exhaust thrust
    @property
    def exhaust_total_thrust(self):
        return self.oxidizer_secondary_exhaust.thrust

    @property
    def engine_dry_mass(self):
        return self.feed_system_mass + self.thrust_chamber.mass + self.oxidizer_secondary_exhaust.mass


if __name__ == '__main__':
    from EngineArguments import arguments as args

    print(GasGeneratorCycle(**args.desgin_arguments, **args.base_arguments, **args.gg_arguments))
