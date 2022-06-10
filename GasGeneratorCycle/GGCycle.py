from scipy.constants import g

from BaseEngineCycle.OpenCycle import OpenEngineCycle
from GasGeneratorCycle.GGComponents import GasGenerator
from BaseEngineCycle.Turbine import Turbine
from dataclasses import dataclass, field


@dataclass
class GasGeneratorCycle(OpenEngineCycle):
    # TODO: dataclasses is dumb so requires defaults for all attributes if super class has defaults
    turbine_maximum_temperature: float = 0  # [K]
    gg_gas_gas_constant: float = 0  # [J/(kg*K)]
    gg_mass_mixture_ratio: float = 0  # [-]
    gg_stay_time: float = 0  # [s]
    gg_structural_factor: float = 0  # [-]
    gg_material_density: float = 0  # [kg/m3]
    gg_yield_strength: float = 0  # [Pa]

    # Override from OpenEngineCycle, is assigned later, not required at init
    turbine_inlet_temperature: float = field(init=False, repr=False, default=None)

    def __post_init__(self):
        super().__post_init__()
        # Different attribute name solely to make it explicit gas generator temperature is limited by the turbine
        self.turbine_inlet_temperature = self.turbine_maximum_temperature

    @property
    def gg_mass_flow(self):
        # Must be equal
        return self.turbine_mass_flow

    # Rewrite Isp to use total mass flow
    @property
    def simple_specific_impulse(self):
        return self.thrust / self.total_mass_flow / g

    @property
    def base_fuel_flow(self):
        return 1 / (self.mass_mixture_ratio + 1) * self.mass_flow

    @property
    def base_oxidizer_flow(self):
        return self.mass_mixture_ratio / (self.mass_mixture_ratio + 1) * self.mass_flow

    @property
    def fuel_flow(self):  # Override EngineCycle flows
        return (1 / (self.mass_mixture_ratio + 1) * self.chamber_mass_flow
                + 1 / (self.gg_mass_mixture_ratio + 1) * self.gg_mass_flow)

    @property
    def oxidizer_flow(self):  # Override EngineCycle flows
        return (self.mass_mixture_ratio / (self.mass_mixture_ratio + 1) * self.chamber_mass_flow
                + self.gg_mass_mixture_ratio / (self.gg_mass_mixture_ratio + 1) * self.gg_mass_flow)

    @property
    def gas_generator(self):
        return GasGenerator(mass_mixture_ratio=self.gg_mass_mixture_ratio, pressure=self.combustion_chamber_pressure,
                            gas_constant=self.gg_gas_gas_constant, stay_time=self.gg_stay_time,
                            turbine_mass_flow=self.turbine.mass_flow_required, safety_factor=self.gg_structural_factor,
                            material_density=self.gg_material_density, yield_strength=self.gg_yield_strength,
                            turbine_temp_limit=self.turbine_maximum_temperature)

    @property
    def pumps_mass(self):
        return self.fuel_pump.mass + self.oxidizer_pump.mass

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
