from BaseEngineCycle.Structure import Structure
from dataclasses import dataclass


@dataclass
class GasGenerator(Structure):
    pressure: float  # [Pa]
    gas_constant: float  # [J/(kg*K)]
    stay_time: float  # [s]
    turbine_mass_flow: float  # [kg/s]
    turbine_temp_limit: float  # [K]

    # TODO: Attributes below don't actually do anything, placeholder for more complex gas generator with
    #  turbine temperature check using CEA. Do not remove mass_mixture_ratio either as it is used as an attribute
    mass_mixture_ratio: float  # [-]

    # oxidizer: Propellant
    # fuel: Propellant
    # fuel_or_oxidizer_rich: Literal['fuel_rich', 'oxidizer_rich']

    @property
    def mass(self):  # Overwrite Structure.mass()
        return (self.safety_factor * 3 / 2 * self.material_density / self.yield_strength
                * self.stay_time * self.pressure / self.gas_density * self.mass_flow)

    @property
    def gas_density(self):
        return self.pressure / (self.gas_constant * self.turbine_temp_limit)

    @property
    def mass_flow(self):
        return self.turbine_mass_flow


