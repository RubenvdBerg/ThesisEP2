import math

from EngineCycles.OpenExpanderCycle.OECycle import OpenExpanderCycle
from dataclasses import dataclass

@dataclass
class SE21D_V3(OpenExpanderCycle):
    """See 'Sippel et al. 2003 - Studies on Expander Bleed Cycle Engines for Launchers' for this rocket's configuration.
        """

    @property
    def coolant_base_state(self):
        cbs = super().coolant_base_state
        cbs.temperature = 32.375
        return cbs

    @property
    def turbine_mass_flow_initial_guess(self):
        return self.base_mass_flow * .015

    @property
    def fuel_pump(self):
        fuel_pump = super().fuel_pump
        fuel_pump._temperature_change = 8.44
        return fuel_pump

    @property
    def oxidizer_pump(self):
        oxidizer_pump = super().oxidizer_pump
        oxidizer_pump._temperature_change = 92.983 - 90
        return oxidizer_pump

    @property
    def turbine(self):
        turbine = super().turbine
        turbine._temperature_change = 369.677 - 506.452
        return turbine

    @property
    def secondary_fuel_pump(self):
        pump = super().secondary_fuel_pump
        pump._temperature_change = 32.375 - 29.44
        return pump


@dataclass
class SE21D_Exact_V3(SE21D_V3):
    cooling_section_pressure_drop: float = 3.36e6

    @property
    def fuel_pump_expected_pressure(self):
        return 8.749e6

    @property
    def oxidizer_pump_expected_pressure(self):
        return 8.749e6

    @property
    def secondary_fuel_pump_expected_pressure(self):
        return 12.109e6

    @property
    def throat_area(self):
        return math.pi * .286**2

    def set_heat_transfer(self):
        self.total_heat_transfer = 111.037e6