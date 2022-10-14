from BaseEngineCycle.HeatExchangerOMECA import RectangularOMECAHeatExchanger
import BaseEngineCycle.EmpiricalRelations as empirical
from dataclasses import dataclass
from typing import Optional





@dataclass
class DeniesHeatExchanger(RectangularOMECAHeatExchanger):
    constant_mach: Optional[float] = 3

    def set_local_mach(self):
        if self.constant_mach:
            self.section_local_mach = self.constant_mach
        else:
            super().set_local_mach()

    @property
    def section_reduction_factor(self):
        return 1

    @property
    def section_hot_gas_convective_heat_transfer_coefficient(self):
        # Changed stagnation_temp input from adiabatic wall to combustion temp
        h_g = empirical.get_hot_gas_convective_heat_transfer_coefficient(
            mass_flow=self.coolant_inlet_flow_state.mass_flow,
            local_diameter=self.section_radius * 2,
            dynamic_viscosity=self.combustion_chamber_flow_state.dynamic_viscosity,
            specific_heat_capacity=self.combustion_chamber_flow_state.specific_heat_capacity,
            prandtl_number=self.combustion_chamber_flow_state.prandtl_number,
            stagnation_temp=self.combustion_temp,
            film_temp=self.section_hot_gas_film_temp,
            mode=self.hot_gas_convective_heat_transfer_coefficient_mode,
            local_mach=self.section_local_mach,
            wall_temp=self.section_hot_side_wall_temp,
            heat_capacity_ratio=self.combustion_chamber_flow_state.heat_capacity_ratio,
            throat_radius_of_curvature=self.thrust_chamber_section.nozzle.conv_throat_long_radius,
            combustion_chamber_pressure=self.combustion_chamber_flow_state.pressure,
            characteristic_velocity=self.characteristic_velocity,
            nozzle_throat_diameter=self.thrust_chamber_section.get_radius(0) * 2,
        )
        return h_g / 0.026 * 0.0195
