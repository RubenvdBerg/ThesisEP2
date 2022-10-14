from dataclasses import dataclass, field
from typing import Optional
from math import pi
import BaseEngineCycle.EmpiricalRelations as empirical
from BaseEngineCycle.HeatExchanger import HeatExchanger


@dataclass
class OMECAHeatExchanger(HeatExchanger):
    maximum_wall_temp: float = field(init=False, repr=False)
    chamber_wall_conductivity: float = 0
    chamber_wall_thickness: float = 0
    characteristic_velocity: Optional[float] = 0
    hot_gas_convective_heat_transfer_coefficient_mode: str = 'ModifiedBartz'
    coolant_heat_transfer_coefficient_mode: str = 'SiederTate'

    # Iterative section variables, not required at init
    section_hot_side_wall_temp: float = field(init=False, repr=False)
    section_cold_side_wall_temp: float = field(init=False, repr=False)
    section_length_to_start_channel: float = field(init=False, repr=False, default=.01)

    def __post_init__(self):
        # Set initial values for additional attributes at start of coolant channel
        self.section_hot_side_wall_temp = 400 if self.is_counter_flow else self.coolant_inlet_flow_state.temperature
        self.section_cold_side_wall_temp = self.coolant_inlet_flow_state.temperature
        super().__post_init__()

    @property
    def hot_wall_temp(self):
        """Overwrite usage of constant maximum_wall_temp to a variable hot_side_wall_temp."""
        return self.section_hot_side_wall_temp

    def determine_heat_flux(self):
        """Overwrite heat flux step in iteration process with wall temp sub-iteration."""
        self.determine_heat_flux_and_iterate_wall_temps()

    def determine_heat_flux_and_iterate_wall_temps(self):
        """Loop until hot and cold side wall temperatures converge."""
        th_w = self.chamber_wall_thickness
        k_w = self.chamber_wall_conductivity
        t_b = self.section_coolant_bulk_temp
        self.set_local_mach()
        self.set_hot_gas_temps()
        t_ad = self.section_hot_gas_adiabatic_wall_temp

        iterations = 0
        while True:
            h_g = self.section_hot_gas_total_heat_transfer_coefficient
            h_c = self.section_coolant_convective_heat_transfer_coefficient
            q = (t_ad - t_b) / (1 / h_g + th_w / k_w + 1 / h_c)
            self.section_heat_flux = q

            new_wall_temp_hot = t_ad - q / h_g if h_g else self.section_hot_side_wall_temp
            new_wall_temp_cold = t_b + q / h_c if h_c else self.section_cold_side_wall_temp

            # Check if iteration is done
            iterations += 1
            error_hot = self.error_is_small(self.section_hot_side_wall_temp, new_wall_temp_hot)
            error_cold = self.error_is_small(self.section_cold_side_wall_temp, new_wall_temp_cold)
            if iterations > self.max_iterations or (error_hot and error_cold):
                self.max_hot_side_wall_temp = max(self.max_hot_side_wall_temp, self.section_hot_side_wall_temp)
                break

            # Update variables
            self.section_hot_side_wall_temp = float(new_wall_temp_hot)
            self.section_cold_side_wall_temp = float(new_wall_temp_cold)

    def set_next_state(self):
        super().set_next_state()
        self.section_length_to_start_channel += self.section_wall_length

    @property
    def section_hot_gas_convective_heat_transfer_coefficient(self):
        """Overwrite parent property with version that allows for mode selection."""
        h_g = empirical.get_hot_gas_convective_heat_transfer_coefficient(
            mass_flow=self.coolant_inlet_flow_state.mass_flow,
            local_diameter=self.section_radius * 2,
            dynamic_viscosity=self.combustion_chamber_flow_state.dynamic_viscosity,
            specific_heat_capacity=self.combustion_chamber_flow_state.specific_heat_capacity,
            prandtl_number=self.combustion_chamber_flow_state.prandtl_number,
            stagnation_temp=self.section_hot_gas_adiabatic_wall_temp,
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
        return h_g

    @property
    def section_coolant_convective_heat_transfer_coefficient(self):
        """Calculate cool coeff. required for wall temp estimation."""
        k_c = self.section_coolant_state.conductivity
        Dh = self.section_coolant_channel_diameter
        re = self.section_coolant_state.get_reynolds(linear_dimension=Dh)
        pr = self.section_coolant_state.prandtl_number
        t_b = self.section_coolant_bulk_temp
        t_w = self.section_cold_side_wall_temp
        e = self.coolant_channel_roughness_height

        h_c = empirical.get_coolant_convective_heat_transfer_coeff(
            k_c, Dh, re, pr, t_b, t_w,
            mode=self.coolant_heat_transfer_coefficient_mode,
            length_to_start_channel=self.section_length_to_start_channel,
        )
        f_rough = empirical.get_roughness_correction(pr, re, e, Dh)
        return h_c * f_rough

    # Properties below only needed for data write, not for calculations
    @property
    def section_coolant_heat_flux(self):
        return self.section_coolant_convective_heat_transfer_coefficient * (self.section_cold_side_wall_temp
                                                                            - self.section_coolant_bulk_temp)

    def init_data(self):
        """Add to parent data keys"""
        super().init_data()
        temps = 'Temperature [K]'
        coeff = 'Heat-Transfer Coefficient [W/(K*m2]'
        flux = 'Heat Flux [W/m2]'
        self.data[temps] = self.data[temps] | {'Cold SideWall': [], 'Hot SideWall': [], }
        self.data[coeff] = self.data[coeff] | {'Coolant': [],}
        self.data[flux] = self.data[flux] | {'Cold Side': [], 'Hot Side': [],}

    def write_data(self):
        """Append new child parameters to data"""
        super().write_data()
        self.data['Heat Flux [W/m2]']['Hot Side'].append(self.section_hot_gas_total_heat_flux)
        self.data['Heat Flux [W/m2]']['Cold Side'].append(self.section_coolant_heat_flux)
        self.data['Temperature [K]']['Hot SideWall'].append(self.section_hot_side_wall_temp)
        self.data['Temperature [K]']['Cold SideWall'].append(self.section_cold_side_wall_temp)
        self.data['Heat-Transfer Coefficient [W/(K*m2]']['Coolant'].append(
            self.section_coolant_convective_heat_transfer_coefficient)

@dataclass
class RectangularOMECAHeatExchanger(OMECAHeatExchanger):
    coolant_channel_fin_thickness: float = 1e-3
    coolant_channel_height_input: float = 1e-3

    @property
    def section_coolant_channel_height(self):
        return self.coolant_channel_height_input

    @property
    def section_coolant_channel_width(self):
        section_circumference = self.section_radius * pi * 2
        return section_circumference / self.number_of_coolant_channels - self.coolant_channel_fin_thickness

    @property
    def section_coolant_channel_area(self):
        return self.section_coolant_channel_height * self.section_coolant_channel_width

    @property
    def section_coolant_channel_diameter(self):
        """Calculate equivalent hydraulic diameter for rectangular channel."""
        h = self.section_coolant_channel_height
        w = self.section_coolant_channel_width
        return 2 * h * w / (h + w)

    def get_fin_correction(self, convective_coeff: float) -> float:
        th_f = self.coolant_channel_fin_thickness
        h_ch = self.section_coolant_channel_height
        w_ch = self.section_coolant_channel_width
        k_w = self.chamber_wall_conductivity
        return empirical.get_fin_correction(convective_coeff, th_f, h_ch, w_ch, k_w)

    @property
    def section_coolant_convective_heat_transfer_coefficient(self):
        h_c = super().section_coolant_convective_heat_transfer_coefficient
        f_fin = self.get_fin_correction(h_c)
        return h_c * f_fin

    def init_data(self):
        super().init_data()
        new_dict = {'Channel Geometry [mm]': {'Height': [],
                                              'Width': [],
                                              'Diameter': [], }}
        # Merge dicts and overwrite with second dict
        self.data = self.data | new_dict

    def write_data(self):
        super().write_data()
        self.data['Channel Geometry [mm]']['Height'].append(self.section_coolant_channel_height * 1e3)
        self.data['Channel Geometry [mm]']['Width'].append(self.section_coolant_channel_width * 1e3)
