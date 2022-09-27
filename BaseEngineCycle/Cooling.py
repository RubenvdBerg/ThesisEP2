import math
from dataclasses import dataclass, field, replace
from typing import Optional, ClassVar
import CoolProp.CoolProp as CoolProp
from BaseEngineCycle.FlowComponent import FlowComponent
from BaseEngineCycle.HeatTransferSection import HeatTransferSection


@dataclass
class CoolingChannelSection(FlowComponent):
    heat_transfer_section: HeatTransferSection = object()
    _total_heat_transfer: Optional[float] = None
    pressure_drop: Optional[float] = None  # [Pa]
    combustion_chamber_pressure: Optional[float] = None  # [Pa]
    _pressure_drop_ratio: float = field(init=False, default=.15)  # [-]
    verbose: bool = True
    _instance_created: ClassVar[bool] = False

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
    def total_heat_transfer(self):
        if self._total_heat_transfer is None:
            return self.heat_transfer_section.total_heat_transfer
        else:
            return self._total_heat_transfer

    # @property
    # def default_inlet_temperature(self):
    #     _temp_in_dict = {'Hydrogen': 20.25, 'Oxygen': 90.15, 'n-Dodecane': 293.15, 'Methane': 111.66}
    #     return _temp_in_dict[self.coolprop_name]

    @property
    def increase_mass_specific_enthalpy(self):
        return self.total_heat_transfer / self.inlet_flow_state.mass_flow

    @property
    def outlet_mass_specific_enthalpy(self):
        return self.inlet_flow_state.mass_specific_enthalpy + self.increase_mass_specific_enthalpy

    @property
    def ideal_outlet_temperature(self):
        return CoolProp.PropsSI('T', 'H', self.outlet_mass_specific_enthalpy, 'P', self.outlet_pressure,
                                self.coolprop_name)

    @property
    def temperature_change(self):
        return self.ideal_outlet_temperature - self.inlet_temperature

    @property
    def bulk_temperature(self):
        return (self.inlet_temperature + self.outlet_temperature) / 2

    def throat_wall_temperature(self, wall_temperature: float, wall_thickness: float = .005, wall_conductivity: float = 350, coolant_channel_diameter:float = .005, number_of_channels: float = 140):
        """Calculates the maximum wall temperature (assumed to be at the throat) to be compared to initial assumption
        T_b = Coolant Bulk Temperature [K]
        q_a_g = Convective Heat Flux Hot-Gas Side [W/m2]
        q_r = Radiative Heat Flux Hot-Gas Side [W/m2]
        h_a_g = Convective Heat Transfer Coefficient Hot-Gas Side [W/(m2 K)]
        h_a_c = Convective Heat Transfer Coefficient Coolant Side [W/(m2 K)]
        L = Wall Thickness [m]
        k = Wall Material Conductivity [W/(m K)]
        """
        T_b = self.bulk_temperature
        h_a_g = self.heat_transfer_section.get_convective_heat_transfer_coefficient(distance_from_throat=0)
        T_w_ad = self.heat_transfer_section.get_adiabatic_wall_temp(distance_from_throat=0)
        q_a_g = h_a_g * (T_w_ad - wall_temperature)
        q_r = self.heat_transfer_section.radiative_heat_transfer_factor * q_a_g
        D = coolant_channel_diameter
        A = math.pi / 4 * D ** 2
        mass_flux = self.mass_flow / A
        bulk_pressure = (self.inlet_pressure + self.outlet_pressure) / 2
        bulk_state = replace(self.inlet_flow_state,
                             pressure=bulk_pressure,
                             temperature=self.bulk_temperature,)
        k_c = bulk_state.conductivity
        pr = bulk_state.prandtl_number
        u = mass_flux / number_of_channels / bulk_state.density
        a = bulk_state.speed_of_sound
        re = bulk_state.get_reynolds(flow_speed=u,linear_dimension=D)
        h_a_c = 0.025 * k_c / D * re ** .8 * pr ** .4 * (self.bulk_temperature / wall_temperature) ** .55
        L = wall_thickness
        k = wall_conductivity
        # q = (T_w_ad - T_b + q_r / h_a_g) / (1 / h_a_g + L / k + 1 / h_a_c)
        q = (T_w_ad - T_b + q_r / h_a_g) / (1 / h_a_g + 1 / h_a_c)
        q2 = q_a_g + q_r
        tw_g = T_w_ad - (q - q_r) / h_a_g
        tw_g2 = T_w_ad - (q2 - q_r) / h_a_g
        tw_c = T_b - q / h_a_c
        tw_c2 = T_b - q2 / h_a_c
        q_a_c = h_a_c * (wall_temperature - T_b)
        return T_w_ad - q / h_a_g

    def coolant_heat_transfer_coefficient(self, diameter: float = 0.01, wall_temperature: float = 0):
        """Returns the Sieder-Tate Heat Transfer Coefficient"""
        D = diameter
        A = math.pi / 4 * D ** 2
        mass_flux = self.mass_flow / A
        bulk_pressure = (self.inlet_pressure + self.outlet_pressure) / 2
        bulk_state = replace(self.inlet_flow_state,
                             pressure=bulk_pressure,
                             temperature=self.bulk_temperature,)
        k = bulk_state.conductivity
        pr = bulk_state.prandtl_number
        u = mass_flux / bulk_state.density
        re = bulk_state.get_reynolds(flow_speed=u,linear_dimension=diameter)
        h_a = 0.025 * k / D * re**.8 * pr**.4 * (self.bulk_temperature / wall_temperature)**.55

        ## Old stuff to see differences in inlet, outlet, and mean properties
        # d = {'in': {}, 'out': {}, 'mean': {}}
        # for name, func in zip(d.keys(), (
        #         lambda x, y={}: getattr(self.inlet_flow_state, x)(**y) if y else getattr(self.inlet_flow_state, x),
        #         lambda x, y={}: getattr(self.outlet_flow_state, x)(**y) if y else getattr(self.outlet_flow_state, x),
        #         lambda x, y={}: ((getattr(self.outlet_flow_state, x)(**y) if y else getattr(self.outlet_flow_state, x))
        #                     + (getattr(self.inlet_flow_state, x)(**y) if y else getattr(self.inlet_flow_state, x))) / 2)):
        #     d[name]['k'] = func("conductivity")
        #     d[name]['Pr'] = func("prandtl_number")
        #     d[name]['u'] = mass_flux / func("density")
        #     d[name]['Re'] = func("get_reynolds", {"flow_speed": d[name]['u'], "linear_dimension": diameter})
        #     d[name]['h_a'] = 0.025 * d[name]['k'] / D * (d[name]['Re'] ** .8) * (d[name][
        #         'Pr'] ** .4) * (self.bulk_temperature / wall_temperature)**.55
        return h_a
