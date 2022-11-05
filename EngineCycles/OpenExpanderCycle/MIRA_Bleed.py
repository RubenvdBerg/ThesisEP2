from OECycle import CoolantBleedCycle
from dataclasses import dataclass, replace
from EngineCycles.BaseEngineCycle.Merger import Merger
from EngineCycles.BaseEngineCycle.Splitter import Splitter
from EngineCycles.BaseEngineCycle.Cooling import CoolingChannelSection
from EngineCycles.BaseEngineCycle.HeatTransferSection import HeatTransferSection


@dataclass
class MIRA(CoolantBleedCycle):
    """See 'Leonardi et al. 2017 - Basic Analysis of a LOX/Methane Expander Bleed Engine' for this rocket's
    configuration
    """
    min_distance_from_throat_heat_transfer_section_2: float = .4

    @property
    def default_turbine_flow_check_state(self):
        return self.secondary_cooling_channel_section.outlet_flow_state

    @property
    def post_fuel_pump_splitter(self):
        return Splitter(inlet_flow_state=self.fuel_pump.outlet_flow_state,
                        mass_flow_fractions=(
                            1 - 0.24925373134328358208955223880597, 0.24925373134328358208955223880597),
                        outlet_flow_names=('cooling', 'chamber'))

    @property
    def post_cooling_splitter(self):
        return Splitter(inlet_flow_state=self.cooling_channel_section.outlet_flow_state,
                        required_outlet_mass_flows=(
                            self.chamber_fuel_flow - self.post_fuel_pump_splitter.outlet_flow_state_chamber.mass_flow,),
                        outlet_flow_names=('chamber', 'secondary_cooling'))

    @property
    def pre_injector_merger(self):
        return Merger(inlet_flow_states=(self.post_fuel_pump_splitter.outlet_flow_state_chamber,
                                         self.post_cooling_splitter.outlet_flow_state_chamber))

    @property
    def injector_inlet_flow_state_fuel(self):
        return self.pre_injector_merger.outlet_flow_state

    @property
    def turbine_inlet_temperature(self):
        return self.secondary_cooling_channel_section.outlet_temperature

    @property
    def secondary_heat_transfer_section(self):
        return HeatTransferSection(**self.convective_heat_transfer_args,
                                   max_distance_section=self.nozzle.div_length,
                                   min_distance_section=self.min_distance_from_throat_heat_transfer_section_2,
                                   radiative_heat_transfer_factor=self.radiative_heat_transfer.radiative_factor)

    @property
    def cooling_inlet_flow_state(self):
        return self.post_fuel_pump_splitter.outlet_flow_state_cooling

    # @property
    # def cooling_channel_section(self):
    #     return replace(super().cooling_channel_section,
    # total_heat_transfer=6.749e6,)

    @property
    def secondary_cooling_inlet_flow_state(self):
        return self.post_cooling_splitter.outlet_flow_state_secondary_cooling

    @property
    def secondary_cooling_channel_section(self):
        return CoolingChannelSection(inlet_flow_state=self.secondary_cooling_inlet_flow_state,
                                     # _total_heat_transfer=self.secondary_heat_transfer_section.total_heat_transfer,
                                     _total_heat_transfer=.302e6,
                                     combustion_chamber_pressure=self.combustion_chamber_pressure,
                                     pressure_drop=self.cooling_section_pressure_drop,
                                     verbose=self.verbose, )


@dataclass
class MIRA_Exact(MIRA):
    @property
    def cooling_channel_section(self):
        return replace(super().cooling_channel_section,
                       _total_heat_transfer=6.749e6,
                       pressure_drop=self.cooling_inlet_flow_state.pressure - 66.6e5)

    @property
    def secondary_cooling_channel_section(self):
        return replace(super().secondary_cooling_channel_section,
                       _total_heat_transfer=.302e6,
                       pressure_drop=self.secondary_cooling_inlet_flow_state.pressure - 66.6e5, )


if __name__ == '__main__':
    import arguments as args
    import math

    # from BaseEngineCycle.Turbine import Turbine
    # from BaseEngineCycle.Pump import Pump
    # from BaseEngineCycle.Propellant import Propellant
    # ox = Propellant(name='oxygen', type= 'oxidizer', main_mass_flow=1, burn_time=1, density=1142.5502246693018, margin_factor=1.1)
    # fu = Propellant(name='methane', type='fuel', main_mass_flow=1, burn_time=1, density=417.61703160433984, margin_factor=1.1)
    # test_ox_pump = Pump(propellant=ox, mass_flow=22.72, pressure_increase=6.24e6 - .3e6, efficiency=.7, specific_power=1, inlet_pressure=.3)
    # test_fu_pump = Pump(propellant=fu, mass_flow=6.7, pressure_increase=7.68e6 - .3e6, efficiency=.7, specific_power=1, inlet_pressure=.3)
    # power = test_fu_pump.power_required + test_ox_pump.power_required
    # test_turbine = Turbine(pump_power_required=power,efficiency=.6,specific_heat_capacity=3067.9,heat_capacity_ratio=1.22, pressure_ratio=14.478260869565217391304347826087, inlet_temperature=600)
    # print(test_turbine.mass_flow_required)
    for EngineClass, pressurething in zip((MIRA_Exact, MIRA), ('exact', 'estimated')):
        engine = EngineClass(**args.change_to_conical_nozzle(args.mira_kwargs, throat_half_angle=math.radians(25)),
                             iterate=True, verbose=True)
        print(engine.cooling_channel_section.throat_wall_temperature(engine.maximum_wall_temperature))
        if pressurething == 'exact':
            engine.thrust_chamber.show_contour()
            engine.theoretical_convective_heat_transfer.show_heat_transfer()
        print(f'Data comparison for MIRA with {pressurething.upper()} pressures')
        data = [('Name', '[Unit]', 1, 1),
                ('Throat Radius', '[m]', engine.nozzle.throat_radius, .1163 / 2),
                ('Exit Radius', '[m]', (engine.exit_area / math.pi) ** .5, .65),
                ('Subs. Length', '[m]', engine.combustion_chamber.length + engine.nozzle.conv_length, .5),
                ('Nozzle Div. Length', '[m]', engine.nozzle.div_length, 1.4),
                ('Heat Transfer', '[MW]', engine.cooling_channel_section.total_heat_transfer * 1e-6, 6.749),
                ('Heat Transfer2', '[MW]', engine.secondary_cooling_channel_section.total_heat_transfer * 1e-6, .302),
                ('Cool Out. Temp.', '[K]', engine.cooling_channel_section.outlet_temperature, 490),
                ('Cool Out. Temp.2', '[K]', engine.secondary_cooling_channel_section.outlet_temperature, 600),
                ('Cool. Mass Flow', '[kg/s]', engine.cooling_channel_section.inlet_mass_flow, 5.03),
                ('Cool. Mass Flow2', '[kg/s]', engine.secondary_cooling_channel_section.inlet_mass_flow, .87),
                ('Main Fuel Flow', '[kg/s]', engine.main_fuel_flow, 6.69),
                ('Main Oxid. Flow', '[kg/s]', engine.main_oxidizer_flow, 22.72),
                ('Turb. Temp.', '[K]', engine.turbine_inlet_temperature, 600),
                ('Turb. FLow', '[kg/s]', engine.turbine_mass_flow, .87),
                ('CC. Fuel Flow', '[kg/s]', engine.chamber_fuel_flow, 5.82),
                ('CC. Oxid. Flow', '[kg/s]', engine.chamber_oxidizer_flow, 22.72),
                ('CC Isp. Vac.', '[s]', engine.chamber_vacuum_specific_impulse, 352.5),
                ('2E Isp. Vac.', '[s]', engine.turbine_exhaust.vacuum_specific_impulse, 141.4),

                # ('', '[]', 1, 1),
                ('Act. Heat Trans.', '[MW]', engine.heat_transfer_section.total_heat_transfer * 1e-6, 6.749),
                ('Act. Heat Trans.2', '[MW]', engine.secondary_heat_transfer_section.total_heat_transfer * 1e-6, .302),
                ('Pumps Power Req.', '[MW]', engine.pump_power_required * 1e-6, 0),
                # ('', '[]', 1, 1),
                ]

        for name, unit, value, expected in data:
            print(f'{name:<18} \t {unit:<7} \t {value:<8.3f} \t {expected:.3f}')
        # print(f'{(125*(.1163/2)**2)**.5:.4f}')
