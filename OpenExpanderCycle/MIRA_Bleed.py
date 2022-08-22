from OECycle import OpenExpanderCycle
from dataclasses import dataclass
from BaseEngineCycle.SplitterMerger import Splitter, Merger
from BaseEngineCycle.Cooling import CoolingChannelSection
from BaseEngineCycle.HeatTransferSection import HeatTransferSection


@dataclass
class MIRA(OpenExpanderCycle):
    min_distance_from_throat_heat_transfer_section_2: float = .4

    @property
    def turbine_mass_flow_initial_guess(self):
        return 5

    @property
    def post_fuel_pump_splitter(self):
        return Splitter(split_ratios=(1 - 0.24925373134328358208955223880597, 0.24925373134328358208955223880597),
                        input_mass_flow=self.main_fuel_flow)

    @property
    def mid_cooling_splitter(self):
        return Splitter(split_ratios=(0.82703777335984095427435387673956, 1 - 0.82703777335984095427435387673956),
                        input_mass_flow=self.cooling_flow)

    @property
    def pre_injector_merger(self):
        return Merger(input_mass_flows=(self.post_fuel_pump_splitter.output_mass_flow_2,
                                        self.mid_cooling_splitter.output_mass_flow_2))

    @property
    def cooling_flow(self):
        return self.post_fuel_pump_splitter.output_mass_flow_1

    @property
    def cooling_flow_2(self):
        return self.mid_cooling_splitter.output_mass_flow_2

    @property
    def chamber_fuel_flow(self):
        return self.pre_injector_merger.output_mass_flow

    @property
    def turbine_inlet_temperature(self):
        return self.cooling_channel_section_2.outlet_temperature

    @property
    def heat_transfer_section_2(self):
        return HeatTransferSection(**self.convective_heat_transfer_args,
                                   max_distance_section=self.nozzle.div_length,
                                   min_distance_section=self.min_distance_from_throat_heat_transfer_section_2,
                                   radiative_heat_transfer_factor=self.radiative_heat_transfer.radiative_factor)

    @property
    def cooling_channel_section(self):
        return CoolingChannelSection(propellant_name=self.fuel.name,
                                     total_heat_transfer=self.heat_transfer_section.total_heat_transfer,
                                     # total_heat_transfer=6.749e6,
                                     outlet_pressure=self.injector.inlet_pressure,
                                     inlet_temperature=self.coolant_inlet_temperature,
                                     mass_flow=self.cooling_flow)

    @property
    def cooling_channel_section_2(self):
        return CoolingChannelSection(propellant_name=self.fuel.name,
                                     # total_heat_transfer=self.heat_transfer_section_2.total_heat_transfer,
                                     total_heat_transfer=.302e6,
                                     outlet_pressure=self.injector.inlet_pressure,
                                     inlet_temperature=self.cooling_channel_section.outlet_temperature,
                                     mass_flow=self.cooling_flow_2)


@dataclass
class MIRA_Exact(MIRA):

    @property
    def cooling_channel_section(self):
        return CoolingChannelSection(propellant_name=self.fuel.name,
                                     # total_heat_transfer=self.heat_transfer_section.total_heat_transfer,
                                     total_heat_transfer=6.749e6,
                                     outlet_pressure=66.6e5,
                                     inlet_temperature=self.coolant_inlet_temperature,
                                     mass_flow=self.cooling_flow,
                                     pressure_drop=10e5)

    @property
    def cooling_channel_section_2(self):
        return CoolingChannelSection(propellant_name=self.fuel.name,
                                     # total_heat_transfer=self.heat_transfer_section_2.total_heat_transfer,
                                     total_heat_transfer=.302e6,
                                     outlet_pressure=66.6e5,
                                     inlet_temperature=self.cooling_channel_section.outlet_temperature,
                                     mass_flow=self.cooling_flow_2,
                                     pressure_drop=0)



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
        engine = EngineClass(**args.change_to_conical_nozzle(args.mira_kwargs, throat_half_angle=math.radians(25)), iterate=True)
        if pressurething == 'exact':
            engine.thrust_chamber.show_contour()
            engine.theoretical_convective_heat_transfer.show_heat_transfer()
        print(f'Data comparison for MIRA with {pressurething.upper()} pressures')
        data = [('Name', '[Unit]', 1, 1),
                ('Throat Radius', '[m]', engine.nozzle.throat_radius, .1163 / 2),
                ('Exit Radius', '[m]', (engine.exit_area / math.pi) ** .5, .65),
                ('Subs. Length', '[m]', engine.combustion_chamber.length + engine.nozzle.conv_length, .5),
                ('Nozzle Div. Length', '[m]', engine.nozzle.div_length, 1.4),
                ('CoolS1 Out. Temp.', '[K]', engine.cooling_channel_section.outlet_temperature, 490),
                ('Turb. In. Temp.', '[K]', engine.cooling_channel_section_2.outlet_temperature, 600),
                ('Heat Transfer', '[MW]', engine.cooling_channel_section.total_heat_transfer * 1e-6, 6.749),
                ('Heat Transfer2', '[MW]', engine.cooling_channel_section_2.total_heat_transfer * 1e-6, .302),
                ('Turb. Mass Flow', '[kg/s]', engine.turbine_mass_flow, .87),
                ('Cool. Mass Flow', '[kg/s]', engine.cooling_flow, 5.03),
                ('Main Fuel Flow', '[kg/s]', engine.main_fuel_flow, 6.69),
                ('Main Oxid. Flow', '[kg/s]', engine.main_oxidizer_flow, 22.72),
                ('Pumps Power Req.', '[MW]', engine.pump_power_required * 1e-6, 0),
                ('Actual Heat Transfer', '[MW]', engine.heat_transfer_section.total_heat_transfer* 1e-6, 6.749),
                ('Actual Heat Transfer2', '[MW]', engine.heat_transfer_section_2.total_heat_transfer* 1e-6, .302),
                ('Cooling Flow 2', '[kg/s]', engine.cooling_flow_2, .87),
                # ('', '[]', 1, 1),
                ]

        for name, unit, value, expected in data:
            print(f'{name:<18} \t {unit:<7} \t {value:<8.3f} \t {expected:.3f}')
        print(f'{(125*(.1163/2)**2)**.5:.4f}')
        print('\n\n\n')
