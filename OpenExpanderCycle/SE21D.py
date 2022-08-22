from OECycle import OpenExpanderCycle
from dataclasses import dataclass
from BaseEngineCycle.Pump import Pump
from BaseEngineCycle.SplitterMerger import Splitter, Merger
from BaseEngineCycle.Cooling import CoolingChannelSection


@dataclass
class SE21D(OpenExpanderCycle):
    fuel_pump2_efficiency: float = 0.75

    @property
    def post_fuel_pump_splitter(self):
        return Splitter(split_ratios=(1 - 0.17500560436393138123552205984393, 0.17500560436393138123552205984393),
                        input_mass_flow=self.main_fuel_flow)
    #
    # @property
    # def post_cooling_splitter(self):
    #     return Splitter(split_ratios=(0.37940710015859460778333536659754, 1 - 0.37940710015859460778333536659754),
    #                     input_mass_flow=self.cooling_flow)

    @property
    def pre_injection_merger(self):
        return Merger(input_mass_flows=(self.post_cooling_splitter.output_mass_flow_1,
                                        self.post_fuel_pump_splitter.output_mass_flow_1))

    @property
    def fuel_pump2(self):
        return Pump(propellant=self.fuel,
                    pressure_increase=self.delta_p_fuel_pump2,
                    efficiency=self.fuel_pump2_efficiency,
                    specific_power=self.fuel_pump_specific_power,
                    mass_flow=self.post_fuel_pump_splitter.output_mass_flow_2,
                    inlet_pressure=self.fuel_pump.outlet_pressure)

    @property
    def cooling_flow(self):
        return self.fuel_pump2.mass_flow

    @property
    def chamber_fuel_flow(self):
        return self.pre_injection_merger.output_mass_flow

    @property
    def delta_p_fuel_pump2(self):
        return .2 * self.combustion_chamber_pressure

    @property
    def delta_p_oxidizer_pump(self):
        return self.injector.inlet_pressure - self.oxidizer_initial_pressure

    @property
    def delta_p_fuel_pump(self):
        return self.injector.inlet_pressure - self.fuel_initial_pressure

    @property
    def pump_power_required(self):
        return self.fuel_pump.power_required + self.fuel_pump2.power_required + self.oxidizer_pump.power_required

@dataclass
class SE21D_PressureExact(SE21D):

    @property
    def delta_p_fuel_pump(self):
        return 8.749e6 - self.fuel_initial_pressure

    @property
    def delta_p_fuel_pump2(self):
        return 12.109e6 - self.fuel_pump.outlet_pressure

    @property
    def delta_p_oxidizer_pump(self):
        return 8.749e6 - self.oxidizer_initial_pressure

    @property
    def cooling_channel_section(self):
        return CoolingChannelSection(propellant_name=self.fuel.name,
                                     total_heat_transfer=self.heat_transfer_section.total_heat_transfer,
                                     outlet_pressure=8.749e6,
                                     inlet_temperature=self.coolant_inlet_temperature,
                                     mass_flow=self.cooling_flow,
                                     pressure_drop=12.109e6-8.749e6)

    def turbine2(self):
        return Turbine(pump_power_required=self.pump_power_required,
                       inlet_temperature=self.turbine_inlet_temperature,
                       efficiency=self.turbine_efficiency,
                       specific_heat_capacity=self.turbine_gas_specific_heat_capacity,
                       heat_capacity_ratio=self.turbine_gas_heat_capacity_ratio,
                       pressure_ratio=self.turbine_pressure_ratio)


if __name__ == '__main__':

    import arguments as args
    from BaseEngineCycle.Turbine import Turbine
    # test_turbine = Turbine(pump_power_required=20.768e6, efficiency=.45, specific_heat_capacity=14515.8, heat_capacity_ratio=1.398, pressure_ratio=27.7033333, inlet_temperature=506)
    # print(test_turbine.mass_flow_required)
    for EngineClass, pressurething in zip((SE21D_PressureExact, SE21D),('exact', 'estimated')):
        engine = EngineClass(**args.change_to_conical_nozzle(args.se_21d_kwargs), iterate=True)
        engine.thrust_chamber.show_contour()
        engine.heat_transfer_section.show_heat_flux()
        print(f'Data comparison for SE21 with {pressurething.upper()} pressures')
        data = [('Name', '[Unit]', 1, 1),
                ('Chamber Diameter', '[m]', engine.combustion_chamber.radius * 2, 0.985),
                ('Chamber Volume', '[m3]', engine.combustion_chamber.volume_incl_convergent, 1.029),
                ('Subs. Length', '[m]', engine.combustion_chamber.length + engine.nozzle.conv_length, 1.582),
                ('Throat Radius', '[m]', engine.nozzle.throat_radius, 0.286),
                ('Nozzle Length', '[m]', engine.nozzle.total_length, 2.548),
                # ('TC Length', '[m]', engine.thrust_chamber.length, 0),
                ('Heat Transfer', '[MW]', engine.heat_transfer_section.total_heat_transfer*1e-6, 111.037),
                ('Turb. In. Temp.', '[K]', engine.turbine_inlet_temperature, 506),
                ('Turb. Mass Flow', '[kg/s]', engine.turbine_mass_flow, 10.174),
                ('Cool. Mass Flow', '[kg/s]', engine.cooling_flow, 16.394),
                ('Main Fuel Flow', '[kg/s]', engine.main_fuel_flow, 93.677),
                ('Main Oxid. Flow', '[kg/s]', engine.main_oxidizer_flow, 456.323),
                ('Fuel Pump1 Power Req.', '[MW]', engine.fuel_pump.power_required*1e-6, 15.443),
                ('Fuel Pump2 Power Req.', '[MW]', engine.fuel_pump2.power_required * 1e-6, 1.01),
                ('Oxid. Pump Power Req.', '[MW]', engine.oxidizer_pump.power_required * 1e-6, 4.315),
                ('Total Pump Power Req.', '[MW]', engine.pump_power_required * 1e-6, 20.768),
                # ('', '[]', 1, 1),
                ]

        for name, unit, value, expected in data:
            print(f'{name:<18} \t {unit:<7} \t {value:<8.3f} \t {expected:.3f}')
        print('\n\n\n')