from OECycle import CoolantBleedCycle
from dataclasses import dataclass
from BaseEngineCycle.Pump import Pump
from BaseEngineCycle.SplitterMerger2 import Splitter, Merger
from BaseEngineCycle.Cooling import CoolingChannelSection


@dataclass
class SE21D(CoolantBleedCycle):
    """See 'Sippel et al. 2003 - Studies on Expander Bleed Cycle Engines for Launchers' for this rocket's configuration.
    """
    fuel_pump2_efficiency: float = 0.75
    
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
    def default_turbine_flow_check_state(self):
        return self.post_cooling_splitter.outlet_flow_state_turbine

    @property
    def post_fuel_pump_splitter(self):
        return Splitter(inlet_flow_state=self.fuel_pump.outlet_flow_state,
                        # mass_flow_fractions=(
                        # 1 - 0.17500560436393138123552205984393, 0.17500560436393138123552205984393),
                        mass_flow_fractions=(
                            1 - 0.2, 0.2),
                        outlet_flow_names=('chamber', 'secondary_pump'))

    # @property
    # def post_cooling_splitter(self):
    #     return Splitter(split_ratios=(0.37940710015859460778333536659754, 1 - 0.37940710015859460778333536659754),
    #                     input_mass_flow=self.cooling_flow)

    @property
    def post_cooling_splitter(self):
        """Splits the flow into: 1. the additional required flow to the main combustion chamber and 2. the rest flow,
        which goes to the turbine
        """
        return Splitter(inlet_flow_state=self.cooling_channel_section.outlet_flow_state,
                        required_outlet_mass_flows=(self.chamber_fuel_flow
                                                    - self.post_fuel_pump_splitter.outlet_flow_state_chamber.mass_flow,),
                        outlet_flow_names=('chamber', 'turbine'))

    @property
    def pre_injection_merger(self):
        return Merger(inlet_flow_states=(self.post_cooling_splitter.outlet_flow_state_chamber,
                                         self.post_fuel_pump_splitter.outlet_flow_state_chamber))

    @property
    def get_injector_inlet_flow_state_fuel(self):
        return self.pre_injection_merger.outlet_flow_state

    @property
    def secondary_fuel_pump(self):
        return Pump(inlet_flow_state=self.post_fuel_pump_splitter.outlet_flow_state_secondary_pump,
                    pressure_increase=self.delta_p_fuel_pump2,
                    efficiency=self.fuel_pump2_efficiency,
                    specific_power=self.fuel_pump_specific_power,
                    _temperature_change=32.375-29.44)

    @property
    def cooling_channel_section(self):
        ccs = super().cooling_channel_section
        ccs.inlet_flow_state = self.secondary_fuel_pump.outlet_flow_state
        return ccs

    @property
    def cooling_inlet_flow_state(self):
        return self.secondary_fuel_pump.outlet_flow_state

    # @property
    # def chamber_fuel_flow(self):
    #     return self.pre_injection_merger.output_mass_flow

    @property
    def delta_p_fuel_pump2(self):
        return .2 * self.combustion_chamber_pressure

    # @property
    # def delta_p_oxidizer_pump(self):
    #     return self.injector.inlet_pressure - self.oxidizer_initial_pressure
    #
    # @property
    # def delta_p_fuel_pump(self):
    #     return self.injector.inlet_pressure - self.fuel_initial_pressure

    @property
    def pump_power_required(self):
        return self.fuel_pump.power_required + self.secondary_fuel_pump.power_required + self.oxidizer_pump.power_required


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
        """Copies cooling channel, while adjusting the pressure drop"""
        ccs = super().cooling_channel_section
        ccs.pressure_drop = 12.109e6 - 8.749e6
        return ccs


if __name__ == '__main__':

    import arguments as args
    from BaseOpenCycle.Turbine import Turbine

    # test_turbine = Turbine(pump_power_required=20.768e6, efficiency=.45, specific_heat_capacity=14515.8, heat_capacity_ratio=1.398, pressure_ratio=27.7033333, inlet_temperature=506)
    # print(test_turbine.mass_flow_required)
    new_args = args.se_21d_kwargs | {'turbine_gas_specific_heat_capacity': None,
                                     'turbine_gas_heat_capacity_ratio': None, }
    for EngineClass, pressurething in zip((SE21D_PressureExact, SE21D), ('exact', 'estimated')):
        # if EngineClass == SE21D:
            # engine = EngineClass(**args.change_to_conical_nozzle(args.se_21d_kwargs), iterate=True)
            engine = EngineClass(**args.change_to_conical_nozzle(new_args), iterate=True)
            # engine.thrust_chamber.show_contour()
            engine.cooling_channel_section.throat_wall_temperature(engine.maximum_wall_temperature)
            # engine.heat_transfer_section.show_heat_flux()
            print(f'Data comparison for SE21 with {pressurething.upper()} pressures')
            data = [('Name', '[Unit]', 1, 1),
                    ('Chamber Diameter', '[m]', engine.combustion_chamber.radius * 2, 0.985),
                    ('Chamber Volume', '[m3]', engine.combustion_chamber.volume_incl_convergent, 1.029),
                    ('Subs. Length', '[m]', engine.combustion_chamber.length + engine.nozzle.conv_length, 1.582),
                    ('Throat Radius', '[m]', engine.nozzle.throat_radius, 0.286),
                    ('Nozzle Length', '[m]', engine.nozzle.total_length, 2.548),
                    ('DeltaP Oxid. Pump', '[MPa]', engine.oxidizer_pump.pressure_increase*1e-6, 8.749-.5),
                    ('DeltaP Fuel Pump', '[MPa]', engine.fuel_pump.pressure_increase*1e-6, 8.749-.3),
                    ('DeltaP Fuel Pump2', '[MPa]', engine.secondary_fuel_pump.pressure_increase*1e-6, 12.109-8.749),
                    ('Fuel Pump1 Power Req.', '[MW]', engine.fuel_pump.power_required * 1e-6, 15.443),
                    ('Fuel Pump2 Power Req.', '[MW]', engine.secondary_fuel_pump.power_required * 1e-6, 1.01),
                    ('Oxid. Pump Power Req.', '[MW]', engine.oxidizer_pump.power_required * 1e-6, 4.315),
                    ('Total Pump Power Req.', '[MW]', engine.pump_power_required * 1e-6, 20.768),
                    # ('TC Length', '[m]', engine.thrust_chamber.length, 0),
                    ('Heat Transfer', '[MW]', engine.heat_transfer_section.total_heat_transfer * 1e-6, 111.037),
                    ('Turb. In. Temp.', '[K]', engine.turbine_inlet_temperature, 506),
                    ('Turb. Mass Flow', '[kg/s]', engine.turbine_mass_flow, 10.174),
                    ('Cool. Mass Flow', '[kg/s]', engine.cooling_channel_section.inlet_mass_flow, 16.394),
                    ('Main Fuel Flow', '[kg/s]', engine.main_fuel_flow, 93.677),
                    ('Main Oxid. Flow', '[kg/s]', engine.main_oxidizer_flow, 456.323),

                                        # ('', '[]', 1, 1),
                    ]

            for name, unit, value, expected in data:
                print(f'{name:<18} \t {unit:<7} \t {value:<8.3f} \t {expected:.3f}')
            print('\n\n\n')
