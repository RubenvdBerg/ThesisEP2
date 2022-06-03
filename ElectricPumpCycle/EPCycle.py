from dataclasses import dataclass

from BaseEngineCycle.EngineCycle import EngineCycle
from BaseEngineCycle.TurboPump import Pump
from ElectricPumpCycle.EPComponents import ElectricMotor, Inverter, Battery, KwakBattery


@dataclass
class ElectricPumpCycle(EngineCycle):
    # TODO: dataclass inheritance is dumb, so all inherited classes can only have default variables if baseclass has any. Solution is to update to 3.10 and use kw_only, probably not worth the hassle
    fuel_specific_heat: float = 1  # J/(kg*K)
    electric_motor_specific_power: float = 1  # W/kg
    inverter_specific_power: float = 1  # W/kg
    battery_specific_power: float = 1  # W/kg
    battery_specific_energy: float = 1  # J/kg
    electric_motor_efficiency: float = 1  # -
    inverter_efficiency: float = 1  # -
    battery_structural_factor: float = 1  # -
    coolant_allowable_temperature_change: float = 1  # K
    verbose: float = 1
    _kwak_fix_cycle_type: str = 'ep'

    def __post_init__(self):
        super().__post_init__()
        self.fuelpumpflow = self.fuel.mass_flow
        self.iterate_coolant_flow()

    # def iterate_coolant_flow2(self):
    #     format_string = '.7f'
    #     # print(f'm_dot_cc_f:{self.fuel.mass_flow}, m_dot_f1:{self.fuelpumpflow}')
    #     # print(f'm_dot_tot1:{self.total_fuelpump_flow}, m_cool:{self.battery.coolant_flow_required}')
    #     self.fuelpumpflow = self.total_fuelpump_flow
    #     # print(f'm_dot_tot2:{self.total_fuelpump_flow}, m_dot_f2:{self.fuelpumpflow}')
    #     while self.fuel_pump.mass_flow * 1.01 <= self.total_fuelpump_flow:
    #         if self.verbose:
    #             print(f'm_fp: {self.fuel_pump.mass_flow:{format_string}}, '
    #                   f'm_c: {self.battery.coolant_flow_required:{format_string}}, '
    #                   f'm_cc,f {self.fuel.mass_flow:{format_string}}, '
    #                   f'm_tfp: {self.total_fuelpump_flow:{format_string}}')
    #         self.fuelpumpflow = self.total_fuelpump_flow
    #         if self.verbose:
    #             print(f'2m_fp: {self.fuel_pump.mass_flow:{format_string}}, '
    #                   f'm_c: {self.battery.coolant_flow_required:{format_string}}, '
    #                   f'm_cc,f {self.fuel.mass_flow:{format_string}}, '
    #                   f'm_tfp: {self.total_fuelpump_flow:{format_string}}')

    def iterate_coolant_flow(self):
        if self.verbose:
            print(f'm_f_cc: {self.fuel.mass_flow}')
            print(
                f'm_fp: {self.total_fuelpump_flow}, m_f_cool_act:{self.actual_battery_coolant_flow}, m_f_cool_req: {self.battery.coolant_flow_required}')
        while self.actual_battery_coolant_flow * 1.001 < self.battery.coolant_flow_required:
            self.fuelpumpflow = self.total_fuelpump_flow
            if self.verbose:
                print(
                    f'm_fp: {self.total_fuelpump_flow}, m_f_cool_act:{self.actual_battery_coolant_flow}, m_f_cool_req: {self.battery.coolant_flow_required}')

    def update_cea(self):
        if self.verbose:
            print('Start reiteration')
        self.cstar, self.cf = self.set_cea()
        self.iterate_coolant_flow()

    # Rewrite of fuel_pump to accommodate for recirculation of battery cooling fuel flow
    @property
    def fuel_pump(self):
        return Pump(propellant=self.fuel,
                    mass_flow=self.fuelpumpflow,
                    pressure_increase=self.delta_p_fuel_pump,
                    efficiency=self.fuel_pump_efficiency,
                    specific_power=self.fuel_pump_specific_power)

    @property
    def total_fuelpump_flow(self):
        return self.battery.coolant_flow_required + self.fuel.mass_flow

    @property
    def electric_motor(self):
        return ElectricMotor(specific_power=self.electric_motor_specific_power,
                             electric_energy_efficiency=self.electric_motor_efficiency,
                             output_power=self.pump_power_required)

    @property
    def inverter(self):
        return Inverter(specific_power=self.inverter_specific_power,
                        electric_energy_efficiency=self.inverter_efficiency,
                        output_power=self.electric_motor.input_power)

    @property
    def battery(self):
        if self.kwak_fix:
            return KwakBattery(specific_power=self.battery_specific_power,
                               specific_energy=self.battery_specific_energy,
                               battery_packing_factor=self.battery_structural_factor,
                               output_power=self.electric_motor.output_power,
                               burn_time=self.burn_time,
                               fuel_specific_heat=self.fuel_specific_heat,
                               coolant_allowable_temperature_change=self.coolant_allowable_temperature_change,
                               inverter_efficiency=self.inverter.electric_energy_efficiency,
                               electric_motor_efficiency=self.electric_motor.electric_energy_efficiency)
        else:
            return Battery(specific_power=self.battery_specific_power,
                           specific_energy=self.battery_specific_energy,
                           battery_packing_factor=self.battery_structural_factor,
                           output_power=self.inverter.input_power,
                           burn_time=self.burn_time,
                           fuel_specific_heat=self.fuel_specific_heat,
                           coolant_allowable_temperature_change=self.coolant_allowable_temperature_change)

    @property
    def actual_battery_coolant_flow(self):
        return self.fuelpumpflow - self.fuel.mass_flow

    @property
    def feed_system_mass(self):
        return self.electric_motor.mass + self.inverter.mass + self.pumps_mass

    @property
    def mass(self):
        return self.props_mass + self.battery.mass + self.tanks_mass + self.pressurant.mass + self.feed_system_mass