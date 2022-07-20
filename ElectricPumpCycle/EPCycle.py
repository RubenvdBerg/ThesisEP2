from dataclasses import dataclass

from BaseEngineCycle.EngineCycle import EngineCycle
from BaseEngineCycle.Pump import Pump
from ElectricPumpCycle.EPComponents import ElectricMotor, Inverter, Battery


@dataclass
class ElectricPumpCycle(EngineCycle):
    # TODO: dataclass inheritance is dumb, so all inherited classes can only have default variables if baseclass has any
    #  Solution is to update to 3.10 and use @dataclass(kw_only=True), probably not worth the hassle
    fuel_specific_heat: float = 0  # J/(kg*K)
    electric_motor_specific_power: float = 0  # W/kg
    inverter_specific_power: float = 0  # W/kg
    battery_specific_power: float = 0  # W/kg
    battery_specific_energy: float = 0  # J/kg
    electric_motor_efficiency: float = 0  # -
    inverter_efficiency: float = 0  # -
    battery_structural_factor: float = 0  # -
    coolant_allowable_temperature_change: float = 0  # K

    _fuelpumpflow: float = None

    def __post_init__(self):
        super().__post_init__()
        self._fuelpumpflow = self.fuel.main_mass_flow
        if self.iterate:
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
            print(f'Start Coolant Flow Iteration')
            print(f'm_f_cc: {self.fuel.main_mass_flow:.5f}')
            print(f'm_fp: {self.total_fuelpump_flow:.5f}, m_f_cool_act:{self.actual_battery_coolant_flow:.5f}, '
                  f'm_f_cool_req: {self.battery.coolant_flow_required:.5f}')
        while abs(self.actual_battery_coolant_flow
                  - self.battery.coolant_flow_required) > self.battery.coolant_flow_required * self.iteration_accuracy:
            self._fuelpumpflow = self.total_fuelpump_flow
            if self.verbose:
                print(f'm_fp: {self.total_fuelpump_flow:.5f}, m_f_cool_act:{self.actual_battery_coolant_flow:.5f}, '
                      f'm_f_cool_req: {self.battery.coolant_flow_required:.5f}')
        if self.verbose:
            print(f'Coolant Flow Set\n')

    def reiterate(self):
        if self.verbose:
            print('Start reiteration')
        self.update_cea()
        self.iterate_coolant_flow()

    # Rewrite of fuel_pump to accommodate for recirculation of battery cooling fuel flow
    @property
    def fuel_pump(self):
        return Pump(propellant=self.fuel,
                    mass_flow=self._fuelpumpflow,
                    pressure_increase=self.delta_p_fuel_pump,
                    efficiency=self.fuel_pump_efficiency,
                    specific_power=self.fuel_pump_specific_power)

    @property
    def total_fuelpump_flow(self):
        return self.battery.coolant_flow_required + self.fuel.main_mass_flow

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
        return Battery(specific_power=self.battery_specific_power,
                       specific_energy=self.battery_specific_energy,
                       battery_packing_factor=self.battery_structural_factor,
                       output_power=self.inverter.input_power,
                       burn_time=self.burn_time,
                       fuel_specific_heat=self.fuel_specific_heat,
                       coolant_allowable_temperature_change=self.coolant_allowable_temperature_change)

    @property
    def actual_battery_coolant_flow(self):
        return self._fuelpumpflow - self.fuel.main_mass_flow

    @property
    def feed_system_mass(self):
        return self.electric_motor.mass + self.inverter.mass + self.pumps_mass

    @property
    def mass(self):
        return self.props_mass + self.battery.mass + self.tanks_mass + self.pressurant.mass + self.feed_system_mass
