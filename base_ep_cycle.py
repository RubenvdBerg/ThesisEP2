from base_engine_cycle import EngineCycle, Pump
from math import log


class ElectricComponent:
    def __init__(self, output_power, electric_energy_efficiency, specific_power):
        self.eta = electric_energy_efficiency  # -
        self.d_p = specific_power  # W/kg
        self.output_power = output_power  # W

    @property
    def input_power(self):
        return self.output_power / self.eta

    @property
    def mass(self):
        return self.output_power / self.d_p


class ElectricMotor(ElectricComponent):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


class Inverter(ElectricComponent):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


class Battery(ElectricComponent):
    def __init__(self, specific_energy, battery_packing_factor, burn_time, fuel_specific_heat,
                 coolant_allowable_temperature_change, **kwargs):
        self.d_e = specific_energy  # J/kg
        self.k_bp = battery_packing_factor  # -
        self.t_b = burn_time  # s
        self.cp_f = fuel_specific_heat  # J/(kg*K)
        self.delta_temp = coolant_allowable_temperature_change  # K
        super().__init__(electric_energy_efficiency=self.eta_e, **kwargs)

    @property
    def eta_e(self):
        return 0.093 * log(self.t_b) + 0.3301

    @property
    def total_energy(self):
        return self.input_power * self.t_b

    @property
    def heat_loss(self):
        return self.total_energy * (1 - self.eta_e)

    @property
    def coolant_flow_required(self):
        return (1 - self.eta_e) * self.total_energy / (self.cp_f * self.delta_temp * self.t_b)

    # Overwrite Super _mass_
    @property
    def mass(self):
        return self.k_bp * self.output_power * max(1 / self.d_p, self.t_b / (self.d_e * self.eta_e))


class KwakBattery(Battery):
    def __init__(self, inverter_efficiency, electric_motor_efficiency, **kwargs):
        self.eta_em = electric_motor_efficiency
        self.eta_inv = inverter_efficiency
        super().__init__(**kwargs)

    @property
    def total_energy(self):
        return self.output_power * self.t_b

    # Overwrite Super _mass_
    @property
    def mass(self):
        return self.k_bp * self.output_power / (self.eta_em * self.eta_inv) * max(1 / self.d_p,
                                                                                  self.t_b / (self.d_e * self.eta_e))


class ElectricPumpCycle(EngineCycle):
    def __init__(self, fuel_specific_heat, electric_motor_specific_power, inverter_specific_power,
                 battery_specific_power, battery_specific_energy, electric_motor_efficiency, inverter_efficiency,
                 battery_structural_factor, coolant_allowable_temperature_change, verbose=False, **kwargs):
        self.cp_f = fuel_specific_heat  # J/(kg*K)
        self.d_em = electric_motor_specific_power  # W/kg
        self.d_inv = inverter_specific_power  # W/kg
        self.d_bat_p = battery_specific_power  # W/kg
        self.d_bat_e = battery_specific_energy  # J/kg
        self.eta_em = electric_motor_efficiency  # -
        self.eta_inv = inverter_efficiency  # -
        self.k_bat = battery_structural_factor  # -
        self.delta_temp = coolant_allowable_temperature_change  # K
        self.verbose = verbose
        super().__init__(**kwargs, kwak_fix_cycle_type='ep')
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
        while self.actual_battery_coolant_flow * 1.0001 < self.battery.coolant_flow_required:
            self.fuelpumpflow = self.total_fuelpump_flow
            if self.verbose:
                print(
                    f'm_fp: {self.total_fuelpump_flow}, m_f_cool_act:{self.actual_battery_coolant_flow}, m_f_cool_req: {self.battery.coolant_flow_required}')

    def reiterate(self):
        if self.verbose:
            print('Start reiteration')
        self.cstar, self.cf = self.set_cea()
        self.iterate_coolant_flow()

    # Rewrite of fuel_pump to accommodate for recirculation of battery cooling fuel flow
    @property
    def fuel_pump(self):
        return Pump(propellant_density=self.fuel.density, mass_flow=self.fuelpumpflow,
                    pressure_change=self.delta_p_fuel_pump, efficiency=self.eta_fp, specific_power=self.d_fp)

    @property
    def total_fuelpump_flow(self):
        return self.battery.coolant_flow_required + self.fuel.mass_flow

    @property
    def electric_motor(self):
        return ElectricMotor(specific_power=self.d_em, electric_energy_efficiency=self.eta_em,
                             output_power=self.pump_power_required)

    @property
    def inverter(self):
        return Inverter(specific_power=self.d_inv, electric_energy_efficiency=self.eta_inv,
                        output_power=self.electric_motor.input_power)

    @property
    def battery(self):
        if self.kwak_fix:
            return KwakBattery(specific_power=self.d_bat_p, specific_energy=self.d_bat_e,
                               battery_packing_factor=self.k_bat,
                               output_power=self.electric_motor.output_power, burn_time=self.t_b,
                               fuel_specific_heat=self.cp_f,
                               coolant_allowable_temperature_change=self.delta_temp,
                               inverter_efficiency=self.inverter.eta, electric_motor_efficiency=self.electric_motor.eta)
        else:
            return Battery(specific_power=self.d_bat_p, specific_energy=self.d_bat_e, battery_packing_factor=self.k_bat,
                           output_power=self.inverter.input_power, burn_time=self.t_b, fuel_specific_heat=self.cp_f,
                           coolant_allowable_temperature_change=self.delta_temp)

    @property
    def actual_battery_coolant_flow(self):
        return self.fuelpumpflow - self.fuel.mass_flow

    @property
    def feed_system_mass(self):
        return self.electric_motor.mass + self.inverter.mass + self.pumps_mass

    @property
    def mass(self):
        return self.props_mass + self.battery.mass + self.tanks_mass + self.pressurant.mass + self.feed_system_mass
