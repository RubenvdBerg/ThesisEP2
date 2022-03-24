from base_engine_cycle import Structure, EngineCycle
from scipy.constants import g


class GasGenerator(Structure):
    def __init__(self, oxidizer, fuel, mass_mixture_ratio, gg_pressure, gg_gas_constant, pump_power_required,
                 temperature_limit, stay_time, heat_capacity_ratio=1.16, pressure_ratio=27, turbine_efficiency=.52,
                 specific_heat=2024.7, **kwargs):
        self.oxidizer = oxidizer  # Propellant Instance
        self.fuel = fuel  # Propellant Instance
        self.mmr = mass_mixture_ratio  # -
        self.p = gg_pressure  # Pa
        self.r_gg = gg_gas_constant  # J/(kg*K)
        self.ppr = pump_power_required  # W
        self.tl = temperature_limit  # K
        self.ts = stay_time  # s
        self.y_gg = heat_capacity_ratio  # -
        self.p_ratio = pressure_ratio  # -
        self.eta = turbine_efficiency  # -
        self.cp = specific_heat  # J/(kg*K)
        super().__init__(**kwargs)

    @property
    def mass(self):
        return self.sf * 3 / 2 * self.md / self.sy * self.ts * self.p / self.gas_density * self.turbine_mass_flow

    @property
    def gas_density(self):
        return self.p / (self.r_gg * self.tl)

    @property
    def turbine_mass_flow(self):
        return self.ppr / (self.eta * self.cp * self.tl * (1 - self.p_ratio ** ((1 - self.y_gg) / self.y_gg)))


class GasGeneratorCycle(EngineCycle):
    def __init__(self, gg_gas_specific_heat, heat_ratio_gg_gas, mass_mixture_ratio_gg, turbine_pressure_ratio,
                 gas_constant_gg_gas, turbine_inlet_temperature, gas_generator_stay_time, turbopump_specific_power,
                 turbine_efficiency, gg_structural_factor, gg_material_density, gg_yield_strength,
                 gg_thrust_contribution, verbose=False, **kwargs):
        self.cp_gg = gg_gas_specific_heat  # J/(kg*K)
        self.y_gg = heat_ratio_gg_gas  # -
        self.mmr_gg = mass_mixture_ratio_gg  # -
        self.k_tu = turbine_pressure_ratio  # -
        self.r_gg = gas_constant_gg_gas  # J/(kg*K)
        self.temp_tu = turbine_inlet_temperature  # K
        self.t_s = gas_generator_stay_time  # s
        self.d_tp = turbopump_specific_power  # W/kg
        self.eta_tu = turbine_efficiency  # -
        self.k_gg = gg_structural_factor  # -
        self.rho_gg = gg_material_density  # kg/m3
        self.s_y_gg = gg_yield_strength  # Pa
        super().__init__(fuel_pump_specific_power=self.d_tp, oxidizer_pump_specific_power=self.d_tp, kwak_fix_cycle_type='gg', **kwargs)
        self.k_thrust_gg = 0.002137040978335770 if self.kwak_fix else gg_thrust_contribution  # -
        self.gg_mass_flow = 0
        self.verbose = verbose
        self.iterate_mass_flow()

    def iterate_mass_flow(self):
        while self.gg_mass_flow*1.0001 < self.gas_generator.turbine_mass_flow:
            if self.verbose:
                print(f'GG:{self.gg_mass_flow:.6f}kg/s')
                print(f'TU1:{self.gas_generator.turbine_mass_flow:.6f}kg/s')
            self.gg_mass_flow = self.gas_generator.turbine_mass_flow
        if self.verbose:
            print(f'Mass Flow Set')

    def reiterate(self):
        if self.verbose:
            print('Start reiteration')
        self.cstar, self.cf = self.set_cea()
        self.iterate_mass_flow()

    @property
    def chamber_mass_flow(self):
        return (1 - self.k_thrust_gg) * self.mass_flow

    @property
    def base_fuel_flow(self):
        return 1 / (self.mmr + 1) * self.mass_flow

    @property
    def base_oxidizer_flow(self):
        return self.mmr / (self.mmr + 1) * self.mass_flow

    @property
    def total_mass_flow(self):
        return self.chamber_mass_flow + self.gg_mass_flow

    # Rewrite Isp to use total mass flow
    @property
    def simple_specific_impulse(self):
        return self.thrust / self.total_mass_flow / g

    @property
    def fuel_flow(self):  # Override EngineCycle flows
        return 1 / (self.mmr + 1) * self.chamber_mass_flow + 1 / (self.mmr_gg + 1) * self.gg_mass_flow

    @property
    def oxidizer_flow(self):  # Override EngineCycle flows
        return self.mmr / (self.mmr + 1) * self.chamber_mass_flow + self.mmr_gg / (self.mmr_gg + 1) * self.gg_mass_flow

    @property
    def gas_generator(self):
        return GasGenerator(oxidizer=self.oxidizer, fuel=self.fuel, mass_mixture_ratio=self.mmr_gg,
                            gg_pressure=self.p_cc, gg_gas_constant=self.r_gg,
                            pump_power_required=self.pump_power_required, temperature_limit=self.temp_tu,
                            stay_time=self.t_s, heat_capacity_ratio=self.y_gg, pressure_ratio=self.k_tu,
                            turbine_efficiency=self.eta_tu, specific_heat=self.cp_gg, safety_factor=self.k_gg,
                            material_density=self.rho_gg, sigma_yield=self.s_y_gg)


    @property
    def pumps_mass(self):
        return self.fuel_pump.mass + self.oxidizer_pump.mass

    @property
    def gg_propellant_mass(self):
        if self.kwak_fix:
            return self.gg_mass_flow * self.t_b
        return self.gg_mass_flow * self.t_b * self.k_p

    @property
    def cc_propellant_mass(self):
        return self.chamber_mass_flow * self.t_b * self.k_p

    @property
    def feed_system_mass(self):
        return self.gas_generator.mass + self.pumps_mass

    @property
    def mass(self):
        return (self.cc_propellant_mass
                + self.gg_propellant_mass
                + self.feed_system_mass
                + self.tanks_mass
                + self.pressurant.mass)
