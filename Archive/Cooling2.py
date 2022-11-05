from dataclasses import dataclass, field
from functools import cached_property
from math import log, isclose
from typing import Optional
from scipy import constants, integrate
from Archive.NISTCoolant import get_heat_capacity
from time import sleep


@dataclass
class Coolant:
    propellant_name: str
    # heat_of_vaporization: float = field(init=False)
    # molar_mass: float = field(init=False)
    # boiling_temperature_1_bar: float = field(init=False)
    # specific_heat_of_vaporization: float = field(init=False)

    def get_CoolProp(self):
        get_CoolProp_properties()
        pass

    @property
    def coolant_name(self):
        p_name = self.propellant_name.upper()
        if 'RP' in p_name:
            return 'RP1'
        if 'H2' in p_name:
            return 'H2'
        if 'O2' in p_name or 'OX' in p_name:
            return 'O2'

    @property
    def heat_of_vaporization(self):  # [J/mol]
        _h_vap_dict = {'H2': 899.08248, 'O2': 3407.787, 'RP1': 41820.0}
        return _h_vap_dict[self.coolant_name]

    @property
    def molar_mass(self):  # [kg/mol]
        _mm_dict = {'H2': 2.01588E-3, 'O2': 15.999E-3, 'RP1': 170.0E-3}
        return _mm_dict[self.coolant_name]

    @property
    def boiling_temperature_1_bar(self):  # [K]
        _boil_dict = {'H2': 20.25, 'O2': 90.15, 'RP1': 489.45}
        return _boil_dict[self.coolant_name]

    @property
    def specific_heat_of_vaporization(self):
        return self.heat_of_vaporization / self.molar_mass

    def liquid_heat_capacity(self, temperature):  # [J/(mol*K)]
        return get_heat_capacity(f'L{self.coolant_name}', temp=temperature)

    def gas_heat_capacity(self, temperature):  # [J/(mol*K)]
        return get_heat_capacity(f'G{self.coolant_name}', temp=temperature)

    def specific_heat_capacity_liquid(self, temperature):  # [J/(kg*K)]
        return self.liquid_heat_capacity(temperature=temperature) / self.molar_mass

    def specific_heat_capacity_gas(self, temperature):  # [J/(kg*K)]
        return self.gas_heat_capacity(temperature=temperature) / self.molar_mass

    def boiling_temperature(self, pressure):  # [K]
        # Clausius-Clapeyron estimation
        T0 = self.boiling_temperature_1_bar
        Dh = self.heat_of_vaporization
        P0 = 1e5  # 1 bar
        Tboil = (T0 ** -1 - (constants.R / Dh) * log(pressure/P0)) ** -1
        return Tboil


@dataclass
class CoolingChannels:
    coolant: Coolant
    total_heat_transfer: float  # [W]
    outlet_pressure: float  # [Pa]
    mass_flow: float  # [kg/s]
    pressure_drop: Optional[float] = None  # [Pa]
    inlet_temperature: Optional[float] = None  # [K]
    _pressure_drop_ratio: float = field(init=False, default=.15)  # [-]
    _outlet_temp_estimate: float = field(init=False, default=200)  # [K]
    verbose: bool = True

    def __post_init__(self):
        # Set optional variables
        if self.inlet_temperature is None:
            self.inlet_temperature = self.default_inlet_temperature
        if self.pressure_drop is None:
            # Humble 1995 p.209 suggest pressure drop to be 10% - 20% of chamber/outlet pressure
            self.pressure_drop = self.outlet_pressure * self._pressure_drop_ratio

        self.iterate_cp()

    @property
    def verbose_message(self):
        return f'T_out -> [K] Estimate: {self._outlet_temp_estimate:.2f} Calculated:{self.outlet_temperature:.2f}\n' \
               f'Average_cp [J/(kg*K)] -> Gas: {self.average_coolant_gas_specific_heat_capacity:.4e} ' \
               f'Liquid: {self.average_coolant_liquid_specific_heat_capacity:.4e}'

    def iterate_cp(self):
        if self.verbose:
            print('Setting specific heat capacities:')
            print(f'T_boil: {self.boiling_temperature:.2f}\n')
            print(self.verbose_message)
        while not isclose(self.outlet_temperature, self._outlet_temp_estimate, rel_tol=0.05):
            self._outlet_temp_estimate = self.outlet_temperature
            print(self.verbose_message)
            sleep(1)
        if self.verbose:
            print('Cp set')

    @property
    def default_inlet_temperature(self):
        _temp_in_dict = {'H2': 20.25, 'O2': 90.15, 'RP1': 293.15}
        return _temp_in_dict[self.coolant.coolant_name]

    @property
    def inlet_pressure(self):
        return self.outlet_pressure + self.pressure_drop

    @property
    def inlet_specific_heat_capacity(self):
        return self.coolant.specific_heat_capacity_liquid(temperature=self.inlet_temperature)

    @cached_property
    def boiling_temperature(self):
        return self.coolant.boiling_temperature(self.inlet_pressure)

    @property
    def start_boiling_enthalpy(self):
        return self.boiling_temperature * self.average_coolant_liquid_specific_heat_capacity

    @property
    def end_boiling_enthalpy(self):
        return self.start_boiling_enthalpy + self.coolant.heat_of_vaporization

    @property
    def average_coolant_liquid_specific_heat_capacity(self):
        min_gas_temp = self.inlet_temperature
        max_gas_temp = self.boiling_temperature if self._outlet_temp_estimate > self.boiling_temperature else self._outlet_temp_estimate
        cumulative_cp, *_ = integrate.quad(func=lambda x: self.coolant.specific_heat_capacity_liquid(temperature=x),
                                           a=min_gas_temp,
                                           b=max_gas_temp)
        return cumulative_cp / (max_gas_temp - min_gas_temp)

    @property
    def average_coolant_gas_specific_heat_capacity(self):
        min_gas_temp = self.boiling_temperature
        max_gas_temp = self._outlet_temp_estimate
        cumulative_cp, *_ = integrate.quad(func=lambda x: self.coolant.specific_heat_capacity_gas(temperature=x),
                                           a=min_gas_temp,
                                           b=max_gas_temp)
        return cumulative_cp / (max_gas_temp - min_gas_temp)

    @property
    def inlet_enthalpy(self):
        if self.boiling_temperature < self.inlet_temperature:
            raise ValueError(
                'The boiling temperature of the coolant is below the inlet temperature, liquid cooling expected')
        return self.inlet_temperature * self.inlet_specific_heat_capacity

    @property
    def enthalpy_increase(self):
        return self.total_heat_transfer / self.mass_flow

    @property
    def outlet_enthalpy(self):
        return self.inlet_enthalpy + self.enthalpy_increase


    @property
    def outlet_temperature(self):
        # Return temperature dependent on enthalpy and phase (or phase transition)
        if self.outlet_enthalpy > self.end_boiling_enthalpy:
            return ((self.outlet_enthalpy - self.end_boiling_enthalpy)
                    / self.average_coolant_gas_specific_heat_capacity + self.boiling_temperature)
        elif self.outlet_enthalpy > self.start_boiling_enthalpy:
            return self.boiling_temperature
        else:
            return self.outlet_enthalpy / self.average_coolant_liquid_specific_heat_capacity

