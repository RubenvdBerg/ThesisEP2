from dataclasses import dataclass, field
from functools import cached_property
from math import log, isclose, isnan
from typing import Optional, Callable
from scipy import constants, integrate, interpolate
import numpy as np
from BaseEngineCycle.NISTCoolant import get_heat_capacity
from time import sleep


@dataclass
class RP1Coolant:
    heat_of_vaporization: float = field(init=False, default=41820.0)  # [J/mol]
    molar_mass: float = field(init=False, default=170.0E-3)  # [kg/mol]
    boil_temp_1_bar: float = field(init=False, default=489.45)  # [K]

    @cached_property
    def specific_heat_of_vaporization(self):  # [J/kg]
        return self.heat_of_vaporization / self.molar_mass

    @cached_property
    def base_data(self) -> np.array:
        return np.genfromtxt('Data/Abdulagatov2011_RP1_Heat_Capacity_Data.csv', delimiter=',')

    @cached_property
    def cp_interpolation_func(self) -> Callable:
        base_array = self.base_data
        ts = base_array[4:, 0]
        ps = base_array[0, 2:]
        cps = base_array[4:, 2:]
        return interpolate.interp2d(ps, ts, cps, bounds_error=True)

    @cached_property
    def cp_interpolation_func_low_t(self) -> Callable:
        base_array = self.base_data
        ts = base_array[1:4, 0]
        ps = base_array[0, 1:]
        cps = base_array[1:4, 1:]
        return interpolate.interp2d(ps, ts, cps, bounds_error=True)

    def get_heat_capacity(self, temperature: float, pressure: float) -> float:
        if temperature < 373.42:
            heat_capacity = float(self.cp_interpolation_func_low_t(pressure, temperature))
        if pressure < 10E6:
            heat_capacity = float(self.cp_interpolation_charlie_taylor)
        else:
            heat_capacity = float(self.cp_interpolation_func(pressure, temperature))
        # if isnan(heat_capacity):
        #     raise ValueError('No heat capacity data available in this range')
        return heat_capacity

    def get_specific_heat_capacity(self, **kwargs) -> float:  # [J/(kg*K)]
        return self.get_heat_capacity(**kwargs) / self.molar_mass

    def get_boiling_temperature(self, pressure: float) -> float:  # [K]
        # Clausius-Clapeyron estimation
        T0 = self.boil_temp_1_bar
        Dh = self.heat_of_vaporization
        P0 = 1e5  # 1 bar
        Tboil = (T0 ** -1 - (constants.R / Dh) * log(pressure/P0)) ** -1
        return Tboil



# To prevent RP1_Coolant from being a mutable default value
def default_factory_rp1_coolant():
    return RP1Coolant()


@dataclass
class RP1CoolingChannels:
    # Similar to Cooling2 CoolingChannels, but assuming boiling point is not reached
    rp1_coolant: RP1Coolant = field(init=False, default_factory=default_factory_rp1_coolant)
    total_heat_transfer: float  # [W]
    outlet_pressure: float  # [Pa]
    mass_flow: float  # [kg/s]
    inlet_temperature: float = 293.15  # [K]
    verbose: bool = True

    # Optional attributes that can be estimated
    pressure_drop: Optional[float] = None  # [Pa]

    # Internal attributes that can be overridden if really needed
    _pressure_drop_ratio: float = field(init=False, default=.15)  # [-]
    _outlet_temp_estimate: float = field(init=False, default=200)  # [K]


    def __post_init__(self):
        # Set optional variables
        if self.pressure_drop is None:
            # Humble 1995 p.209 suggest pressure drop to be 10% - 20% of chamber/outlet pressure
            self.pressure_drop = self.outlet_pressure * self._pressure_drop_ratio

        self.iterate_cp()

    @property
    def verbose_message(self):
        return f'T_out -> [K] Estimate: {self._outlet_temp_estimate:.2f} Calculated:{self.outlet_temperature:.2f}\n' \
               f'Average_cp [J/(kg*K)]: {self.average_coolant_specific_heat_capacity:.4e}'

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
    def inlet_pressure(self):
        return self.outlet_pressure + self.pressure_drop

    def get_specific_heat_capacity(self, temperature: float) -> float:
        return self.rp1_coolant.get_specific_heat_capacity(temperature=temperature, pressure=self.inlet_pressure)

    @property
    def inlet_specific_heat_capacity(self):
        return self.get_specific_heat_capacity(temperature=self.inlet_temperature)

    @cached_property
    def boiling_temperature(self):
        return self.rp1_coolant.get_boiling_temperature(pressure=self.inlet_pressure)

    @property
    def average_coolant_specific_heat_capacity(self):
        min_temp = self.inlet_temperature
        max_temp = self.boiling_temperature if self._outlet_temp_estimate > self.boiling_temperature else self._outlet_temp_estimate
        cumulative_cp, *_ = integrate.quad(func=lambda x: self.get_specific_heat_capacity(temperature=x),
                                           a=min_temp,
                                           b=max_temp)
        return cumulative_cp / (max_temp - min_temp)

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
        outlet_temp = self.outlet_enthalpy / self.average_coolant_specific_heat_capacity
        if outlet_temp > self.boiling_temperature:
            raise ValueError(f'RP1 Assumed to stay boiling temperature, however a higher outlet temperature was found')
        return outlet_temp

if __name__ == '__main__':
    # rp1_coolant = RP1Coolant()
    print(rp1_coolant.get_heat_capacity(temperature=633., pressure=9E6))
    # coolch = RP1CoolingChannels(total_heat_transfer=1E6,outlet_pressure=10E6,mass_flow=10)
    # print(coolch.boiling_temperature)