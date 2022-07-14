from dataclasses import dataclass, field
from typing import Optional
from functools import cached_property

@dataclass
class Coolant:
    propellant_name: str

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
        _mm_dict = {'H2': 2.01588E-3, 'O2': 15.999E-3, 'RP1': 173.0E-3}
        return _mm_dict[self.coolant_name]

    @property
    def boiling_temperature_1_bar(self):  # [K]
        _boil_dict = {'H2': 20.25, 'O2': 90.15, 'RP1': 489.45}
        return _boil_dict[self.coolant_name]

    @property
    def liquid_heat_capacity(self):  # [J/(mol*K)]
        # Values for H2 and O2 taken from NIST at 20.25 K, 90.15 K respectively and 100 bar
        # RP1 value taken from Abdulagatov2011 at 334.15 K and 100 bar
        _cp_l_dict = {'H2': 15.56, 'O2': 53.02, 'RP1': 2.151E3*self.molar_mass}
        return _cp_l_dict[self.coolant_name]

    @property
    def gas_heat_capacity(self):  # [J/(mol*K)]
        # Values for H2 and O2 Taken from NIST at 400 K
        _cp_g_dict = {'H2': 29.18, 'O2': 30.10, 'RP1': 346.}
        return _cp_g_dict[self.coolant_name]

    @property
    def default_inlet_temperature(self):
        _temp_in_dict = {'H2': 20.25, 'O2': 90.15, 'RP1': 293.15}
        return _temp_in_dict[self.coolant.coolant_name]

    @property
    def specific_heat_of_vaporization(self):  # [J/kg]
        return self.heat_of_vaporization / self.molar_mass

    @property
    def liquid_specific_heat_capacity(self):  # [J/(kg*K)]
        return self.liquid_heat_capacity / self.molar_mass

    @property
    def gas_specific_heat_capacity(self):  # [J/(kg*K)]
        return self.gas_heat_capacity / self.molar_mass

    def boiling_temperature(self, pressure):  # [K]
        # Clausius-Clapeyron estimation
        temp_0 = self.boiling_temperature_1_bar
        delta_h = self.heat_of_vaporization
        p_0 = 1e5  # 1 bar
        temp_boil = (temp_0 ** -1 - (constants.R / delta_h) * log(pressure/p_0)) ** -1
        return temp_boil


@dataclass
class Simple_CoolingChannels:
    coolant: Coolant
    total_heat_transfer: float  # [W]
    outlet_pressure: float  # [Pa]
    mass_flow: float  # [kg/s]
    pressure_drop: Optional[float] = None  # [Pa]
    inlet_temperature: Optional[float] = None  # [K]
    _pressure_drop_ratio: float = field(init=False, default=.15)  # [-]
    verbose: bool = True

    def __post_init__(self):
        # Set optional variables
        if self.inlet_temperature is None:
            self.inlet_temperature = self.coolant.default_inlet_temperature
        if self.pressure_drop is None:
            # Humble 1995 p.209 suggest pressure drop to be 10% - 20% of chamber/outlet pressure
            self.pressure_drop = self.outlet_pressure * self._pressure_drop_ratio

    @property
    def inlet_pressure(self):  # [Pa]
        return self.outlet_pressure + self.pressure_drop

    @cached_property
    def boiling_temperature(self):  # [K]
        return self.coolant.boiling_temperature(self.inlet_pressure)

    @property
    def start_boiling_enthalpy(self):  # [J/kg]
        return self.boiling_temperature * self.coolant.specific_heat_of_vaporization

    @property
    def end_boiling_enthalpy(self):  # [J/kg]
        return self.start_boiling_enthalpy + self.coolant.heat_of_vaporization

    @property
    def inlet_enthalpy(self):  # [J/kg]
        if self.boiling_temperature < self.inlet_temperature:
            raise ValueError(
                'The boiling temperature of the coolant is below the inlet temperature, liquid expected at inlet')
        return self.inlet_temperature * self.coolant.liquid_heat_capacity

    @property
    def enthalpy_increase(self):  # [J/kg]
        return self.total_heat_transfer / self.mass_flow

    @property
    def outlet_enthalpy(self):  # [J/kg]
        return self.inlet_enthalpy + self.enthalpy_increase

    @property
    def outlet_temperature(self):  # [K]
        # Return temperature dependent on enthalpy and phase (or phase transition)
        if self.outlet_enthalpy > self.end_boiling_enthalpy:
            return ((self.outlet_enthalpy - self.end_boiling_enthalpy)
                    / self.coolant.gas_specific_heat_capacity + self.boiling_temperature)
        elif self.outlet_enthalpy > self.start_boiling_enthalpy:
            return self.boiling_temperature
        else:
            return self.outlet_enthalpy / self.coolant.liquid_specific_heat_capacity
