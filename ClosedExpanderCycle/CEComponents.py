from dataclasses import dataclass
from BaseEngineCycle.Cooling2 import Coolant


@dataclass
class EXturbine:
    pass


@dataclass
class EXCoolingChannels:
    coolant: Coolant
    max_outlet_temperature: float  # [K]
    total_heat_transfer: float  # [W]
    outlet_pressure: float  # [Pa]
    inlet_temperature: float  # [K]
    mass_flow: float  # [kg/s]
    pressure_drop: Optional[float] = None  # [Pa]
    _pressure_drop_ratio: float = .15  # [-]

    @property
    @cache
    def inlet_pressure(self):
        if self.pressure_drop is None:
            # Humble 1995 p.209 suggest pressure drop to be 10% - 20% of chamber/outlet pressure
            return self.outlet_pressure * (1 + self._pressure_drop_ratio)
        else:
            return self.outlet_pressure - self.pressure_drop

    @property
    @cache
    def inlet_enthalpy(self):
        if self.coolant.boiling_temperature(self.inlet_pressure) < self.inlet_temperature:
            raise ValueError(
                'The boiling temperature of the coolant is below the inlet temperature, liquid cooling expected')
        return self.inlet_temperature * self.coolant.specific_heat_capacity_liquid

    @property
    @cache
    def enthalpy_increase(self):
        return self.total_heat_transfer / self.mass_flow

    @property
    @cache
    def outlet_enthalpy(self):
        return self.inlet_enthalpy + self.enthalpy_increase

    @property
    @cache
    def outlet_temperature(self):
        cp_g = self.coolant.specific_heat_capacity_gas
        cp_l = self.coolant.specific_heat_capacity_liquid
        temp_boil = self.coolant.boiling_temperature(self.inlet_pressure)
        h_vap_start = self.coolant.start_boiling_enthalpy(self.inlet_pressure)
        h_vap_end = self.coolant.end_boiling_enthalpy(self.inlet_pressure)

        if self.outlet_enthalpy > h_vap_end:
            return (self.outlet_enthalpy - h_vap_end) / cp_g + temp_boil
        elif self.outlet_enthalpy > h_vap_start:
            return temp_boil
        else:
            return self.outlet_enthalpy / cp_l

    @property
    @cache
    def coolant_outlet_temperature(self):
        return self.inlet_temperature + self.temperature_increase
