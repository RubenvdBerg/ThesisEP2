from EngineCycles.Abstract.EngineCycle import EngineCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from dataclasses import dataclass
from scipy import constants
from KwakFix.KwakFixComponents import KwakBattery, KwakPump, KwakTank, KwakPropellant
from EngineComponents.Abstract.FlowState import ManualFlowState


@dataclass
class KwakEngineCycle(EngineCycle):
    pressurant_heat_capacity_ratio: float = 1.667,
    pressurant_molar_mass: float = 0.00399733779,
    manual_oxidizer_density: float = 0
    manual_fuel_density: float = 0

    def calc_pump_outlet_pressures(self):
        pass

    @property
    def pressurant_initial_state(self):
        super_state = super().pressurant_initial_state
        return ManualFlowState(propellant_name=super_state.propellant_name,
                               temperature=super_state.temperature,
                               pressure=super_state.pressure,
                               mass_flow=super_state.mass_flow,
                               type=super_state.type,
                               _molar_mass=self.pressurant_molar_mass,
                               _heat_capacity_ratio=self.pressurant_heat_capacity_ratio, )

    @property
    def oxidizer(self):
        return KwakPropellant(main_flow_state=self.oxidizer_main_flow_state,
                              burn_time=self.burn_time,
                              margin_factor=self.propellant_margin_factor,
                              manual_propellant_density=self.manual_oxidizer_density)

    @property
    def fuel(self):
        return KwakPropellant(main_flow_state=self.fuel_main_flow_state,
                              burn_time=self.burn_time,
                              margin_factor=self.propellant_margin_factor,
                              manual_propellant_density=self.manual_fuel_density)

    @property
    def oxidizer_tank(self):
        return KwakTank(inlet_flow_state=self.oxidizer_main_flow_state,
                        propellant_volume=self.oxidizer.volume,
                        max_acceleration=self.max_acceleration,
                        ullage_factor=self.ullage_volume_factor,
                        pressurant_tank_volume=self.pressurant_tank.volume,
                        structure_material=self.oxidizer_tank_material,
                        safety_factor=self.tanks_structural_factor,
                        manual_propellant_density=self.manual_oxidizer_density, )

    @property
    def fuel_tank(self):
        return KwakTank(inlet_flow_state=self.fuel_main_flow_state,
                        propellant_volume=self.fuel.volume,
                        max_acceleration=self.max_acceleration,
                        ullage_factor=self.ullage_volume_factor,
                        pressurant_tank_volume=None,
                        structure_material=self.fuel_tank_material,
                        safety_factor=self.tanks_structural_factor,
                        manual_propellant_density=self.manual_fuel_density, )

    @property
    def oxidizer_pump(self):
        return KwakPump(inlet_flow_state=self.oxidizer_tank.outlet_flow_state,
                        expected_outlet_pressure=self.oxidizer_pump_outlet_pressure,
                        efficiency=self.oxidizer_pump_efficiency,
                        specific_power=self.oxidizer_pump_specific_power,
                        manual_propellant_density=self.manual_oxidizer_density, )

    @property
    def fuel_pump(self):
        return KwakPump(inlet_flow_state=self.fuel_tank.outlet_flow_state,
                        expected_outlet_pressure=self.fuel_pump_outlet_pressure,
                        efficiency=self.fuel_pump_efficiency,
                        specific_power=self.fuel_pump_specific_power,
                        manual_propellant_density=self.manual_fuel_density, )


@dataclass
class KwakFixGasGeneratorCycle(GasGeneratorCycle, KwakEngineCycle):
    _iteration_done: bool = False

    def set_initial_values(self):
        super().set_initial_values()
        r = self.gg_base_flow_state.specific_gas_constant
        self.gg_base_flow_state._density = self.gg_pressure / (r * self.turbine_maximum_temperature)

    def iterate_flow(self):
        m_tu = self.turbine.mass_flow_required
        m_gg_o = m_tu * self.gg_mass_mixture_ratio / (1 + self.gg_mass_mixture_ratio)
        m_gg_f = m_tu * 1 / (1 + self.gg_mass_mixture_ratio)
        self._iterative_turbine_mass_flow += m_gg_f + m_gg_o
        m_tu = self.turbine.mass_flow_required
        m_gg_o = m_tu * self.gg_mass_mixture_ratio / (1 + self.gg_mass_mixture_ratio)
        m_gg_f = m_tu * 1 / (1 + self.gg_mass_mixture_ratio)
        self.mo = self.main_oxidizer_flow - m_gg_o
        self.mf = self.main_fuel_flow - m_gg_f
        self._iterative_turbine_mass_flow = m_tu
        self._iteration_done = True

    @property
    def turbine_mass_flow_initial_guess(self):
        return 0

    @property
    def overall_specific_impulse(self):
        """Calculate specific impulse without accounting for turbine exhaust thrust contribution."""
        return self.chamber_mass_flow * self.ideal_thrust_coefficient * self.characteristic_velocity / self.total_mass_flow / constants.g

    @property
    def chamber_mass_flow(self):
        if self._iteration_done:
            return self.mf + self.mo
        else:
            return self.base_mass_flow

    @property
    def chamber_fuel_flow(self):
        if self._iteration_done:
            return self.mf
        else:
            return super().chamber_fuel_flow

    @property
    def chamber_oxidizer_flow(self):
        if self._iteration_done:
            return self.mo
        else:
            return super().chamber_oxidizer_flow

    @property
    def gg_propellant_mass(self):
        return self.gg_mass_flow * self.burn_time


@dataclass
class KwakFixElectricPumpCycle(ElectricPumpCycle, KwakEngineCycle):
    battery_fix_on: bool = True

    @property
    def battery(self):
        if self.battery_fix_on:
            return KwakBattery(specific_power=self.battery_specific_power,
                               specific_energy=self.battery_specific_energy,
                               battery_packing_factor=self.battery_structural_factor,
                               output_power=self.pumps_power_required,
                               burn_time=self.burn_time,
                               inverter_efficiency=self.inverter.electric_energy_efficiency,
                               electric_motor_efficiency=self.electric_motor.electric_energy_efficiency)
        else:
            return super().battery

    @property
    def fuel_pump(self):
        return KwakPump(inlet_flow_state=self.pre_fuel_pump_merger.outlet_flow_state,
                        expected_outlet_pressure=self.fuel_pump_outlet_pressure,
                        efficiency=self.fuel_pump_efficiency,
                        specific_power=self.fuel_pump_specific_power,
                        manual_propellant_density=self.manual_fuel_density)
