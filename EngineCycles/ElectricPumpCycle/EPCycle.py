from dataclasses import dataclass, replace, field
from typing import Optional
from EngineCycles.BaseEngineCycle.EngineCycle import EngineCycle
from EngineCycles.BaseEngineCycle.Pump import Pump
from EngineCycles.ElectricPumpCycle.EPComponents import ElectricMotor, Inverter, Battery, BatteryCooler
from EngineCycles.BaseEngineCycle.Merger import Merger
from EngineCycles.BaseEngineCycle.Splitter import Splitter
from EngineCycles.BaseEngineCycle.FlowState import FlowState, DefaultFlowState


@dataclass
class ElectricPumpCycle(EngineCycle):
    # TODO: dataclass inheritance is dumb, so all inherited classes can only have default variables if baseclass has any
    #  Solution is to update to 3.10 and use @dataclass(kw_only=True), probably not worth the hassle
    electric_motor_specific_power: float = 0  # W/kg
    inverter_specific_power: float = 0  # W/kg
    battery_specific_power: float = 0  # W/kg
    battery_specific_energy: float = 0  # J/kg
    electric_motor_efficiency: float = 0  # -
    inverter_efficiency: float = 0  # -
    battery_structural_factor: float = 0  # -
    battery_coolant_temperature_change: float = 0  # K
    electric_motor_heat_loss_factor: float = 0
    electric_motor_magnet_temp_limit: float = 0
    electric_motor_ox_leak_factor: float = 0
    battery_coolant_specific_heat_capacity: Optional[float] = None  # J/(kg*K)

    # Iteration attribute not required at init
    _iterative_battery_cooler_outlet_flow_state: FlowState = field(init=False, repr=False, default=DefaultFlowState())

    def __post_init__(self):
        """"Initial iterative_flow_state """
        super().__post_init__()
        if self.iterate:
            self.iterate_coolant_flow()
        self.electric_motor.calc_cooling()

    def set_initial_values(self):
        super().set_initial_values()
        self._iterative_battery_cooler_outlet_flow_state = replace(self.fuel_tank.outlet_flow_state, mass_flow=0)

    def iterate_coolant_flow(self):
        if self.verbose:
            print(f'Start Battery Coolant Flow Iteration')
            print(f'm_f_cc: {self.main_fuel_flow:.5f}')
            print(self.verbose_iteration_message)

        while self.battery_flow_error_larger_than_accuracy():
            self._iterative_battery_cooler_outlet_flow_state = self.battery_cooler.outlet_flow_state
            if self.verbose:
                print(self.verbose_iteration_message)
        if self.verbose:
            print(f'Battery Coolant Flow Set\n')

    def battery_flow_error_larger_than_accuracy(self):
        error = abs(self.actual_battery_coolant_flow - self.battery_cooler.coolant_flow_required)
        margin = self.battery_cooler.coolant_flow_required * self.iteration_accuracy
        return error > margin

    @property
    def verbose_iteration_message(self):
        return f'm_fp: {self.fuel_pump.inlet_mass_flow:.5f}, m_f_cool_act:{self.actual_battery_coolant_flow:.5f}, ' \
               f'm_f_cool_req: {self.battery_cooler.coolant_flow_required:.5f}'

    def reiterate(self):
        if self.verbose:
            print('Start reiteration')
        self.update_cea()
        self.iterate_coolant_flow()

    @property
    def pre_fuel_pump_merger(self):
        return Merger(inlet_flow_states=(self.fuel_tank.outlet_flow_state, self._iterative_battery_cooler_outlet_flow_state))

    @property
    def post_fuel_pump_splitter(self):
        return Splitter(inlet_flow_state=self.fuel_pump.outlet_flow_state,
                        required_outlet_mass_flows=(self.fuel_tank.outlet_mass_flow,),
                        outlet_flow_names=('chamber', 'battery'))

    # Rewrite of fuel_pump to accommodate for recirculation of battery cooling fuel flow (instead of actual split flow)
    @property
    def fuel_pump(self):
        return Pump(inlet_flow_state=self.pre_fuel_pump_merger.outlet_flow_state,
                    expected_outlet_pressure=self.fuel_pump_outlet_pressure,
                    efficiency=self.fuel_pump_efficiency,
                    specific_power=self.fuel_pump_specific_power,
                    propellant_density=self.fuel_density,)

    @property
    def electric_motor(self):
        return ElectricMotor(specific_power=self.electric_motor_specific_power,
                             electric_energy_efficiency=self.electric_motor_efficiency,
                             output_power=self.pumps_power_required,
                             electric_heat_loss_factor=self.electric_motor_heat_loss_factor,
                             oxidizer_pump_inlet_flow_state=self.oxidizer_pump.inlet_flow_state,
                             oxidizer_leakage_factor=self.electric_motor_ox_leak_factor,
                             magnet_temp_limit=self.electric_motor_magnet_temp_limit, )

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
                       burn_time=self.burn_time,)

    @property
    def battery_cooler(self):
        return BatteryCooler(inlet_flow_state=self.post_fuel_pump_splitter.outlet_flow_state_battery,
                             outlet_pressure_required=self.fuel_tank.outlet_pressure,
                             coolant_allowable_temperature_change=self.battery_coolant_temperature_change,
                             coolant_specific_heat_capacity=self.battery_coolant_specific_heat_capacity,
                             power_heat_loss=self.battery.power_heat_loss,)

    @property
    def cooling_inlet_flow_state(self):
        return self.post_fuel_pump_splitter.outlet_flow_state_chamber

    @property
    def actual_battery_coolant_flow(self):
        return self.post_fuel_pump_splitter.outlet_flow_state_battery.mass_flow

    @property
    def feed_system_mass(self):
        return self.electric_motor.mass + self.inverter.mass + self.pumps_mass

    @property
    def mass(self):
        return super().mass + self.battery.mass + self.inverter.mass + self.electric_motor.mass

if __name__ == '__main__':
    import arguments as args
    print(ElectricPumpCycle(thrust=100e3, burn_time=300, combustion_chamber_pressure=1e6, is_frozen=True, **args.base_arguments, **args.ep_arguments))