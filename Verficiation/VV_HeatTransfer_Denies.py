from BaseEngineCycle.EngineCycle import EngineCycle
from BaseEngineCycle.Nozzle import ConicalNozzle
from BaseEngineCycle.Injector import Injector
from BaseEngineCycle.ThrustChamber import ThrustChamber
from BaseEngineCycle.CombustionChamber import CombustionChamber
from BaseEngineCycle.FlowState import ManualFlowState
from BaseEngineCycle.HeatTransferSection2 import HeatExchanger
from dataclasses import replace
import arguments as args
from math import radians, pi
import pandas as pd
from typing import Optional


def test_heat_transfer(engine_kwargs: dict,
                       number_of_coolant_channels: float,
                       chamber_wall_thickness: float,
                       chamber_wall_conductivity: float,
                       heat_class: Optional[HeatExchanger] = None,
                       throat_area: Optional[float] = None,
                       coolant_mass_flow: Optional[float] = None,
                       coolant_inlet_temp: Optional[float] = None,
                       coolant_inlet_pressure: Optional[float] = None,
                       **heat_exchanger_kwargs):
    engine = EngineCycle(**engine_kwargs, **args.duel_pump_kwargs)

    if heat_class is None:
        heat_class = HeatExchanger
    if throat_area is None:
        throat_area = engine.throat_area
    if coolant_inlet_pressure is None:
        coolant_inlet_pressure = engine.cooling_inlet_flow_state.pressure
    if coolant_inlet_temp is None:
        coolant_inlet_temp = engine.cooling_inlet_flow_state.temperature
    if coolant_mass_flow is None:
        coolant_mass_flow = engine.cooling_inlet_flow_state.mass_flow

    injector = Injector(material_density=0, safety_factor=0, yield_strength=0)
    nozzle = ConicalNozzle(throat_area=throat_area,
                           expansion_ratio=engine_kwargs['expansion_ratio'],
                           conv_chamber_bend_ratio=engine_kwargs['convergent_chamber_bend_ratio'],
                           conv_throat_bend_ratio=engine_kwargs['convergent_throat_bend_ratio'],
                           conv_half_angle=engine_kwargs['convergent_half_angle'],
                           div_throat_half_angle=engine_kwargs['divergent_throat_half_angle'],
                           area_ratio_chamber_throat=engine_kwargs['area_ratio_chamber_throat'])
    chamber = CombustionChamber(material_density=0, safety_factor=0, yield_strength=0,
                                throat_area=throat_area,
                                combustion_chamber_pressure=engine_kwargs['combustion_chamber_pressure'],
                                convergent_volume_estimate=nozzle.conv_volume_estimate,
                                area_ratio_chamber_throat=engine_kwargs['area_ratio_chamber_throat'],
                                propellant_mix=engine.propellant_mix_name,
                                characteristic_length=engine_kwargs['chamber_characteristic_length'],
                                verbose=True)
    manual_cc_flow_state = ManualFlowState(propellant_name='Combustion_Chamber_Gas',
                                           temperature=engine.combustion_temperature,
                                           pressure=engine_kwargs['combustion_chamber_pressure'],
                                           mass_flow=engine.chamber_mass_flow,
                                           type='combusted',
                                           _specific_heat_capacity=engine.cc_hot_gas_specific_heat_capacity,
                                           _heat_capacity_ratio=engine.cc_hot_gas_heat_capacity_ratio,
                                           _prandtl_number=engine.cc_hot_gas_prandtl_number,
                                           _dynamic_viscosity=engine.cc_hot_gas_dynamic_viscosity,
                                           )

    thrustchamber = ThrustChamber(injector=injector, chamber=chamber, nozzle=nozzle,
                                  heat_capacity_ratio=manual_cc_flow_state.heat_capacity_ratio)

    heattransfer = heat_class(coolant_inlet_flow_state=replace(engine.cooling_inlet_flow_state,
                                                               mass_flow=coolant_mass_flow,
                                                               pressure=coolant_inlet_pressure,
                                                               temperature=coolant_inlet_temp, ),
                              number_of_coolant_channels=number_of_coolant_channels,
                              radiative_factor=engine.radiative_heat_transfer.radiative_factor,
                              thrust_chamber=thrustchamber,
                              combustion_chamber_flow_state=manual_cc_flow_state,
                              chamber_wall_conductivity=chamber_wall_conductivity,
                              chamber_wall_thickness=chamber_wall_thickness,
                              **heat_exchanger_kwargs)
    thrustchamber.show_contour(distance_from_injector=True)
    heattransfer.plot_all()
    heattransfer.plot_geometry()


if __name__ == '__main__':
    from BaseEngineCycle.HeatTransferSection2 import DetailedHeatExchanger
    original = False
    test_heat_transfer(engine_kwargs=args.denies_kwargs,
                       throat_area=0.001433726,
                       heat_class=DetailedHeatExchanger,
                       number_of_coolant_channels=64 if original else 72,
                       chamber_wall_thickness=4.2e-3 if original else 1e-3,
                       chamber_wall_conductivity=295 if original else 365,
                       coolant_mass_flow=0.763461538462 if original else .76,
                       coolant_inlet_temp=110,
                       coolant_inlet_pressure=60e5,
                       coolant_heat_transfer_coefficient_mode='DittusBoelter',
                       counter_flow=True,
                       verbose=False,
                       hot_gas_convective_heat_transfer_coefficient_mode='Bartz2',
                       amount_of_sections=286,
                       iteration_accuracy=1e-6)
