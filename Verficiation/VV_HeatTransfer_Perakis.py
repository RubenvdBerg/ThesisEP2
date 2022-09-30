from BaseEngineCycle.EngineCycle import EngineCycle
from BaseEngineCycle.Nozzle import ConicalNozzle
from BaseEngineCycle.Injector import Injector
from BaseEngineCycle.ThrustChamber import ThrustChamber
from BaseEngineCycle.CombustionChamber import CombustionChamber
from BaseEngineCycle.FlowState import ManualFlowState
from BaseEngineCycle.HeatTransferSection2 import HeatExchanger
from dataclasses import replace
import arguments as args
from math import radians
import pandas as pd

main_kwargs = args.change_to_conical_nozzle(args.tcd1_kwargs, throat_half_angle=radians(25))
engine = EngineCycle(**main_kwargs, **args.duel_pump_kwargs)
injector = Injector(material_density=0, safety_factor=0, yield_strength=0)
nozzle = ConicalNozzle(throat_area=engine.throat_area,
                       expansion_ratio=main_kwargs['expansion_ratio'],
                       conv_chamber_bend_ratio=main_kwargs['convergent_chamber_bend_ratio'],
                       conv_throat_bend_ratio=main_kwargs['convergent_throat_bend_ratio'],
                       conv_half_angle=main_kwargs['convergent_half_angle'],
                       div_throat_half_angle=main_kwargs['divergent_throat_half_angle'],
                       area_ratio_chamber_throat=main_kwargs['area_ratio_chamber_throat'])
chamber = CombustionChamber(material_density=0, safety_factor=0, yield_strength=0,
                            throat_area=engine.throat_area,
                            combustion_chamber_pressure=engine.combustion_chamber_pressure,
                            convergent_volume_estimate=nozzle.conv_volume_estimate,
                            area_ratio_chamber_throat=main_kwargs['area_ratio_chamber_throat'],
                            propellant_mix=engine.propellant_mix_name,
                            characteristic_length=main_kwargs['chamber_characteristic_length'],
                            verbose=True)
manual_cc_flow_state = ManualFlowState(propellant_name='Combustion_Chamber_Gas',
                                       temperature=engine.combustion_temperature,
                                       pressure=engine.combustion_chamber_pressure,
                                       mass_flow=engine.chamber_mass_flow,
                                       type='combusted',
                                       _specific_heat_capacity=engine.cc_hot_gas_specific_heat_capacity,
                                       _heat_capacity_ratio=engine.cc_hot_gas_heat_capacity_ratio,
                                       _prandtl_number=engine.cc_hot_gas_prandtl_number,
                                       _dynamic_viscosity=engine.cc_hot_gas_dynamic_viscosity,
                                       )

thrustchamber = ThrustChamber(injector=injector, chamber=chamber, nozzle=nozzle,
                              heat_capacity_ratio=manual_cc_flow_state.heat_capacity_ratio)

heattransfer = HeatExchanger(coolant_inlet_flow_state=replace(engine.cooling_inlet_flow_state,
                                                              mass_flow=10,
                                                              temperature=30,
                                                              pressure=150e5),
                             coolant_channel_diameter=.0025,
                             number_of_coolant_channels=138,
                             radiative_factor=engine.radiative_heat_transfer.radiative_factor,
                             thrust_chamber=thrustchamber,
                             combustion_chamber_flow_state=manual_cc_flow_state,
                             amount_of_sections=200,
                             verbose=False,
                             counter_flow=False,
                             chamber_wall_conductivity=350,
                             chamber_wall_thickness=1e-3,
                             )

heattransfer.plot_all()
