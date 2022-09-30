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

hyprob_kwargs = args.hyprob_kwargs
engine_kwargs = {key: value for key, value in hyprob_kwargs.items() if key != 'throat_area'}
engine = EngineCycle(**engine_kwargs, **args.duel_pump_kwargs)

injector = Injector(material_density=0, safety_factor=0, yield_strength=0)
nozzle = ConicalNozzle(throat_area=hyprob_kwargs['throat_area'],
                       expansion_ratio=hyprob_kwargs['expansion_ratio'],
                       conv_chamber_bend_ratio=hyprob_kwargs['convergent_chamber_bend_ratio'],
                       conv_throat_bend_ratio=hyprob_kwargs['convergent_throat_bend_ratio'],
                       conv_half_angle=hyprob_kwargs['convergent_half_angle'],
                       div_throat_half_angle=hyprob_kwargs['divergent_throat_half_angle'],
                       area_ratio_chamber_throat=hyprob_kwargs['area_ratio_chamber_throat'])
chamber = CombustionChamber(material_density=0, safety_factor=0, yield_strength=0,
                            throat_area=hyprob_kwargs['throat_area'],
                            combustion_chamber_pressure=hyprob_kwargs['combustion_chamber_pressure'],
                            convergent_volume_estimate=nozzle.conv_volume_estimate,
                            area_ratio_chamber_throat=hyprob_kwargs['area_ratio_chamber_throat'],
                            propellant_mix=engine.propellant_mix_name,
                            characteristic_length=hyprob_kwargs['chamber_characteristic_length'],
                            verbose=True)
manual_cc_flow_state = ManualFlowState(propellant_name='Combustion_Chamber_Gas',
                                       temperature=engine.combustion_temperature,
                                       pressure=hyprob_kwargs['combustion_chamber_pressure'],
                                       mass_flow=engine.chamber_mass_flow,
                                       type='combusted',
                                       _specific_heat_capacity=engine.cc_hot_gas_specific_heat_capacity,
                                       _heat_capacity_ratio=engine.cc_hot_gas_heat_capacity_ratio,
                                       _prandtl_number=engine.cc_hot_gas_prandtl_number,
                                       _dynamic_viscosity=engine.cc_hot_gas_dynamic_viscosity,
                                       )

thrustchamber = ThrustChamber(injector=injector, chamber=chamber, nozzle=nozzle,
                              heat_capacity_ratio=manual_cc_flow_state.heat_capacity_ratio)

# thrustchamber.show_contour(distance_from_injector=True)

channel_height = 2.45716897722*1e-3
channel_width = 4.53032220526*1e-3
channel_area = channel_width * channel_height
channel_equivalent_diameter = 2 * (channel_area / pi)**.5

Dh = 0.003186198562223674

heattransfer = HeatExchanger(coolant_inlet_flow_state=replace(engine.cooling_inlet_flow_state,
                                                              mass_flow=1.92,
                                                              pressure=15583600.000000002,
                                                              temperature=112.384, ),
                             coolant_channel_diameter=Dh,
                             number_of_coolant_channels=96,
                             radiative_factor=engine.radiative_heat_transfer.radiative_factor,
                             thrust_chamber=thrustchamber,
                             combustion_chamber_flow_state=manual_cc_flow_state,
                             amount_of_sections=500,
                             verbose=False,
                             wall_conductivity=365,
                             wall_thickness=.9e-3,
                             )

# print(heattransfer.data)
heattransfer.plot_all()