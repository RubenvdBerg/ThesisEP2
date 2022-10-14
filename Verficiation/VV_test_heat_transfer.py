from dataclasses import replace
from typing import Optional

import arguments as args
from BaseEngineCycle.CombustionChamber import CombustionChamber
from BaseEngineCycle.EngineCycle import EngineCycle
from BaseEngineCycle.FlowState import ManualFlowState
from BaseEngineCycle.HeatExchanger import HeatTransferPlots
from BaseEngineCycle.HeatExchangerOMECA import OMECAHeatExchanger
from BaseEngineCycle.Injector import Injector
from BaseEngineCycle.Nozzle import ConicalNozzle
from BaseEngineCycle.ThrustChamber import ThrustChamberSection


def test_heat_transfer(engine_kwargs: dict,
                       make_plots: bool = False,
                       heat_class: OMECAHeatExchanger = OMECAHeatExchanger,
                       throat_area: Optional[float] = None,
                       coolant_mass_flow: Optional[float] = None,
                       coolant_inlet_temp: Optional[float] = None,
                       coolant_inlet_pressure: Optional[float] = None,
                       **heat_exchanger_kwargs):
    engine = EngineCycle(**engine_kwargs, **args.duel_pump_kwargs)

    if heat_class is None:
        heat_class = OMECAHeatExchanger
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

    thrustchambersection = ThrustChamberSection(injector=injector, chamber=chamber, nozzle=nozzle,
                                                heat_capacity_ratio=manual_cc_flow_state.heat_capacity_ratio,
                                                min_distance=engine.distance_from_throat_start_cooling,
                                                max_distance=engine.distance_from_throat_end_cooling,
                                                min_distance_expansion_ratio=engine.expansion_ratio_start_cooling,
                                                max_distance_expansion_ratio=engine.expansion_ratio_end_cooling)

    if issubclass(heat_class, OMECAHeatExchanger):
        heat_exchanger_kwargs['characteristic_velocity'] = engine.characteristic_velocity

    heattransfer = heat_class(coolant_inlet_flow_state=replace(engine.cooling_inlet_flow_state,
                                                               mass_flow=coolant_mass_flow,
                                                               pressure=coolant_inlet_pressure,
                                                               temperature=coolant_inlet_temp, ),
                              radiative_factor=engine.radiative_heat_transfer.radiative_factor,
                              thrust_chamber_section=thrustchambersection,
                              combustion_chamber_flow_state=manual_cc_flow_state,
                              _save_data=True if make_plots else heat_exchanger_kwargs.pop('_save_data', True),
                              **heat_exchanger_kwargs)
    if make_plots:
        plots = HeatTransferPlots(heattransfer.data)
        return heattransfer, plots
    return heattransfer