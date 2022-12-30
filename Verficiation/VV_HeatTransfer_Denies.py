from EngineArguments import arguments as args
from Verficiation.VV_test_heat_transfer import test_heat_transfer
from Archive.HeatExchangerDenies import DeniesHeatExchanger

original = False
heatt, plots = test_heat_transfer(engine_kwargs=args.denies_kwargs,
                                  throat_area=0.001433726,
                                  heat_class=DeniesHeatExchanger,
                                  number_of_coolant_channels=64 if original else 72,
                                  chamber_wall_thickness=4.2e-3 if original else 1e-3,
                                  chamber_wall_conductivity=295 if original else 365,
                                  coolant_mass_flow=0.763461538462 if original else .76,
                                  coolant_inlet_temp=110,
                                  coolant_inlet_pressure=60e5,
                                  coolant_heat_transfer_coefficient_mode='DittusBoelter',
                                  is_counter_flow=True,
                                  verbose=False,
                                  make_plots=True,
                                  hot_gas_convective_heat_transfer_coefficient_mode='Bartz',
                                  amount_of_sections=284,
                                  iteration_accuracy=1e-6)
plots.plot_all()
