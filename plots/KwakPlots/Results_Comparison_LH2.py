from isp_plot import plot_all_ratio_plots_data, get_engines, plot_all_ratio_from_data, open_data_dict

# a = get_engines('ep',
#                 thrust=100e3,
#                 burn_time=1200,
#                 # burn_time=300,
#                 is_frozen=True,
#                 verbose=True,
#                 kwak=False,
#                 errors=True,
#                 fuel_name='RP1_NASA',
#                 mass_mixture_ratio=5.5,
#                 p_cc_range=(10, 11),
#                 expansion_ratio=100,
#                 battery_coolant_temperature_change=40,)
# print(a)
plot_all_ratio_plots_data(short=False,
                          is_frozen=True,
                          names=('EP', 'GG', 'OE'),
                          savefig=True,
                          savedata=True,
                          kwak=False,
                          burn_times=(300,),
                          thrusts=(100e3,),
                          p_cc_range=(10, 11),
                          isp_single_figure=True,
                          is_ratio=True,
                          fuel_name='LH2_NASA',
                          mass_mixture_ratio=5.5,
                          expansion_ratio_end_cooling=10,
                          exit_pressure_forced=0.002e6,
                          errors=True,
                          verbose=True,
                          )
# path = r'C:\Users\rvand\PycharmProjects\ThesisEP2\plots\KwakPlots\all_plots_data\20230411-172752_Normal'
# data_dict = open_data_dict(path)
#
# plot_all_ratio_from_data(data_dict,
#                          savefig=True,
#                          isp_single_figure=True,
#                          is_ratio=True,
#                          title=False,
#                          ylims={'dv': (0.6, 0.9),
#                                 'mr': (0.30, 0.75),
#                                 'm0': (0.99, 1.17),
#                                 'isp': (None, (1.0, 1.02))})
