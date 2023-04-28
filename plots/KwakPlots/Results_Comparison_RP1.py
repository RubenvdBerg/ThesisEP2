from isp_plot import plot_all_ratio_plots_data, get_engines, plot_all_ratio_from_data, open_data_dict
from Easy_Plots import easy_plot, double_input_plot
# plot_all_ratio_plots_data(short=False,
#                           is_frozen=True,
#                           names=('EP', 'GG','OE'),
#                           savefig=True,
#                           savedata=True,
#                           kwak=False,
#                           isp_single_figure=True,
#                           # burn_times=(1200,),
#                           burn_times=(300, 750, 1200),
#                           thrusts=(100e3,),
#                           # p_cc_range=(10, 11),
#                           is_ratio=True,
#                           fuel_name='RP1_NASA',
#                           mass_mixture_ratio=2.45,
#                           exit_pressure_forced=.002e6,
#                           expansion_ratio_end_cooling=2,
#                           errors=False,
#                           verbose=True,
#                           maximum_wall_temperature=600,
#                           )
# base_path = r'C:\Users\rvand\PycharmProjects\ThesisEP2\plots\KwakPlots\all_plots_data'
# # path = base_path + r'\RP1_20mbar_10eps_end_750_Final'
# path = base_path + r'\RP1_20mbar_10eps_end_750_Ft100_Only_Final'
# data_dict = open_data_dict(path)
#
# plot_all_ratio_from_data(data_dict,
#                          savefig=False,
#                          isp_single_figure=True,
#                          is_ratio=True,
#                          title=False,
#                          ylims={'dv': (0.70, 1.00),
#                                 'mr': (0.35, 0.80),
#                                 'm0': (0.99, 1.044),
#                                 'isp': (None, (1.0, 1.06))})
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle

# DO NOT CHANGE!!!!!
engine_kwargs = {
    'thrust': 100e3,
    'burn_time': 300,
    'exit_pressure_forced': 0.002e6,
    'expansion_ratio_end_cooling': 10,
    'combustion_chamber_pressure': 10e6,
    'maximum_wall_temperature': 900,
}

# easy_plot(
#     engine_classes=[ElectricPumpCycle],
#     input_prefix='',
#     input_attribute='burn_time',
#     input_range=(300, 1200),
#     output_attribute='battery.power_mass',
#     # output_attribute='battery.energy_mass',
#     **engine_kwargs
# )

def make_rp1_graphs(savefig:bool=False):

    for tb in [300, 750, 1200]:
        savefig_path = f'Final_Plots/V2_RP1_m0_{tb}s' if savefig else None
        double_input_plot(engine_classes=[ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle],
                          input1_attribute='combustion_chamber_pressure',
                          input2_attribute='burn_time',
                          input1_range=(3e6,10e6),
                          input2_range=(tb,tb),
                          output_attribute='initial_mass',
                          num1=8,
                          num2=1,
                          input1_prefix='M',
                          input2_prefix='',
                          output_prefix='',
                          savefig=savefig_path,
                          **engine_kwargs)
    for attribute in ['mass_ratio', 'change_in_velocity']:
        savefig_path = f'Final_Plots/V2_RP1_{attribute}' if savefig else None
        double_input_plot(engine_classes=[ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle],
                          input1_attribute='combustion_chamber_pressure',
                          input2_attribute='burn_time',
                          input1_range=(3e6, 10e6),
                          input2_range=(1200, 300),
                          output_attribute=attribute,
                          num1=8,
                          num2=3,
                          input1_prefix='M',
                          input2_prefix='',
                          output_prefix='',
                          savefig=savefig_path,
                          **engine_kwargs)
    double_input_plot(engine_classes=[ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle],
                      input1_attribute='combustion_chamber_pressure',
                      input2_attribute='burn_time',
                      input1_range=(3e6, 10e6),
                      input2_range=(1200, 1200),
                      output_attribute='overall_specific_impulse',
                      num1=8,
                      num2=1,
                      input1_prefix='M',
                      input2_prefix='',
                      output_prefix='',
                      savefig=f'Final_Plots/V2_RP1_overall_specific_impulse' if savefig else None,
                      **engine_kwargs)


if __name__ == '__main__':
    make_rp1_graphs(savefig=False)
