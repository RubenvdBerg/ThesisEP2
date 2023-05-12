from plots.KwakPlots.Easy_Plots import easy_plot, double_input_plot, double_output_plot
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from plots.KwakPlots.Results_Comparison_RP1 import engine_kwargs

double_input_plot(
    engine_classes=[ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle],
    output_attribute='overall_specific_impulse',
    input2_attribute='burn_time',
    input2_range=(300,300),
    input1_attribute='combustion_chamber_pressure',
    input1_range=(3e6, 20e6),
    input1_prefix='M',
    num1=9,
    num2=1,
    savefig=r'C:\Users\rvand\PycharmProjects\ThesisEP2\plots\KwakPlots\Final_Plots\V3_RP1_overall_specific_impulse',
    **engine_kwargs,
)

# double_input_plot(
#     engine_classes=[GasGeneratorCycle, OpenExpanderCycle],
#     input1_attribute='combustion_chamber_pressure',
#     input1_range=(3e6, 10e6),
#     input2_attribute='burn_time',
#     input2_range=(300, 1200),
#     output_attribute='turbine_propellant_mass',
#     num1=8,
#     num2=3,
#     **engine_kwargs
# )
# double_output_plot(
#     engine_classes=[ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle],
#     input_attribute='combustion_chamber_pressure',
#     input_range=(3e6,10e6),
#     input_prefix='M',
#     output1_attribute='final_mass',
#     output2_attribute='tanks_plus_pressurant',
#     legend_kwargs={'loc':'upper center', 'ncol':2},
#     # ylims=((595,645),(480,530)),
#     # twinx=False,
#     **engine_kwargs,
# )
#
# double_output_plot(
#     engine_classes=[ElectricPumpCycle, GasGeneratorCycle],
#     input_attribute='combustion_chamber_pressure',
#     input_range=(3e6,10e6),
#     input_prefix='M',
#     output1_attribute='feed_system_mass',
#     output2_attribute='total_thrust_chamber_mass',
#     legend_kwargs={'loc':'upper center', 'ncol':2},
#     # ylims=((595,645),(480,530)),
#     # twinx=False,
#     **engine_kwargs,
# )
#
# double_output_plot(
#     engine_classes=[GasGeneratorCycle,OpenExpanderCycle],
#     input_attribute='combustion_chamber_pressure',
#     input_range=(3e6,10e6),
#     input_prefix='M',
#     output1_attribute='turbine_propellant_mass',
#     output2_attribute='chamber_propellant_mass',
#     legend_kwargs={'loc':'upper center', 'ncol':2},
#     **engine_kwargs,
# )
