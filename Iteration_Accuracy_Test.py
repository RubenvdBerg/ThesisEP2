from EngineCycles.ElectricPumpCycle import ElectricPumpCycle
from EngineCycles.GasGeneratorCycle import GasGeneratorCycle
from EngineCycles.OpenExpanderCycle import OpenExpanderCycle
from EngineArguments.default_arguments import get_default_kwargs
for EngineClass in [ElectricPumpCycle, GasGeneratorCycle, OpenExpanderCycle]:
    print(EngineClass.__name__)
    data = []
    for i in [2,3,4,14]:
        ec = EngineClass(thrust=100e3,
                         burn_time=900,
                         combustion_chamber_pressure=10e6,
                         expansion_ratio=100,
                         iteration_accuracy=10**(i*-1),
                         verbose=False,
                         **get_default_kwargs(EngineClass))
        data.append([i, ec.initial_mass, ec.final_mass, ec.ideal_delta_v])

    def get_change(old_val, new_val):
        return (new_val - old_val)/old_val
    def get_change_row(row):
        data1 = [f'{get_change(data[-1][i],row[i]):.3e}' for i in range(1,len(row))]
        print(row[0], *data1)

    for row in data:
        get_change_row(row)