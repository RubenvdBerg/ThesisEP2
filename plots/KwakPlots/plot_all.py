from isp_plot import plot_all_ratio_plots
from fuel_flow_plots import plot_fuel_flow_battery
from time import strftime
import os
import json

def plot_all(savefig: bool, kwak: bool, verbose:bool=True, default_vals:bool= True):
    kwargs = {'savefig': savefig,
              'default_vals': default_vals,
              'verbose': verbose,
              'kwak': kwak}

    suffix = 'Fixed' if kwargs['kwak'] else 'Normal'
    dir_name = strftime("%Y%m%d-%H%M%S") + '_' + suffix
    if kwargs['savefig']:
        os.mkdir(dir_name)
    plot_all_ratio_plots(short=False, dir_name=dir_name, **kwargs)
    plot_fuel_flow_battery(savefig=kwargs['savefig'], dir_name=dir_name, kwak=kwargs['kwak'])
    with open(dir_name + '/' + 'Notes.txt', 'w') as f:
        f.write(json.dumps(kwargs))


if __name__ == '__main__':
    plot_all(savefig=True, kwak=True)
    # plot_all(savefig=True, kwak=False)