from EngineFunctions.BaseFunctions import multi_legend
import matplotlib.pyplot as plt
import json
import numpy as np

with open('data/own_data.json', 'r') as file:
    own_data = json.load(file)


with open('data/omeca_data.json', 'r') as file:
    omeca_data = json.load(file)


xplotVals = omeca_data['Plot Vals [cm]']['x']
yplotVals = omeca_data['Plot Vals [cm]']['y']
hgvals = omeca_data['Heat-Transfer Coefficient [W/(K*m2]']['Hot Gas']
hcvals = omeca_data['Heat-Transfer Coefficient [W/(K*m2]']['Coolant']
Tvals = omeca_data['Coolant State']['T']
T0vals = omeca_data['Coolant State']['T0']
pvals = omeca_data['Coolant State']['p']
p0vals = omeca_data['Coolant State']['p0']
Twvals = omeca_data['Temperature [K]']['Hot SideWall']
TwChannelvals = omeca_data['Temperature [K]']['Cold SideWall']
qvals = omeca_data['Heat Flux [W/m2]']['Total']
channelHeightvals = omeca_data['Channel Geometry [mm]']['Height']
channelWidthvals = omeca_data['Channel Geometry [mm]']['Width']
channelDhvals = omeca_data['Channel Geometry [mm]']['Diameter']

own_x_vals= own_data['Plot Vals [cm]']['x']
own_y_vals = own_data['Plot Vals [cm]']['y']


def plot_contour(ax, x_vals=xplotVals, y_vals=yplotVals, linestyle='-.'):
    ax3 = ax.twinx()
    ax3.plot(x_vals, y_vals, color='grey', linestyle=linestyle)
    ax3.set_ylim([min(y_vals) * .9, max(y_vals) * 1.4])
    ax3.get_yaxis().set_visible(False)


# Heat Transfer Coefficients
_, ax = plt.subplots()
plot_contour(ax)
ax2 = ax.twinx()
ax.plot(xplotVals, np.array(hgvals) / 1000, color='r', label='Hot Gas D.', linestyle='-.')
ax2.plot(xplotVals, np.array(hcvals) / 1000, color='b', label='Coolant D.', linestyle='-.')

plot_contour(ax, x_vals=own_x_vals, y_vals=own_y_vals, linestyle='-')
variable_name = 'Heat-Transfer Coefficient [W/(K*m2]'
axes = (ax, ax2)
for name, color, axis in zip(['Hot Gas', 'Coolant'], ['r', 'b'], axes):
    values = [hc * 1e-3 for hc in own_data[variable_name][name]]
    axis.plot(own_x_vals, values, color=color, label=name)
    axis.set_ylabel(name + r' Heat-Transfer Coefficient [$kW/(K\cdot m^2$)]')

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height * .9])
multi_legend((ax, ax2), bbox_to_anchor=(.5, 1.1), loc='lower center', ncol=4)
ax.set_xlabel(r'$x$ coordinate [cm]')
ax.set_title(variable_name.split('[')[0])
plt.show()

# Coolant State
_, ax = plt.subplots()
plot_contour(ax)
plot_contour(ax, x_vals=own_x_vals, y_vals=own_y_vals, linestyle='-')
ax2 = ax.twinx()
for variable, axis, color, vals in zip(['T', 'T0', 'p', 'p0'],
                                       [ax, ax, ax2, ax2],
                                       ['orange', 'red', 'lightblue', 'blue'],
                                       [Tvals, T0vals, pvals, p0vals]):
    axis.plot(own_x_vals, own_data['Coolant State'][variable], label=variable, color=color, linestyle='-')
    axis.plot(xplotVals, vals, label=variable + ' D.', color=color, linestyle='-.')

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height * .9])
multi_legend((ax, ax2), bbox_to_anchor=(.5, 1.1), loc='lower center', ncol=4)
ax.set_ylabel(r'Temperature [$K$]')
ax2.set_ylabel(r'Pressure [$Pa$]')
ax.set_xlabel(r'$x$ coordinate [cm]')
ax.set_title('Coolant Bulk State')
plt.show()

# Wall Temps
_, ax = plt.subplots()

ax2 = ax.twinx()
ax2.set_ylabel(r'Heat Flux [$MW/m^2$]')

plot_contour(ax)
plot_contour(ax, x_vals=own_x_vals, y_vals=own_y_vals, linestyle='-')

hot_tws = own_data['Temperature [K]']['Hot SideWall']
cold_tws = own_data['Temperature [K]']['Cold SideWall']
heat_fluxs = [q * 1e-6 for q in own_data['Heat Flux [W/m2]']['Total']]

ax2.plot(own_x_vals, heat_fluxs, color='darkorange', label='Heat Flux')
ax2.plot(xplotVals, np.array(qvals) * 1e-6, color='darkorange', label='Heat Flux D.', linestyle='-.')

ax.plot(own_x_vals, hot_tws, color='r', label='Hot Side')
ax.plot(xplotVals, Twvals, color='r', label='Hot Side D.', linestyle='-.')

ax.plot(own_x_vals, cold_tws, color='b', label='Cold Side')
ax.plot(xplotVals, TwChannelvals, color='b', label='Cold Side D.', linestyle='-.')

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height * .9])
multi_legend((ax, ax2), bbox_to_anchor=(.5, 1.1), loc='lower center', ncol=3)
ax.set_ylim([min(TwChannelvals) * .7, max(Twvals) * 1.2])
ax.set_ylabel(r'Temperature [$K$]')
ax.set_xlabel(r'$x$ coordinate [cm]')
ax.set_title('Wall Temperatures')
plt.show()

# Geometry
_, ax = plt.subplots()

ax2 = ax.twinx()
ax2.set_ylabel(r'Radius [$cm$]')
ax2.plot(xplotVals, yplotVals, label='Radius D.', color='black')
ax2.plot(own_x_vals, own_y_vals, label='Radius', color='black')

geo_data = own_data['Channel Geometry [mm]']

ax.plot(xplotVals, np.array(channelHeightvals) * 1e3, color='blue', label='Height D.', linestyle='-.')
ax.plot(own_x_vals, geo_data['Height'] , color='blue', label='Height')
ax.plot(xplotVals, np.array(channelWidthvals) * 1e3, color='g', label='Width D.', linestyle='-.')
ax.plot(own_x_vals, geo_data['Width'] , color='g', label='Width')
ax.plot(xplotVals, np.array(channelDhvals) * 1e3, color='r', label='Eq. Dia. D.', linestyle='-.')
ax.plot(own_x_vals, geo_data['Diameter'], color='r', label='Eq. Dia.')


box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height * .9])
multi_legend((ax, ax2), bbox_to_anchor=(.5, 1.1), loc='lower center', ncol=4)
ax2.set_ylim([min(yplotVals)*.9, max(yplotVals)*1.5])
ax.set_ylabel(r'Length [$mm$]')
ax.set_xlabel(r'$x$ coordinate [cm]')
ax.set_title('Channel Geometry')
plt.show()