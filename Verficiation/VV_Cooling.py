from BaseEngineCycle.Cooling2 import CoolingChannels

# for propellant_name in ('LH2_NASA', 'LO2_NASA', 'RP1_NASA'):
#     coolant = Coolant(propellant_name=propellant_name)
#     channels = CoolingChannels(coolant=coolant, total_heat_transfer=5E9, outlet_pressure=1E6, mass_flow=1000)
channels = CoolingChannels(propellant_name='LH2_NASA', total_heat_transfer=5E9, outlet_pressure=1E6, mass_flow=1000)
print()
