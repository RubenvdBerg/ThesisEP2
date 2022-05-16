from scipy.interpolate import interp1d
from rocketcea.cea_obj import CEA_Obj
from base_gg_cycle import GasGeneratorCycle
from base_ep_cycle import ElectricPumpCycle
from base_engine_cycle import Nozzle
from arguments import base_arguments, ep_arguments, gg_arguments
# temperature = [297,
#                600,
#                800,
#                900,
#                1000,
#                1050,
#                1100,
#                1150,
#                1200,
#                1300,
#                1373]
# ultimate_tensile_strength = [
#     733,
#     689,
#     617,
#     468,
#     273,
#     212,
#     154,
#     113,
#     79,
#     50,
#     27
# ]
# sigma_ult_function = interp1d(temperature, ultimate_tensile_strength)
# cea = CEA_Obj(oxName='LO2_NASA', fuelName='LH2_NASA')
# psia_to_pascal = 6894.75728
#
#
# def pascal_to_psia(pressure_in_pascal: float):
#     return pressure_in_pascal / psia_to_pascal
#
#
# def rankine_to_kelvin(temperature_in_rankine: float):
#     return temperature_in_rankine * 5 / 9
#
#
# is_frozen = True
# frozen, frozenAtThroat = (1, 1) if is_frozen else (0, 0)
# mmr = 5.1
# chamber_pressure = 10E6
# area_ratio = 45
# isp_vac, c_star, t_comb = cea.get_IvacCstrTc(Pc=pascal_to_psia(chamber_pressure), MR=mmr, eps=area_ratio, frozen=frozen,
#                                              frozenAtThroat=frozenAtThroat)
# temps = cea.get_Temperatures(Pc=pascal_to_psia(chamber_pressure), MR=mmr, eps=area_ratio, frozen=frozen,
#                              frozenAtThroat=frozenAtThroat)
# print(*(rankine_to_kelvin(temp) for temp in temps))
# is_frozen = False
# temps = cea.get_Temperatures(Pc=pascal_to_psia(chamber_pressure), MR=mmr, eps=area_ratio, frozen=frozen,
#                              frozenAtThroat=frozenAtThroat)
# print(*(rankine_to_kelvin(temp) for temp in temps))

if __name__ == '__main__':
    ep = ElectricPumpCycle(**base_arguments,**ep_arguments)
    gg = GasGeneratorCycle(**base_arguments,**gg_arguments)
