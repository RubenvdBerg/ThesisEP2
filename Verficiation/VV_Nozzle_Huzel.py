from EngineComponents.Base.Nozzle import BellNozzle, ConicalNozzle
from EngineComponents.Base.CombustionChamber import CombustionChamber
from EngineComponents.Base.ThrustChamber import ThrustChamber
from math import radians

# def inch_to_meter(length_in_inch: float) -> float:
#     return length_in_inch * 0.0254
inch_to_meter = 0.0254
# inch_to_meter = 1
throat_area = 487 * inch_to_meter**2
expansion_ratio = 14
area_ratio_chamber_throat = 1.6


bnozzle = BellNozzle(throat_area=throat_area, expansion_ratio=expansion_ratio, area_ratio_chamber_throat=area_ratio_chamber_throat,
                     conv_chamber_bend_ratio=0, conv_throat_bend_ratio=1.5, conv_half_angle=radians(20),
                     div_throat_half_angle=radians(27.4), div_exit_half_angle=radians(9.8))
cnozzle = ConicalNozzle(throat_area=throat_area, expansion_ratio=expansion_ratio, area_ratio_chamber_throat=area_ratio_chamber_throat,
                        conv_chamber_bend_ratio=0, conv_throat_bend_ratio=1.5, conv_half_angle=radians(20),
                        div_throat_half_angle=radians(15))
cc = CombustionChamber(1, 1, 1, throat_area=throat_area,
                       characteristic_length=45*inch_to_meter,
                       combustion_chamber_pressure=1E6,
                       area_ratio_chamber_throat=area_ratio_chamber_throat,
                       convergent_volume_estimate=cnozzle.conv_volume_estimate)
tc = ThrustChamber(injector=1, chamber=cc, nozzle=bnozzle, heat_capacity_ratio=1.14)
tc2 = ThrustChamber(injector=1, chamber=cc, nozzle=cnozzle, heat_capacity_ratio=1.14)
tc2.show_contour()
tc.show_contour()

print('Parameter [unit]   Expected  Actual ')
print(f'CC Length [inch]  \t 18.17 \t {cc.length / inch_to_meter:.2f}')
print(f'CC Radius [inch]  \t 15.75 \t {cc.radius/ inch_to_meter:.2f}')
print(f'Throat Radius [inch] 12.45 \t {bnozzle.throat_radius / inch_to_meter:.2f}')
print(f'Exit Radius [inch]\t 46.70 \t {bnozzle.exit_radius/ inch_to_meter:.2f}')
print(f'Conv Length [inch]\t 12.40 \t {bnozzle.conv_length/ inch_to_meter:.2f}')
print(f'Div Length [inch] \t 102.4 \t {bnozzle.div_length/ inch_to_meter:.1f}')
print(f'CDiv Length [inch]\t 128.4 \t {cnozzle.div_length/ inch_to_meter:.1f}')
print(f'Conv Volume [inch3]\t 7760 \t {cnozzle.conv_volume_estimate/ inch_to_meter**3:.0f}')
print(f'CConv Volume[inch3]\t 7760 \t {cnozzle.conv_volume_estimate/ inch_to_meter**3:.0f}')
