from BaseEngineCycle.Nozzle import BellNozzle, ConicalNozzle
from BaseEngineCycle.ThrustChamber import ThrustChamber
from BaseEngineCycle.CombustionChamber import Injector, CombustionChamber
from math import pi, radians

# Variables from LE5A engine as given in Kakuma2000
exit_radius = 2.
expansion_ratio = 130
exit_area = pi * exit_radius ** 2.
throat_area = exit_area / expansion_ratio
throat_radius = (throat_area / pi) ** .5
chamber_radius = throat_radius * 2
cc = CombustionChamber(1, 1, 1, throat_area=throat_area, combustion_chamber_pressure=3.98E6, propellant_mix='LOX/LH2')

inj = Injector(1, 1, 1, 1, 1, False)
bnozzle = BellNozzle(throat_area=throat_area, expansion_ratio=expansion_ratio, chamber_radius=cc.radius,
                     conv_chamber_bend_ratio=1.0, conv_throat_bend_ratio=.8, conv_half_angle=radians(30),
                     div_throat_half_angle=radians(35), div_exit_half_angle=radians(5))
cnozzle = ConicalNozzle(throat_area=throat_area, expansion_ratio=expansion_ratio, chamber_radius=cc.radius,
                        conv_chamber_bend_ratio=1.0, conv_throat_bend_ratio=.8, conv_half_angle=radians(30),
                        div_throat_half_angle=radians(15))
tc = ThrustChamber(nozzle=bnozzle, chamber=cc, injector=inj, heat_capacity_ratio=1.14)
tc.show_contour()
tc2 = ThrustChamber(nozzle=cnozzle, chamber=cc, injector=inj, heat_capacity_ratio=1.14)
tc2.show_contour()
print(ConicalNozzle(throat_area=.1 ** 2 * pi,
                    expansion_ratio=50,
                    chamber_radius=.2,
                    conv_half_angle=radians(20),
                    div_throat_half_angle=radians(30),
                    conv_chamber_bend_ratio=1,
                    conv_throat_bend_ratio=.8).div_length)
print(BellNozzle(throat_area=.1 ** 2 * pi,
                 expansion_ratio=50,
                 chamber_radius=.2,
                 conv_half_angle=radians(20),
                 div_throat_half_angle=radians(30),
                 div_exit_half_angle=radians(5),
                 conv_chamber_bend_ratio=1,
                 conv_throat_bend_ratio=.8).div_length)
