from rocketcea.cea_obj import CEA_Obj

pa_to_psia = 1.4503773800721814532099090803672E-4
ft_to_m = 0.3048


# noinspection PyPep8Naming
def get_cstar_cf(chamber_pressure, mixture_ratio=2.45, fuel_name='RP1_NASA', ox_name='LO2_NASA',
                 exit_pressure=0.002E6, isfrozen=0):
    p_ratio = chamber_pressure / exit_pressure
    chamber_pressure = chamber_pressure * pa_to_psia
    exit_pressure = exit_pressure * pa_to_psia
    frozen, frozenAtThroat = (1, 1) if isfrozen else (0, 0)
    cea = CEA_Obj(oxName=ox_name, fuelName=fuel_name)
    eps = cea.get_eps_at_PcOvPe(Pc=chamber_pressure, MR=mixture_ratio, PcOvPe=p_ratio, frozen=frozen,
                                frozenAtThroat=frozenAtThroat)
    _, c_star_ft, _ = cea.get_IvacCstrTc(Pc=chamber_pressure, MR=mixture_ratio, eps=eps, frozen=frozen,
                                         frozenAtThroat=frozenAtThroat)
    c_star_m = c_star_ft * ft_to_m
    if isfrozen:
        c_f = cea.getFrozen_PambCf(Pamb=exit_pressure, Pc=chamber_pressure, MR=mixture_ratio, eps=eps,
                                   frozenAtThroat=frozenAtThroat)[0]
    else:
        c_f = cea.get_PambCf(Pamb=exit_pressure, Pc=chamber_pressure, MR=mixture_ratio, eps=eps)[0]
    return c_star_m, c_f

if __name__ == '__main__':
    print(get_cstar_cf(10E6))