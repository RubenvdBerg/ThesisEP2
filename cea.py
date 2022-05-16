from rocketcea.cea_obj import CEA_Obj

pa_to_psia = 1.4503773800721814532099090803672E-4
ft_to_m = 0.3048
rankine_to_kelvin = 5 / 9


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


def get_full_output(chamber_pressure, mixture_ratio=2.45, fuel_name='RP1_NASA', ox_name='LO2_NASA',
                    exit_pressure=0.002E6, isfrozen=0):
    p_ratio = chamber_pressure / exit_pressure
    chamber_pressure = chamber_pressure * pa_to_psia
    frozen, frozenAtThroat = (1, 1) if isfrozen else (0, 0)
    cea = CEA_Obj(oxName=ox_name, fuelName=fuel_name)
    return cea.get_full_cea_output(Pc=chamber_pressure, MR=mixture_ratio, PcOvPe=p_ratio, frozen=frozen,
                                   frozenAtThroat=frozenAtThroat, short_output=1, show_transport=1,
                                   show_mass_frac=False)


# (1764.9487103940155, 1.8773714044430267)


def get_cea_values(chamber_pressure, mixture_ratio=2.45, fuel_name='RP1_NASA', ox_name='LO2_NASA',
                   exit_pressure=0.002E6, isfrozen=0):
    p_ratio = chamber_pressure / exit_pressure
    chamber_pressure = chamber_pressure * pa_to_psia
    exit_pressure = exit_pressure * pa_to_psia
    frozen, frozenAtThroat = (1, 1) if isfrozen else (0, 0)
    cea = CEA_Obj(oxName=ox_name, fuelName=fuel_name, )
    kwargs = {'Pc': chamber_pressure, 'MR': mixture_ratio, 'frozen': frozen, 'frozenAtThroat': frozenAtThroat}

    # print(cea.get_full_cea_output(**kwargs, PcOvPe=p_ratio, eps=None, short_output=1))
    eps = cea.get_eps_at_PcOvPe(**kwargs, PcOvPe=p_ratio)
    kwargs['eps'] = eps
    kwargs2 = {i: kwargs[i] for i in kwargs if i != 'frozenAtThroat'}
    _, c_star_ft, _ = cea.get_IvacCstrTc(**kwargs)
    c_star_m = c_star_ft * ft_to_m

    if isfrozen:
        c_f = cea.getFrozen_PambCf(**kwargs, Pamb=exit_pressure)[0]
    else:
        c_f = cea.get_PambCf(Pamb=exit_pressure, Pc=chamber_pressure, MR=mixture_ratio, eps=eps)[0]
    temps_R = cea.get_Temperatures(**kwargs)
    temps_K = [temp * rankine_to_kelvin for temp in temps_R]
    transport_cc = cea.get_Chamber_Transport(**kwargs2)
    transport_th = cea.get_Throat_Transport(**kwargs2)
    transport_ex = cea.get_Exit_Transport(**kwargs2)
    transport = list(zip(transport_cc, transport_th, transport_ex))
    transport = transport_unit_fixer(transport)
    mw_gamma_cc = cea.get_Chamber_MolWt_gamma(Pc=chamber_pressure, MR=mixture_ratio, eps=eps)
    mw_gamma_th = cea.get_Throat_MolWt_gamma(**kwargs2)
    mw_gamma_ex = cea.get_exit_MolWt_gamma(**kwargs2)
    mw_gamma = list((zip(mw_gamma_cc, mw_gamma_th, mw_gamma_ex)))
    mw_gamma = mw_gamma_unit_fixer(mw_gamma)
    return c_star_m, c_f, temps_K, transport, mw_gamma, eps


def mw_gamma_unit_fixer(mw_gamma_list):
    mw_gamma_list[0] = [x * 1E-3 for x in mw_gamma_list[0]]  # g/mol to kg/mol (not SI, but required)
    mw_gamma_list[1] = [x * 1 for x in mw_gamma_list[1]]
    return mw_gamma_list


def transport_unit_fixer(transport_list):
    transport_list[0] = [x * 4.184 * 1E3 for x in transport_list[0]]   # From cal/gK to J/(kg K)
    transport_list[1] = [x * 1E-4 for x in transport_list[1]]   # From miliPoise to Pascal second
    transport_list[2] = [x * 4.184*1E-1 for x in transport_list[2]]  # From mcal/(K s cm) to J/(K s m)
    transport_list[3] = [x * 1 for x in transport_list[3]]   # From [-] to [-]
    return transport_list


if __name__ == '__main__':
    chamber_pressure, mixture_ratio = 1E6, 2.45
    # print(get_cea_values(chamber_pressure, mixture_ratio))
    cea_output = get_cea_values(chamber_pressure, mixture_ratio)
    _, _, _, transport, mw_gamma = cea_output
    print(transport)
    print(mw_gamma)
    # time_0 = perf_counter()
    # for i in range(100):
    #     get_cea_values(chamber_pressure, mixture_ratio)
    # time_zero = perf_counter() - time_0
    # time_1 = perf_counter()
    # for i in range(100):
    #     get_full_output(chamber_pressure, mixture_ratio)
    # time_one = perf_counter() - time_1
    # print(time_zero*1e3, time_one*1e3)
