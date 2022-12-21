from EngineCycles.BaseEngineCycle.EngineCycle import EngineCycle
from EngineCycles.GasGeneratorCycle.GGCycle import GasGeneratorCycle, GasGeneratorCycle_DoubleTurbine
from EngineCycles.ElectricPumpCycle.EPCycle import ElectricPumpCycle
from EngineCycles.OpenExpanderCycle.OECycle import CoolantBleedCycle, OpenExpanderCycle, OpenExpanderCycle_DoubleTurbine
from EngineCycles.BaseEngineCycle.FlowState import FlowState, ManualFlowState
from EngineCycles.BaseEngineCycle.Pump import Pump
from EngineCycles.BaseEngineCycle.Cooling import CoolingChannelSection
from EngineCycles.ElectricPumpCycle.EPComponents import ElectricComponent, Battery
from EngineCycles.BaseOpenCycle.Turbine import Turbine
from numpy import log10, isclose
from typing import Optional
from PIL import Image, ImageDraw, ImageFont


def format_si(value: float, unit: str, digits: int = 5, force_prefix: Optional[dict] = None):
    if force_prefix is None:
        force_prefix = {'g/s': 'k', 'Pa': 'M', 'W': 'M'}

    si_prefixes = ('Y', 'Z', 'E', 'P', 'T', 'G', 'M', 'k', '', 'd', 'c', 'm', '\u03BC', 'n', 'p', 'f', 'a', 'z', 'y')
    si_index = 8
    n_before_comma = log10(abs(value)) // 1 + 1
    if force_prefix is not None and unit in force_prefix:
        prefix = force_prefix[unit]
        x = si_index - si_prefixes.index(prefix)
    else:
        x = round(n_before_comma / 3 - 1)
        if (1 > x > -2):
            x = 0
        si_index -= x
        prefix = si_prefixes[si_index]
    value *= 10 ** (-3 * x)
    decimals = int(digits - (n_before_comma - (3 * x)))
    format = min(digits, decimals)
    if digits == decimals:
        format = digits - 1
    val = f'{value:.{format}f}'
    if len(val) > digits + 1:
        val = f'{value:.{format - 1}f}'
    return rf'{val} {prefix}{unit}'


def make_schematic(engine: EngineCycle):
    switcher = {
        ElectricPumpCycle: (get_ep_components_coordinates, 'EP'),
        GasGeneratorCycle: (get_gg_components_coordinates, 'GG'),
        CoolantBleedCycle: (get_cb_components_coordinates, 'CB'),
        OpenExpanderCycle: (get_oe_components_coordinates, 'OE'),
        GasGeneratorCycle_DoubleTurbine: (get_gg2_components_coordinates, 'GG2'),
        OpenExpanderCycle_DoubleTurbine: (get_oe2_components_coordinates, 'GG2'),
    }
    for base_class, value in switcher.items():
        if issubclass(type(engine), base_class):
            get_comps_coords, name = value
    comps, coords = get_comps_coords(engine)
    fontsize = 42
    myfont = ImageFont.truetype(r'C:\Users\rvand\PycharmProjects\ThesisEP2\plots\Imaging\DejaVuSans.ttf', fontsize)
    image_path = rf'C:\Users\rvand\PycharmProjects\ThesisEP2\plots\Imaging\Schematics\{name}_Cycle.png'
    strings = tuple(format_values(comps))
    with Image.open(image_path) as img:
        drawer = ImageDraw.Draw(img)
        for coord, string_row in zip(coords, strings):
            coord_now = list(coord)
            for string in string_row:
                drawer.text(coord_now, string, fill=(0, 0, 0), font=myfont)
                coord_now[1] = coord_now[1] + int(fontsize * 1.15)
        img.show()


def get_ambient_pressure_string(ambient_pressure: Optional[float]):
    if ambient_pressure is None:
        p_a = u'p\u2091'
    elif isclose(ambient_pressure, 0):
        p_a = '0.0'
    elif isclose(ambient_pressure, 101325):
        p_a = '1 atm'
    else:
        p_a = format_si(ambient_pressure, 'Pa')
    return u'p\u2090' + f' = {p_a}'


def get_base_comps_coords(engine: EngineCycle, x_y1: tuple, x_y2: tuple):
    from math import isclose
    ambient_string = get_ambient_pressure_string(engine.ambient_pressure)

    components = [
        engine.fuel_pump.inlet_flow_state,
        engine.fuel_pump,
        engine.cooling_channel_section.outlet_flow_state,
        engine.fuel_pump.outlet_flow_state,
        engine,
        engine.oxidizer_pump.inlet_flow_state,
        'placeholder',
        engine.oxidizer_pump,
        engine.oxidizer_pump.outlet_flow_state,
        ManualFlowState(propellant_name='CombustionGas',
                        temperature=engine.combustion_temperature,
                        pressure=engine.combustion_chamber_pressure,
                        mass_flow=engine.chamber_mass_flow,
                        type='combusted'),
        (engine.chamber_thrust, engine.chamber_specific_impulse, engine.expansion_ratio),
        engine.propellant_mix_name.split('/')[0],  # Oxidizer
        engine.propellant_mix_name.split('/')[1],  # Fuel
        f"{ambient_string: >20}",
    ]

    dy3 = 193
    dy2 = 144
    x1, y1 = x_y1
    x2, y2 = x_y2
    coords = [
        (x1, y1),
        (x1, y1 + dy3),
        (x1, y1 + dy3 + dy2),
        (x1, y1 + 2 * dy3 + dy2),
        (x1, y1 + 3 * dy3 + dy2),
        (x2, y2),
        (x2, y2 + dy3),
        (x2, y2 + dy3 + dy2),
        (x2, y2 + dy3 + 2 * dy2),
        (x2, y2 + 2 * dy3 + 2 * dy2),
        (x2, y2 + 3 * dy3 + 2 * dy2),
        (x2 - 10, y2 - 450),
        (x2 - 960, y2 - 450),
        (x2 - 550, y2 + 1030), ]
    return components, coords


def get_ep_components_coordinates(engine: ElectricPumpCycle):
    ambient_string = get_ambient_pressure_string(engine.ambient_pressure)

    components = (
        engine.fuel_tank.outlet_flow_state,
        engine.fuel_pump,
        engine.battery_cooler.outlet_flow_state,
        engine.cooling_channel_section.outlet_flow_state,
        engine.post_fuel_pump_splitter.outlet_flow_state_chamber,
        engine,
        engine.oxidizer_pump.inlet_flow_state,
        engine.oxidizer_pump,
        engine.battery,
        engine.oxidizer_pump.outlet_flow_state,
        ManualFlowState(propellant_name='CombustionGas',
                        temperature=engine.combustion_temperature,
                        pressure=engine.combustion_chamber_pressure,
                        mass_flow=engine.chamber_mass_flow,
                        type='combusted'),
        (engine.chamber_thrust, engine.chamber_specific_impulse, engine.expansion_ratio),
        f"{ambient_string: >20}",
        engine.propellant_mix_name.split('/')[0],  # Oxidizer
        engine.propellant_mix_name.split('/')[1],  # Fuel
        engine.electric_motor,
        engine.inverter,
        format_si(engine.thrust, 'N'),
        format_si(engine.overall_specific_impulse, 's'),
    )

    dy1, dy2, dy3 = 80, 144, 193
    x1, y1 = 196, 840
    x2, y2 = x1 + (1295 - 116), y1 - dy3
    coordinates = (
        (x1 - 122, y1 - dy3 - 65),
        (x1, y1),
        (x1, y1 + dy2),
        (x1, y1 + dy2 + dy3),
        (x1, y1 + dy2 + dy3 * 2),
        (x1, y1 + dy2 + dy3 * 3),
        (x2, y2),
        (x2, y2 + dy3),
        (x2, y2 + dy2 + dy3),
        (x2, y2 + dy2 * 2 + dy3),
        (x2, y2 + dy2 * 2 + dy3 * 2),
        (x2, y2 + dy2 * 2 + dy3 * 3),
        (x1 + 554, y1 + 948),  # Ambient
        (x1 + 1194, y1 - 580),  # Ox
        (x1 + 123, y1 - 580),  # Fu
        (x1 + 389, y1 - 605),  # Em
        (x1 + 804, y1 - 605),  # Inv
        (x1 + 464, y1 - 780),
        (x1 + 929, y1 - 780)
    )

    return components, coordinates


def get_gg_components_coordinates(engine: GasGeneratorCycle):
    x_y2 = (1420, 640)
    base_components, base_coords = get_base_comps_coords(engine, x_y1=(315, 645), x_y2=x_y2)
    # Find and replace placeholder with battery component
    i = base_components.index('placeholder')
    base_components[i] = engine.turbine

    # Mass Flows after Splitters
    m_f_gg = engine.post_fuel_pump_splitter.outlet_flow_states['gg'].mass_flow
    m_o_gg = engine.post_oxidizer_pump_splitter.outlet_flow_states['gg'].mass_flow
    m_o_ch = engine.post_oxidizer_pump_splitter.outlet_flow_states['main'].mass_flow

    ambient_string = get_ambient_pressure_string(engine.ambient_pressure)

    components = (engine.fuel_pump.inlet_flow_state,
                  engine.fuel_pump,
                  engine.fuel_pump.outlet_flow_state,
                  format_si(m_f_gg * 1e3, 'g/s'),
                  engine.cooling_channel_section.outlet_flow_state,
                  engine,
                  engine.oxidizer_pump.inlet_flow_state,
                  engine.turbine,
                  engine.oxidizer_pump,
                  engine.oxidizer_pump.outlet_flow_state,
                  format_si(m_o_gg * 1e3, 'g/s'),
                  format_si(m_o_ch * 1e3, 'g/s'),
                  ManualFlowState(propellant_name='CombustionGas',
                                  temperature=engine.combustion_temperature,
                                  pressure=engine.combustion_chamber_pressure,
                                  mass_flow=engine.chamber_mass_flow,
                                  type='combusted'),
                  (engine.chamber_thrust, engine.chamber_specific_impulse, engine.expansion_ratio),
                  engine.propellant_mix_name.split('/')[0],  # Oxidizer
                  engine.propellant_mix_name.split('/')[1],  # Fuel
                  f"{ambient_string: >20}",
                  engine.gas_generator.outlet_flow_state,
                  (engine.secondary_exhaust.thrust, engine.secondary_exhaust.specific_impulse,
                   engine.secondary_exhaust.expansion_ratio),
                  engine.turbine.outlet_flow_state,
                  format_si(engine.thrust, 'N'),
                  format_si(engine.overall_specific_impulse, 's'),
                  )

    dy1 = 80
    dy2 = 144
    dy3 = 193
    x1, y1 = 320, 635
    x2, y2 = 1540, 520
    x3, y3 = 725, -60

    coords = ((x1, y1),
              (x1, y1 + dy3),
              (x1, y1 + dy3 + dy2),
              (x1, y1 + dy3 * 2 + dy2 - 5),
              (x1, y1 + dy3 * 2 + dy2 + dy1),
              (x1, y1 + dy3 * 3 + dy2 + dy1),
              (x2, y2),
              (x2, y2 + dy3),
              (x2, y2 + dy3 + dy2 + 3),
              (x2, y2 + dy3 + dy2 * 2),
              (x2, y2 + dy3 * 2 + dy2 * 2),
              (x2, y2 + dy3 * 2 + dy2 * 2 + dy1),
              (x2, y2 + dy3 * 2 + dy2 * 2 + dy1 * 2 + 5),
              (x2, y2 + dy3 * 3 + dy2 * 2 + dy1 * 2),
              (1415, 200),
              (460, 200),
              (870, 1670),
              (879, 753),
              (825, 50),
              (825, 50 + dy3),
              (x3, y3),
              (x3 + 465, y3),
              )
    coords = tuple((x, y + 120) for x, y in coords)
    return components, coords


def get_gg2_components_coordinates(engine: GasGeneratorCycle_DoubleTurbine):
    # Mass Flows after Splitters
    m_f_gg = engine.post_fuel_pump_splitter.outlet_flow_states['gg'].mass_flow
    m_o_gg = engine.post_oxidizer_pump_splitter.outlet_flow_states['gg'].mass_flow
    m_o_ch = engine.post_oxidizer_pump_splitter.outlet_flow_states['main'].mass_flow

    ambient_string = get_ambient_pressure_string(engine.ambient_pressure)

    components = (
        engine.fuel_pump.inlet_flow_state,
        engine.fuel_pump,
        engine.fuel_pump.outlet_flow_state,
        format_si(m_f_gg * 1e3, 'g/s'),
        engine.gas_generator.outlet_flow_state,
        engine.cooling_channel_section.outlet_flow_state,
        engine,
        # Second Column
        engine.oxidizer_pump.inlet_flow_state,
        engine.oxidizer_pump,
        engine.oxidizer_pump.outlet_flow_state,
        format_si(m_o_gg * 1e3, 'g/s'),
        format_si(m_o_ch * 1e3, 'g/s'),
        ManualFlowState(propellant_name='CombustionGas',
                        temperature=engine.combustion_temperature,
                        pressure=engine.combustion_chamber_pressure,
                        mass_flow=engine.chamber_mass_flow,
                        type='combusted'),
        (engine.chamber_thrust, engine.chamber_specific_impulse, engine.expansion_ratio),
        # Others
        (engine.fuel_secondary_exhaust.thrust, engine.fuel_secondary_exhaust.specific_impulse,
         engine.fuel_secondary_exhaust.expansion_ratio),
        engine.fuel_turbine.outlet_flow_state,
        engine.fuel_turbine,
        (engine.oxidizer_secondary_exhaust.thrust, engine.oxidizer_secondary_exhaust.specific_impulse,
         engine.oxidizer_secondary_exhaust.expansion_ratio),
        engine.oxidizer_turbine.outlet_flow_state,
        engine.oxidizer_turbine,
        engine.propellant_mix_name.split('/')[0],  # Oxidizer
        engine.propellant_mix_name.split('/')[1],  # Fuel
        f"{ambient_string: >20}",
        format_si(engine.thrust, 'N'),
        format_si(engine.overall_specific_impulse, 's'),
    )

    dy1 = 80
    dy2 = 144
    dy3 = 193
    x1, y1 = 320, 740
    x2, y2 = 1610, y1
    x3, y3 = 785, -60
    x4, y4 = 568, 105

    coords = (
        (x1, y1),
        (x1, y1 + dy3),
        (x1, y1 + dy3 + dy2),
        (x1, y1 + dy3 * 2 + dy2 - 10),
        (x1, y1 + dy3 * 2 + dy2 + dy1),
        (x1, y1 + dy3 * 3 + dy2 + dy1),
        (x1, y1 + dy3 * 4 + dy2 + dy1),
        (x2, y2),
        (x2, y2 + dy3),
        (x2, y2 + dy3 + dy2 + 3),
        (x2, y2 + dy3 * 2 + dy2),
        (x2, y2 + dy3 * 2 + dy2 + dy1),
        (x2, y2 + dy3 * 2 + dy2 + dy1 * 2 + 5),
        (x2, y2 + dy3 * 3 + dy2 + dy1 * 2),
        (x4, y4),  # Fuel Scnd Thrust
        (x4, y4 + dy3),  # Fuel Turb Flow
        (x4 - 120, y4 + dy3 * 2),  # Fuel Turbine
        (x4 + 672, y4),  # Ox Scnd Thrust
        (x4 + 672, y4 + dy3),  # Ox Turb Flow
        (x4 + 802, y4 + dy3 * 2),  # Ox Turbine
        (1765, 350),  # Ox Name
        (190, 350),  # Fuel Name
        (870, 1770),  # Ambient P
        (x3, y3),  # Thrust Total
        (x3 + 465, y3),  # Isp Total
    )
    coords = tuple((x, y + 120) for x, y in coords)
    return components, coords


def get_cb_components_coordinates(engine: CoolantBleedCycle):
    x_y2 = (1420, 640)
    base_components, base_coords = get_base_comps_coords(engine, x_y1=(315, 645), x_y2=x_y2)
    # Find and replace placeholder with battery component
    i = base_components.index('placeholder')
    base_components[i] = engine.turbine

    m_tu = engine.post_cooling_splitter.outlet_flow_states['turbine'].mass_flow
    m_ch = engine.post_cooling_splitter.outlet_flow_states['chamber'].mass_flow
    tau = f'{m_ch / (m_tu + m_ch): >9.2%}'
    components = (*base_components,
                  format_si(m_tu * 1e3, 'g/s', 4),
                  format_si(m_ch * 1e3, 'g/s', 4),
                  tau,
                  (engine.secondary_exhaust.thrust, engine.secondary_exhaust.specific_impulse,
                   engine.secondary_exhaust.expansion_ratio),
                  engine.turbine.outlet_flow_state,
                  format_si(engine.thrust, 'N'),
                  format_si(engine.overall_specific_impulse, 's'),
                  )
    dy1 = 80
    dy3 = 193
    x1, y1 = 940, 748
    x2, y2 = 825, 47
    x2, y2 = x_y2[0] - 595, x_y2[1] - 593
    x3, y3 = 745, -60
    coords = (*base_coords,
              (x1, y1),
              (x1, y1 + dy1),
              (x1, y1 + dy1 * 2),
              (x2, y2),
              (x2, y2 + dy3),
              (x3, y3),
              (x3 + 465, y3),
              )
    coords = tuple((x, y + 120) for x, y in coords)
    return components, coords


def get_oe_components_coordinates(engine: OpenExpanderCycle):
    x_y2 = (1620, 635)
    dy1 = 80
    dy3 = 193
    base_components, base_coords = get_base_comps_coords(engine, x_y1=(315, 635), x_y2=x_y2)
    # Find and replace placeholder with battery component
    i = base_components.index('placeholder')
    base_components[i] = engine.turbine
    # Change fuel_tank coordinates
    i_rp1 = base_components.index(engine.propellant_mix_name.split('/')[1])
    base_coords[i_rp1] = (435, 190)
    # Switchs positions of fuel pump and cooling outlet
    ccs_out, fp_out = base_components[2], base_components[3]
    base_components[3], base_components[2] = ccs_out, fp_out
    # Insert secondary fuel pump outlet and add new coordinates (for next component)
    base_components.insert(4, engine.secondary_fuel_pump.outlet_flow_state)
    old_x, old_y = base_coords[4]
    base_coords.insert(5, (old_x, old_y + dy3))

    m_tu = engine.post_cooling_splitter.outlet_flow_states['turbine'].mass_flow
    m_ch = engine.post_cooling_splitter.outlet_flow_states['chamber'].mass_flow
    tau = f'{m_ch / (m_tu + m_ch): >9.2%}'
    m_ch2 = engine.pre_cooling_splitter.outlet_flow_states['chamber'].mass_flow
    components = (*base_components,
                  (engine.secondary_exhaust.thrust, engine.secondary_exhaust.specific_impulse,
                   engine.secondary_exhaust.expansion_ratio),
                  engine.turbine.outlet_flow_state,
                  format_si(m_tu * 1e3, 'g/s', 4),
                  tau,
                  format_si(m_ch2 * 1e3, 'g/s', 5),
                  format_si(m_ch * 1e3, 'g/s', 5),
                  engine.pre_injection_merger.outlet_flow_state,
                  engine.secondary_fuel_pump,
                  format_si(engine.thrust, 'N'),
                  format_si(engine.overall_specific_impulse, 's'),
                  )
    x1, y1 = 1140, 758
    x2, y2 = x_y2[0] - 595, x_y2[1] - 593
    x3, y3 = 670, 1000
    x4, y4 = 825, -70
    coords = (*base_coords,
              (x2, y2),
              (x2, y2 + dy3),
              (x1, y1),
              (x1, y1 + dy1),
              (x3, y3),
              (x3, y3 + dy1),
              (x3, y3 + 250),
              (x3, 250),
              (x4, y4),
              (x4 + 465, y4),
              )
    coords = tuple((x, y + 130) for x, y in coords)
    return components, coords


def get_oe2_components_coordinates(engine: OpenExpanderCycle_DoubleTurbine):
    m_tu = engine.post_cooling_splitter.outlet_flow_states['turbine'].mass_flow
    m_ch = engine.post_cooling_splitter.outlet_flow_states['chamber'].mass_flow
    m_ch2 = engine.pre_cooling_splitter.outlet_flow_states['chamber'].mass_flow

    components = (
        # Fuel Col
        engine.fuel_pump.inlet_flow_state,
        engine.fuel_pump,
        engine.secondary_fuel_pump,
        engine.fuel_pump.outlet_flow_state,
        engine.secondary_fuel_pump.outlet_flow_state,
        engine.cooling_channel_section.outlet_flow_state,
        engine,
        # Ox Col
        engine.oxidizer_pump.inlet_flow_state,
        engine.oxidizer_pump,
        engine.oxidizer_turbine,
        engine.fuel_turbine,
        engine.oxidizer_pump.outlet_flow_state,
        ManualFlowState(propellant_name='CombustionGas',
                        temperature=engine.combustion_temperature,
                        pressure=engine.combustion_chamber_pressure,
                        mass_flow=engine.chamber_mass_flow,
                        type='combusted'),
        (engine.chamber_thrust, engine.chamber_specific_impulse, engine.expansion_ratio),
        # Fuel Col 2
        format_si(m_ch2 * 1e3, 'g/s', 5),
        format_si(m_ch * 1e3, 'g/s', 5),
        engine.pre_injection_merger.outlet_flow_state,
        # Top cols
        (engine.fuel_secondary_exhaust.thrust, engine.fuel_secondary_exhaust.specific_impulse,
         engine.fuel_secondary_exhaust.expansion_ratio),
        engine.fuel_turbine.outlet_flow_state,
        (engine.oxidizer_secondary_exhaust.thrust, engine.oxidizer_secondary_exhaust.specific_impulse,
         engine.oxidizer_secondary_exhaust.expansion_ratio),
        engine.oxidizer_turbine.outlet_flow_state,
        # Other
        format_si(m_tu * 1e3, 'g/s', 4),
        engine.propellant_mix_name.split('/')[0],  # Oxidizer
        engine.propellant_mix_name.split('/')[1],  # Fuel
        f"{ambient_string: >20}",
        format_si(engine.thrust, 'N'),
        format_si(engine.overall_specific_impulse, 's'),
    )
    dy1 = 80
    dy2 = 144
    dy3 = 193
    x1, y1 = 300, 550
    x2, y2 = 1400, y1
    x3, y3 = 500, y1 + dy2 + dy3 + 10
    x4, y4 = 600, 100
    x5, y5 = 700, 50
    coords = (
        (x1, y1),
        (x1, y1 + dy3),
        (x1, y1 + dy3 + dy2),
        (x1, y1 + dy3 + dy2 * 2),
        (x1, y1 + dy3 * 2 + dy2 * 2),
        (x1, y1 + dy3 * 3 + dy2 * 2),
        (x1, y1 + dy3 * 4 + dy2 * 2),
        (x2, y2),
        (x2, y2 + dy2),
        (x2, y2 + dy2 * 2),
        (x2, y2 + dy2 * 3),
        (x2, y2 + dy3 + dy2 * 3),
        (x2, y2 + dy3 * 2 + dy2 * 3),
        (x2, y2 + dy3 * 3 + dy2 * 3),
        # Fuel col 2
        (x3, y3),
        (x3, y3 + dy1),
        (x3, y3 + dy1 * 2),
        # Top cols
        (x4, y4),
        (x4, y4 + dy3),
        (x4 + 300, y4),
        (x4 + 300, y4 + dy3),
        # Other
        (100, 100),  # Turbine Mass Flow
        (100, 100),  # Ox Name
        (100, 100),  # Fu Name
        (100, 100),  # Pa
        (100, 100),  # Thrust
        (100, 100),  # Isp

    )

    return components, coords


def format_power_comp(power: float, efficiency: float, **kwargs):
    eta_f = ' >13.2f'
    return format_si(power, 'W', **kwargs), f'{efficiency:{eta_f}}'


def format_values(components: tuple):
    for component in components:
        if isinstance(component, FlowState):
            pressure = format_si(component.pressure, 'Pa')
            massflow = format_si(component.mass_flow * 1e3, 'g/s')
            temperature = format_si(component.temperature, 'K')
            yield massflow, pressure, temperature
        elif isinstance(component, Pump):
            power, efficiency = format_power_comp(component.power_required,
                                                  component.efficiency)
            yield power, efficiency
        elif isinstance(component, EngineCycle):
            power, efficiency = format_power_comp(component.total_heat_transfer,
                                                  component.expansion_ratio_end)
            yield power, efficiency

        elif isinstance(component, tuple):
            thrust = format_si(component[0], 'N')
            isp = format_si(component[1], 's')
            eps = format_si(component[2], '')
            yield thrust, isp, eps
        elif isinstance(component, ElectricComponent):
            digits = 5 if isinstance(component, Battery) else 5
            power, efficiency = format_power_comp(component.output_power,
                                                  component.electric_energy_efficiency,
                                                  digits=digits)
            yield power, efficiency

        elif isinstance(component, Turbine):
            power, efficiency = format_power_comp(component.power_required,
                                                  component.efficiency)
            yield power, efficiency
        elif isinstance(component, str):
            yield (component,)
        else:
            raise ValueError('Component type not in format list')


if __name__ == '__main__':
    import arguments as args
    from EngineCycles.Functions.IRTFunctions import get_characteristic_velocity
    from math import radians

    design_args = {'thrust': 1000e3,
                   'burn_time': 300,
                   'combustion_chamber_pressure': 10e6,
                   'is_frozen': True,
                   'ambient_pressure': None,
                   'exit_pressure_forced': 50000, }

    se_21d_kwargs = args.base_arguments_o | {
        'thrust': 1947.03e3,
        'combustion_chamber_pressure': 6.649e6,
        'fuel_initial_temperature': 21,
        'oxidizer_initial_temperature': 90,
        'expansion_ratio': 12.52,
        'mass_mixture_ratio': 5.5,
        'fuel_name': 'LH2_NASA',
        'burn_time': 100,
        'expansion_ratio_end_cooling': 5,
        'chamber_characteristic_length': 4.0,
        'fuel_pump_efficiency': .7,
        'secondary_fuel_pump_efficiency': .75,
        'oxidizer_pump_efficiency': .76,
        'fuel_initial_pressure': .3e6,
        'oxidizer_initial_pressure': .5e6,
        'turbine_maximum_temperature': 506.452,
        'turbopump_specific_power': 13.5E3,
        'area_ratio_chamber_throat': (.985 / 2) ** 2 / .286 ** 2,
        'ambient_pressure': 101325,
        'divergent_throat_half_angle': radians(15),
        'specific_impulse_correction_factor': .99,

        'oxidizer_turbine_efficiency': .45,
        'fuel_turbine_efficiency': .45,
        'oxidizer_shaft_mechanical_efficiency': .99,
        'fuel_shaft_mechanical_efficiency': .99,
        'oxidizer_exhaust_exit_pressure_forced': .04e6,
        'fuel_exhaust_exit_pressure_forced': .04e6,
        'oxidizer_turbine_outlet_pressure_forced': .3e6,
        'fuel_turbine_outlet_pressure_forced': .3e6,
        'oxidizer_secondary_specific_impulse_correction_factor': .98,
        'fuel_secondary_specific_impulse_correction_factor': .98,
        # 'turbine_pressure_ratio': 27.7033333,
        # 'exhaust_expansion_ratio': 1.655,
    }
    j2_kwargs = {'fuel_name': 'LH2_NASA',
                 'thrust': 1023e3,
                 'combustion_chamber_pressure': 5.4e6,
                 'expansion_ratio': 27.5,
                 'mass_mixture_ratio': 5.5,
                 'area_ratio_chamber_throat': 1.58,
                 'chamber_characteristic_length': 0.62,
                 'turbine_maximum_temperature': 922,
                 'gg_pressure': 4.7e6,
                 # 'fuel_pump_outlet_pressure': 8.62,
                 # 'oxidizer_pump_outlet_pressure': 7.64,
                 'fuel_pump_efficiency': 73e-2,
                 'oxidizer_pump_efficiency': 80e-2,
                 'fuel_turbine_pressure_ratio': 7.2,
                 'oxidizer_turbine_pressure_ratio': 2.65,
                 'fuel_turbine_efficiency': .6,
                 'oxidizer_turbine_efficiency': .47,
                 'burn_time': 475,
                 }
    j2_total_kwargs = args.base_arguments_o | j2_kwargs | {
        'ambient_pressure': 0,
        'fuel_exhaust_expansion_ratio': 4,
        'oxidizer_exhaust_expansion_ratio': 4,
        'convergent_throat_bend_ratio': 0.4,
        'convergent_chamber_bend_ratio': 0.5,
    }

    cycle_list = (
        # (ElectricPumpCycle, args.ep_arguments),
        # (GasGeneratorCycle, args.gg_arguments),
        # (CoolantBleedCycle, args.cb_arguments),
        # (OpenExpanderCycle, args.oe_arguments),
        (OpenExpanderCycle_DoubleTurbine, se_21d_kwargs),
        (GasGeneratorCycle_DoubleTurbine, j2_total_kwargs),
    )

    for Cycle, extra_args in cycle_list:
        complete_args = args.base_arguments | extra_args | design_args
        engine = Cycle(**complete_args)
        make_schematic(engine)
