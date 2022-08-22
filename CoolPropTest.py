from CoolProp.CoolProp import PropsSI, FluidsList
# print(FluidsList())
mm = PropsSI('Hydrogen', 'molemass')
cp = PropsSI('C', 'T', 50., 'P', 1.150E6, 'Hydrogen')
h1 = PropsSI('H', 'T', 32.375, 'P', 12.019E6, 'Hydrogen')
h2 = PropsSI('H', 'T', 506.219, 'P', 8.749E6, 'Hydrogen')
print(PropsSI('C', 'T', 506.452, 'P', 8.311E6, 'Hydrogen'))
print(h1, h2, h2-h1)
h_flow = 16.394
print(f'{(h2-h1)*h_flow*1e-6:.3f} MW')
m1 = PropsSI('H', 'T', 117, 'P', 7.68E6, 'Methane')
m2 = PropsSI('H', 'T', 490, 'P', 6.66E6, 'Methane')
m3 = PropsSI('H', 'T', 600, 'P', 6.66E6, 'Methane')
m_flow1 = 5.03
m_flow2 = .87
print(f'{(m2-m1)*m_flow1*1e-6:.3f} MW')
print(f'{(m3-m2)*m_flow2*1e-6:.3f} MW')


def in_out(temps, pressures, prop_name):
    t_in, t_out = temps
    p_in, p_out = pressures
    cp_in = PropsSI('C', 'T', t_in, 'P', p_in, prop_name)
    cp_out = PropsSI('C', 'T', t_out, 'P', p_out, prop_name)
    y_in = cp_in / PropsSI('O', 'T', t_in, 'P', p_in, prop_name)
    y_out = cp_out / PropsSI('O', 'T', t_out, 'P', p_out, prop_name)
    rho_in = PropsSI('D', 'T', t_in, 'P', p_in, prop_name)
    rho_out = PropsSI('D', 'T', t_out, 'P', p_out, prop_name)
    print(f'Cp [kJ/kgK], Gamma [-], rho [kg/m3] for {prop_name}:')
    print(f'   \t IN \t\t OUT')
    print(f'Cp \t {cp_in:.1f} \t {cp_out:.1f}')
    print(f'Y  \t {y_in:.4f} \t {y_out:.4f}')
    print(f'D  \t {rho_in:.3f} \t {rho_out:.3f}')

t_in_h2, t_out_h2 = 506.452, 369.677
p_in_h2, p_out_h2 = 8.311E6, .3E6
in_out((t_in_h2, t_out_h2), (p_in_h2, p_out_h2), 'Hydrogen')
t_in_m, t_out_m = 600, 470.5
p_in_m, p_out_m = 6.66E6,.46E6
in_out((t_in_m, t_out_m), (p_in_m, p_out_m), 'Methane')

print(PropsSI('D', 'T', 115, 'P', .3E6, 'Methane'))
print(PropsSI('D', 'T', 90, 'P', .3E6, 'Oxygen'))