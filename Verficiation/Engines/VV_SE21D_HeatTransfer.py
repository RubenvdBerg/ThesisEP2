from CoolProp.CoolProp import PropsSI

h1 = PropsSI('H',
              'T', 32.375,
              'P', 12.109e6,
              'Hydrogen')
h2 = PropsSI('H',
              'T', 506.219,
              'P', 8.749e6,
              'Hydrogen')
dh = h2 - h1

P_q = dh * 16.394
print(P_q*1e-3)