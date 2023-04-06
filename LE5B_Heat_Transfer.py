from CoolProp.CoolProp import PropsSI

h1 = PropsSI('H', 'T', 21, 'P', 68.8e5, 'Hydrogen')
h2 = PropsSI('H', 'T', 409, 'P', 36.9e5, 'Hydrogen')
m_flow = .66
q = (h2 - h1)/m_flow
print(q*1e-6)