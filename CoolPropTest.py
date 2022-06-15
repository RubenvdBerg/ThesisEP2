from CoolProp.CoolProp import PropsSI, FluidsList
# print(FluidsList())
mm = PropsSI('Hydrogen', 'molemass')
cp = PropsSI('C', 'T', 50., 'P', 1.150E6, 'Hydrogen')
print(cp*mm)