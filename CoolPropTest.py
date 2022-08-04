from CoolProp.CoolProp import PropsSI, FluidsList
# print(FluidsList())
mm = PropsSI('Hydrogen', 'molemass')
cp = PropsSI('C', 'T', 50., 'P', 1.150E6, 'Hydrogen')
h1 = PropsSI('H', 'T', 32.375, 'P', 12.019E6, 'Hydrogen')
h2 = PropsSI('H', 'T', 506.219, 'P', 8.749E6, 'Hydrogen')
print(PropsSI('C', 'T', 506.452, 'P', 8.311E6, 'Hydrogen'))
print(h1, h2, h2-h1)
print(f'{(h2-h1)*16.394*1e-6:.3f} MW')