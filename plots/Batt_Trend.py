import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from itertools import tee


with open(r'C:\Users\rvand\PycharmProjects\ThesisEP2\data\Data\Tiede_2022_Batt_Trend.csv', 'r') as f:
    df = pd.read_csv(f)
    df = df.drop([0])
    df = df.astype(float)


def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

fig, ax = plt.subplots()

for (x,y) in pairwise(df):
    if 'Unnamed' in y:
        if x == 'Trend':
            ax.plot(list(df[x]), list(df[y]), label='Trend of Best', linestyle='--', color='black')
        else:
            ax.plot(list(df[x]), list(df[y]), label=f'{x}', linestyle='', marker='o')


spacing = 50
minorLocator = MultipleLocator(spacing)
ax.yaxis.set_minor_locator(minorLocator)
spacing = 10
minorLocator = MultipleLocator(spacing)
ax.xaxis.set_minor_locator(minorLocator)
ax.grid(which = 'both')

ax.set_xlabel(r'Time [Year]')
ax.set_ylabel(r'Specific Energy [Wh/kg]')
ax.set_ylim(0,500)
ax.set_xlim(1940, 2030)
plt.legend()
plt.savefig('Tiede_Batt_Trend.png', dpi=800)
plt.show()