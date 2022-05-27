class Rocket:
    def __init__(self, density, thrust, specific_impulse):
        self.density = density
        self.thrust = thrust
        self.specific_impulse = specific_impulse

    @property
    def mass_flow(self):
        return self.thrust / self.specific_impulse

    @property
    def pump(self):
        return Pump(density=self.density, mass_flow=self.mass_flow)

    @property
    def deltav(self):
        return self.pump.power * self.specific_impulse


class Pump:
    def __init__(self, density, mass_flow):
        self.density = density
        self.mass_flow = mass_flow

    @property
    def power(self):
        return self.density * self.mass_flow


class BetterRocket:
    def __init__(self, pump: Pump, thrust: float, specific_impulse: float):
        self.pump = pump
        self.thrust = thrust
        self.specific_impulse = specific_impulse

    @property
    def mass_flow(self):
        return self.thrust / self.specific_impulse

    @property
    def deltav(self):
        return self.pump.power * self.specific_impulse


@dataclass
class EvenBetterRocket:
    pump: Pump
    thrust: float
    specific_impulse: float

    @property
    def mass_flow(self):
        return self.thrust / self.specific_impulse

    @property
    def deltav(self):
        return self.pump.power * self.specific_impulse
