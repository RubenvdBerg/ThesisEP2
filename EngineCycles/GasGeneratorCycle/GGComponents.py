import warnings

from EngineCycles.BaseEngineCycle.Structure import PressureStructure
from EngineCycles.BaseEngineCycle.FlowState import ManualFlowState, DefaultFlowState, FlowState
from EngineCycles.Functions.CEAFunctions import get_cea_dict_gg, get_gas_generator_mmr
from dataclasses import dataclass, field
import scipy.constants as constants
from typing import Optional

@dataclass
class GasGenerator(PressureStructure):
    """Handles mass estimation of the Gas Generator and merging of inlet flows. mass_mixture_ratio is required for
    clarity, so it is available as attribute for other components, as well as inheritance

    Manually input gas generator conditions or calculate them with CEA. If calculated with CEA the combustion temperature is matched to the turbine temp limit.
    """
    stay_time: float = 0  # [s]
    pressure: float = 0  # [Pa]
    turbine_temp_limit: float = 0  # [K]
    oxidizer_inlet_flow_state: FlowState = DefaultFlowState()
    fuel_inlet_flow_state: FlowState = DefaultFlowState()
    is_frozen: bool = True
    oxidizer_name: str = None
    fuel_name: str = None
    specific_heat_capacity: Optional[float] = None  # [J/(kgK)]
    heat_capacity_ratio: Optional[float] = None  # [-]
    molar_mass: Optional[float] = None  # [kg/mol]
    mass_mixture_ratio: Optional[float] = None  # [-]
    _chamber_temperature: float = field(init=False, repr=False, default=None)

    def __post_init__(self):
        self.set_cea()
        self.check_temp()

    @property
    def cea_kwargs(self):
        return {'fuelName': self.fuel_name,
                'oxName': self.oxidizer_name,
                'Pc': self.pressure,
                'frozen': 1 if self.is_frozen else 0,
                'frozenAtThroat': 1 if self.is_frozen else 0}

    def set_cea(self):
        attributes = {'heat_capacity_ratio': 'y_cc', 'specific_heat_capacity': 'cp_cc', 'molar_mass': 'mm_cc', '_chamber_temperature': 'T_C'}
        if any(getattr(self, name) is None for name in attributes):
            cea_dict = get_cea_dict_gg(MR=self.mass_mixture_ratio, **self.cea_kwargs)
            for attribute_name, value in attributes.items():
                if getattr(self, attribute_name) is None:
                    setattr(self, attribute_name, cea_dict[value])

    def check_temp(self):
        if self._chamber_temperature > self.turbine_temp_limit*1.05:
            warnings.warn(f'The manually provided mixture ratio [{self.mass_mixture_ratio}] for the gas generator leads to a combustion temperature [{self._chamber_temperature}] that is higher than the turbine inlet temperature [{self.turbine_temp_limit}].')

    @property
    def outlet_temperature(self):
        """Gas generator operated at limit of turbine, temperature wise"""
        return self.turbine_temp_limit

    @property
    def outlet_pressure(self):
        return self.pressure

    @property
    def outlet_mass_flow(self) -> float:
        return sum((self.oxidizer_inlet_flow_state.mass_flow, self.fuel_inlet_flow_state.mass_flow))

    @property
    def outlet_flow_state(self) -> FlowState:
        return ManualFlowState(propellant_name='ExhaustGas',
                               temperature=self.outlet_temperature,
                               pressure=self.outlet_pressure,
                               mass_flow=self.outlet_mass_flow,
                               type='combusted',
                               _specific_heat_capacity=self.specific_heat_capacity,
                               _heat_capacity_ratio=self.heat_capacity_ratio,
                               _molar_mass=self.molar_mass, )

    @property
    def mass(self):  # Overwrite Structure.mass()
        return (self.safety_factor * 3 / 2 * self.material_density / self.yield_strength
                * self.stay_time * self.pressure / self.gas_density * self.outlet_mass_flow)

    @property
    def gas_density(self):
        """Estimate gas generator exhaust gas density with ideal gas law"""
        return self.pressure / (self.outlet_flow_state.specific_gas_constant * self.turbine_temp_limit)
