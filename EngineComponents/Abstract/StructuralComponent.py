from dataclasses import dataclass
from EngineComponents.Abstract.Material import Material


@dataclass
class StructuralComponent:
    structure_material: Material
    safety_factor: float  # [-]
