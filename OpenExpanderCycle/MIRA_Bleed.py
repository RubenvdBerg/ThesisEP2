from OECycle import OpenExpanderCycle
from dataclasses import dataclass
from BaseEngineCycle.Pump import Pump
from BaseEngineCycle.SplitterMerger import Splitter, Merger


@dataclass
class MIRA(OpenExpanderCycle):
