from dataclasses import dataclass


@dataclass
class Splitter:
    input_mass_flow: float
    split_ratios: tuple

    def __post_init__(self):
        if sum(self.split_ratios) != 1:
            raise ValueError('The sum of all split ratios should be one')
        for i, split_ratio in enumerate(self.split_ratios):
            setattr(self,
                    f'output_mass_flow_{i + 1}',
                    split_ratio * self.input_mass_flow)


@dataclass
class Merger:
    input_mass_flows: tuple

    @property
    def output_mass_flow(self):
        return sum(input_mass_flows)
