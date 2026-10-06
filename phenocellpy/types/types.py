import json
from dataclasses import dataclass, field

from phenocellpy.types.base import Validated
from phenocellpy.types.validators import (not_none, positive, non_negative, fraction, boolean, integer,
                                          non_empty_str, function, function_or_false, sequence)

from typing import Union, Tuple, List, Dict   

# Runtime state (time in phase, current volumes, ...) belongs to the Phase/Phenotype objects, not to these configs.
#
# Time-related values (`phase_duration`) are in the units of the owning PhenotypeConfig's `time_unit`. `dt` only
# lives in PhenotypeConfig; it is passed down to the phases when the runtime objects are built.
#
# Phase parameters left as `None` take the value of the phase class' `default_config()` (e.g., the volume change
# rates of :class:`Phases.Apoptosis`, the duration of :class:`Phases.NecrosisSwell`). If that is None too, the phase
# class or :class:`CellVolumes` computes its own default.


@dataclass()
class VolumeRatesConfig(Validated):
    """
    Volume change rates of a phase, in 1/`time_unit`. Used by the `CellVolume` model.

    :param cytoplasm_volume_change_rate: Change rate for the cytoplasmic volume. `>= 0`
    :param nuclear_volume_change_rate: Change rate for the nuclear volume. `>= 0`
    :param fluid_change_rate: Rate of change of the cell fluid part. `>= 0`
    :param calcification_rate: Rate of calcification of the cell. `>= 0`
    """
    cytoplasm_volume_change_rate: Union[float, None] = field(default=None, metadata={"checks": (non_negative,)})
    nuclear_volume_change_rate: Union[float, None] = field(default=None, metadata={"checks": (non_negative,)})
    fluid_change_rate: Union[float, None] = field(default=None, metadata={"checks": (non_negative,)})
    calcification_rate: Union[float, None] = field(default=None, metadata={"checks": (non_negative,)})

    @classmethod
    def from_dict(cls, data: Dict):
        cls._reject_unknown_keys(data)
        return cls(**data)

@dataclass()
class VolumeConfig(Validated):
    """
    Volume parameters of a phase. `None` means "use the :class:`CellVolumes` default".

    :param target_fluid_fraction: Fraction of the cell volume it will attempt to keep as fluid. In [0, 1]
    :param nuclear_solid_target: Nuclear solid volume the volume model will tend towards. `>= 0`
    :param cytoplasm_solid_target: Cytoplasm solid volume the volume model will tend towards. `>= 0`
    :param relative_rupture_volume: Proportion of the initial volume that causes cell lysis. `> 0`
    :param nuclear_fluid: Initial volume of the fluid part of the nucleus. `>= 0`
    :param nuclear_solid: Initial volume of the solid part of the nucleus. `>= 0`
    :param cytoplasm_fluid: Initial volume of the fluid part of the cytoplasm. `>= 0`
    :param cytoplasm_solid: Initial volume of the solid part of the cytoplasm. `>= 0`
    :param target_cytoplasm_to_nuclear_ratio: Cytoplasmic to nuclear volumes ratio the volume model will tend
        towards. `>= 0`
    :param calcified_fraction: Initial fraction of the cell that is calcified. In [0, 1]
    :param simulated_cell_volume: Volume of the simulated cell (e.g., a CompuCell3D or Tissue Forge cell). `> 0`
    :param rates: Volume change rates
    """
    target_fluid_fraction: Union[float, None] = field(default=None, metadata={"checks": (fraction,)})
    nuclear_solid_target: Union[float, None] = field(default=None, metadata={"checks": (non_negative,)})
    cytoplasm_solid_target: Union[float, None] = field(default=None, metadata={"checks": (non_negative,)})
    relative_rupture_volume: Union[float, None] = field(default=None, metadata={"checks": (positive,)})
    nuclear_fluid: Union[float, None] = field(default=None, metadata={"checks": (non_negative,)})
    nuclear_solid: Union[float, None] = field(default=None, metadata={"checks": (non_negative,)})
    cytoplasm_fluid: Union[float, None] = field(default=None, metadata={"checks": (non_negative,)})
    cytoplasm_solid: Union[float, None] = field(default=None, metadata={"checks": (non_negative,)})
    target_cytoplasm_to_nuclear_ratio: Union[float, None] = field(default=None, metadata={"checks": (non_negative,)})
    calcified_fraction: Union[float, None] = field(default=None, metadata={"checks": (fraction,)})
    simulated_cell_volume: Union[float, None] = field(default=None, metadata={"checks": (positive,)})
    rates: VolumeRatesConfig = field(default_factory=VolumeRatesConfig)

    def validate(self):
        if not isinstance(self.rates, VolumeRatesConfig):
            raise TypeError(f"'rates' must be a VolumeRatesConfig. Got {type(self.rates).__name__}.")

    @classmethod
    def from_dict(cls, data: Dict):
        cls._reject_unknown_keys(data)
        return cls(**{**data, "rates": VolumeRatesConfig.from_dict(data.get("rates", {}))})


@dataclass()
class TimingConfig(Validated):
    """
    Timing parameters of a phase. `None` means "use the phase class' default".

    :param phase_duration: Expected duration of the phase, in units of the phenotype's `time_unit`. In the case of
        the stochastic transition, the transition rate will be `dt/phase_duration`. `> 0`
    :param fixed_duration: Sets the transition from this phase to the next to be deterministic (True) or stochastic
        (False)
    """
    phase_duration: Union[float, None] = field(
        default=None, metadata={"checks": (positive,)})
    fixed_duration: Union[bool, None] = field(
        default=None, metadata={"checks": (boolean,)})

    @classmethod
    def from_dict(cls, data: Dict):
        cls._reject_unknown_keys(data)
        return cls(**data)


@dataclass()
class EventConfig(Validated):
    """
    Events at the exit of a phase. `None` means "use the phase class' default".

    :param division_at_phase_exit: If the simulated cell should divide when leaving this phase
    :param removal_at_phase_exit: If the simulated cell should be removed (i.e., it dies, is killed, leaves the
        simulated domain) when leaving this phase
    """
    division_at_phase_exit: Union[bool, None] = field(
        default=None, metadata={"checks": (boolean,)})
    removal_at_phase_exit: Union[bool, None] = field(
        default=None, metadata={"checks": (boolean,)})

    @classmethod
    def from_dict(cls, data: Dict):
        cls._reject_unknown_keys(data)
        return cls(**data)


@dataclass()
class FunctionsConfig(Validated):
    """
    User-defined functions of a phase and their arguments. All functions must be *args functions. For the entry,
    exit, arrest, and transition functions `None` means "use the phase class' default function" and `False` means "no
    function" (where the phase class supports it). Functions can't be stored in JSON, so they can only be set from
    Python.

    :param entry_function: Called immediately when entering this phase. Some phases have a default entry function,
        such as Apoptosis. Must have no return
    :param exit_function: Called just before exiting this phase. Some phases have a default exit function. Must have
        no return
    :param arrest_function: Returns if the cell should exit the current phase early and the cell cycle and enter
        senescence
    :param check_transition_to_next_phase_function: Returns if the cell should advance to the next phase in the
        phenotype. If None the phase uses a deterministic or stochastic transition depending on `fixed_duration`
    :param user_phase_time_step: User-defined function to be executed with the time-step
    :param *_args: Args (list or tuple) for the function of the same name. Required if the function is defined
    """
    entry_function: object = field(default=None, metadata={"checks": (function_or_false,)})
    entry_function_args: Union[List, Tuple, None] = field(default=None, metadata={"checks": (sequence,)})
    exit_function: object = field(default=None, metadata={"checks": (function_or_false,)})
    exit_function_args: Union[List, Tuple, None] = field(default=None, metadata={"checks": (sequence,)})
    arrest_function: object = field(default=None, metadata={"checks": (function_or_false,)})
    arrest_function_args: Union[List, Tuple, None] = field(default=None, metadata={"checks": (sequence,)})
    check_transition_to_next_phase_function: object = field(default=None, metadata={"checks": (function_or_false,)})
    check_transition_to_next_phase_function_args: Union[List, Tuple, None] = field(default=None,
                                                                             metadata={"checks": (sequence,)})
    user_phase_time_step: object = field(default=None, metadata={"checks": (function,)})
    user_phase_time_step_args: Union[List, Tuple, None] = field(default=None, metadata={"checks": (sequence,)})

    def validate(self):
        for name in ("entry_function", "exit_function", "arrest_function", "check_transition_to_next_phase_function",
                     "user_phase_time_step"):
            if callable(getattr(self, name)) and getattr(self, f"{name}_args") is None:
                raise ValueError(f"'{name}' is defined but '{name}_args' is not. Expected a list or tuple.")

    @classmethod
    def from_dict(cls, data: Dict):
        cls._reject_unknown_keys(data)
        return cls(**data)


@dataclass()
class PhaseConfig(Validated):
    """
    :param name: Descriptive name of the phase (e.g., S, G, M, necrotic swelling). `None` means "use the phase class'
        default"
    :param index: Position of the phase in the phenotype's phase list. `>= 0`
    :param kind: Name of the :mod:`phenocellpy.phases` class to build this phase with (e.g., "Ki67Positive"). `None`
        means "the phenotype's class for this position" (:class:`Phases.Phase` for a generic :class:`Phenotype`)
    :param previous_phase_index: Index of the phase preceding this phase in the phenotype
    :param next_phase_index: Index of the phase the cycle goes to when this phase is exited
    """
    name: Union[str, None] = field(default=None, metadata={"checks": (non_empty_str,)})
    index: int = field(default=0, metadata={"checks": (not_none, integer, non_negative)})
    kind: Union[str, None] = field(default=None, metadata={"checks": (non_empty_str,)})
    previous_phase_index: Union[int, None] = field(default=None, metadata={"checks": (integer,)})
    next_phase_index: Union[int, None] = field(default=None, metadata={"checks": (integer,)})

    timing: TimingConfig = field(default_factory=TimingConfig)
    volume: VolumeConfig = field(default_factory=VolumeConfig)
    events: EventConfig = field(default_factory=EventConfig)
    functions: FunctionsConfig = field(default_factory=FunctionsConfig)

    _groups = {"timing": TimingConfig, "volume": VolumeConfig, "events": EventConfig, "functions": FunctionsConfig}

    def validate(self):
        for name, expected in self._groups.items():
            value = getattr(self, name)
            if not isinstance(value, expected):
                raise TypeError(f"Phase {self._label()}: '{name}' must be a {expected.__name__}. "
                                f"Got {type(value).__name__}.")

    def _label(self):
        """Identifies the phase in error messages. The name is only filled from the phase class when it is built."""
        return f"'{self.name}'" if self.name is not None else f"with index {self.index}"

    @classmethod
    def from_dict(cls, data: Dict):
        cls._reject_unknown_keys(data)
        kwargs = {key: data[key] for key in ("name", "index", "kind", "previous_phase_index", "next_phase_index")
                  if key in data}
        kwargs.update({name: group.from_dict(data.get(name, {})) for name, group in cls._groups.items()})
        return cls(**kwargs)


@dataclass()
class PhenotypeConfig(Validated):
    """
    :param name: Name for the phenotype
    :type str
    :param dt: time-step duration in units of `time_unit`. `> 0`. Shared by all phases
    :type float
    :param time_unit: Time unit
    :type str
    :param space_unit: Space unit
    :type space_unit: str
    :param phases: The different phases of the phenotype. `phases[i].index` must be `i`
    :type tuple of PhaseConfig (a list is accepted and converted)
    :param senescent_phase: Special outside-of-phenotype-order senescent phase. If None the default
        :class:`Phases.SenescentPhase` is used. If False the phenotype won't have a senescent phase
    :type PhaseConfig, None, or False
    :param starting_phase_index: Which phase to start the phenotype model at. -1 (random) is currently disabled and
        falls back to 0
    :type int
    :param user_phenotype_time_step: User-defined function to be executed with the time-step
    :type function
    :param user_phenotype_time_step_args: args for `user_phenotype_time_step`
    :type list or tuple
    """
    name: str = field(default="unnamed", metadata={"checks": (not_none, non_empty_str)})
    dt: float = field(default=1.0, metadata={"checks": (not_none, positive)})
    time_unit: str = field(default="min", metadata={"checks": (not_none, non_empty_str)})
    space_unit: str = field(default="micrometer", metadata={"checks": (not_none, non_empty_str)})
    phases: Tuple[PhaseConfig, ...] = ()
    senescent_phase: Union[PhaseConfig, bool, None] = None
    starting_phase_index: Union[int, None] = field(default=0, metadata={"checks": (integer,)})
    user_phenotype_time_step: object = field(default=None, metadata={"checks": (function,)})
    user_phenotype_time_step_args: Union[List, Tuple, None] = field(default=(None,), metadata={"checks": (sequence,)})

    def validate(self):
        self.phases = tuple(self.phases)

        if not self.phases:
            raise ValueError(f"Phenotype '{self.name}' needs at least one phase.")

        for position, phase in enumerate(self.phases):
            if not isinstance(phase, PhaseConfig):
                raise TypeError(f"Phenotype '{self.name}': phases must be PhaseConfig. "
                                f"Got {type(phase).__name__} at position {position}.")
            # the phase index is used to look up the phase in the phenotype's phase list
            if phase.index != position:
                where = f"phase '{phase.name}' is at position {position} but has" if phase.name is not None else \
                    f"the phase at position {position} has"
                raise ValueError(f"Phenotype '{self.name}': {where} index {phase.index}. Phase indices must match their "
                                 f"position in `phases`.")
            # same lookup when the phase is exited
            if phase.next_phase_index is None or not 0 <= phase.next_phase_index < len(self.phases):
                raise ValueError(f"Phenotype '{self.name}': phase {phase._label()} has `next_phase_index` "
                                 f"{phase.next_phase_index}. It must be in [0, {len(self.phases) - 1}].")

        if isinstance(self.senescent_phase, PhaseConfig):
            if self.senescent_phase.index < len(self.phases):
                raise ValueError(f"Phenotype '{self.name}': senescent phase index {self.senescent_phase.index} "
                                 f"collides with the regular phase indices [0, {len(self.phases) - 1}].")
        elif self.senescent_phase is not None and self.senescent_phase is not False:
            raise ValueError(f"Phenotype '{self.name}': `senescent_phase` must be a PhaseConfig, False, or None. "
                             f"Got {self.senescent_phase!r}.")

        if self.starting_phase_index is not None and self.starting_phase_index != -1 and \
                not 0 <= self.starting_phase_index < len(self.phases):
            raise ValueError(f"Phenotype '{self.name}': `starting_phase_index` must be -1 or in "
                             f"[0, {len(self.phases) - 1}]. Got {self.starting_phase_index}.")

        if self.user_phenotype_time_step is not None and self.user_phenotype_time_step_args is None:
            raise ValueError(f"Phenotype '{self.name}': `user_phenotype_time_step` is defined but "
                             f"`user_phenotype_time_step_args` is not. Expected a list or tuple.")

    @classmethod
    def from_dict(cls, data: Dict):
        """
        Builds the config from a dict (e.g., loaded from JSON). Missing keys take the dataclass defaults.
        `"senescent_phase": false` disables the senescent phase; `null` or omitting it uses the default one.
        """
        cls._reject_unknown_keys(data)
        kwargs = {key: data[key] for key in ("name", "dt", "time_unit", "space_unit", "starting_phase_index",
                                             "user_phenotype_time_step", "user_phenotype_time_step_args")
                  if key in data}

        kwargs["phases"] = tuple(PhaseConfig.from_dict(phase_data) for phase_data in data.get("phases", ()))

        senescent = data.get("senescent_phase")
        kwargs["senescent_phase"] = PhaseConfig.from_dict(senescent) if isinstance(senescent, dict) else senescent

        return cls(**kwargs)

    @classmethod
    def from_json_file(cls, filename: str):
        with open(filename, "r") as f:
            data = json.load(f)

        return cls.from_dict(data)


if __name__ == "__main__":
    # run as `python -m phenocellpy.types.types <file.json>`. Running this file directly puts this directory on
    # sys.path and this module shadows the standard library's `types`.
    import sys

    print(PhenotypeConfig.from_json_file(sys.argv[1]))
