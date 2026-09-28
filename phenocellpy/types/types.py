import json
from dataclasses import dataclass, field

from phenocellpy.types.base import Validated
from phenocellpy.types.validators import (not_none, positive, non_negative, fraction, boolean, integer,
                                          non_empty_str)

# Configs are frozen: checks only run at construction, so allowing `cfg.dt = -1` afterwards would bypass them.
# Runtime state (time in phase, current volumes, ...) belongs to the Phase/Phenotype objects, not to these configs.
#
# Time-related values (`phase_duration`) are in the units of the owning PhenotypeConfig's `time_unit`. `dt` only
# lives in PhenotypeConfig; it is passed down to the phases when the runtime objects are built.


@dataclass(frozen=True)
class VolumeConfig(Validated):
    """
    Volume parameters of a phase. `None` means "use the :class:`CellVolumes` default".
    """
    target_fluid_fraction: float | None = field(default=None, metadata={"checks": (fraction,)})
    nuclear_solid_target: float | None = field(default=None, metadata={"checks": (non_negative,)})
    cytoplasm_solid_target: float | None = field(default=None, metadata={"checks": (non_negative,)})
    relative_rupture_volume: float | None = field(default=None, metadata={"checks": (positive,)})

    @classmethod
    def from_dict(cls, data: dict):
        cls._reject_unknown_keys(data)
        return cls(**data)


@dataclass(frozen=True)
class TimingConfig(Validated):
    """
    :param phase_duration: Expected duration of the phase, in units of the phenotype's `time_unit`. `> 0`
    :param fixed_duration: If True the phase lasts exactly `phase_duration`, otherwise the exit is stochastic
    """
    phase_duration: float = field(default=10.0, metadata={"checks": (not_none, positive)})
    fixed_duration: bool = field(default=False, metadata={"checks": (not_none, boolean)})

    @classmethod
    def from_dict(cls, data: dict):
        cls._reject_unknown_keys(data)
        return cls(**data)


@dataclass(frozen=True)
class EventConfig(Validated):
    division_at_phase_exit: bool = field(default=False, metadata={"checks": (not_none, boolean)})
    removal_at_phase_exit: bool = field(default=False, metadata={"checks": (not_none, boolean)})

    @classmethod
    def from_dict(cls, data: dict):
        cls._reject_unknown_keys(data)
        return cls(**data)


@dataclass(frozen=True)
class PhaseConfig(Validated):
    """
    :param name: Name of the phase
    :param index: Position of the phase in the phenotype's phase list. `>= 0`
    """
    name: str = field(default="unnamed", metadata={"checks": (not_none, non_empty_str)})
    index: int = field(default=0, metadata={"checks": (not_none, integer, non_negative)})

    timing: TimingConfig = field(default_factory=TimingConfig)
    volume: VolumeConfig = field(default_factory=VolumeConfig)
    events: EventConfig = field(default_factory=EventConfig)

    def validate(self):
        for name, expected in (("timing", TimingConfig), ("volume", VolumeConfig), ("events", EventConfig)):
            value = getattr(self, name)
            if not isinstance(value, expected):
                raise TypeError(f"Phase '{self.name}': '{name}' must be a {expected.__name__}. "
                                f"Got {type(value).__name__}.")

    @classmethod
    def from_dict(cls, data: dict):
        cls._reject_unknown_keys(data)
        kwargs = {key: data[key] for key in ("name", "index") if key in data}
        return cls(
            **kwargs,
            timing=TimingConfig.from_dict(data.get("timing", {})),
            volume=VolumeConfig.from_dict(data.get("volume", {})),
            events=EventConfig.from_dict(data.get("events", {})),
        )


def default_senescent_phase() -> PhaseConfig:
    """
    Config equivalent of the default :class:`Phases.SenescentPhase`. Its duration (60 days) assumes `time_unit` is
    minutes.
    """
    return PhaseConfig(
        name="senescent",
        index=9999,
        timing=TimingConfig(phase_duration=60 * 24 * 60, fixed_duration=True),
    )


@dataclass(frozen=True)
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
    :param senescent_phase: Special outside-of-phenotype-order senescent phase. None disables it
    :type PhaseConfig or None
    """
    name: str = field(default="unnamed", metadata={"checks": (not_none, non_empty_str)})
    dt: float = field(default=1.0, metadata={"checks": (not_none, positive)})
    time_unit: str = field(default="min", metadata={"checks": (not_none, non_empty_str)})
    space_unit: str = field(default="micrometer", metadata={"checks": (not_none, non_empty_str)})
    phases: tuple[PhaseConfig, ...] = ()
    senescent_phase: PhaseConfig | None = field(default_factory=default_senescent_phase)

    def validate(self):
        # frozen, so bypass __setattr__. A tuple keeps the phase list itself immutable too
        object.__setattr__(self, "phases", tuple(self.phases))

        if not self.phases:
            raise ValueError(f"Phenotype '{self.name}' needs at least one phase.")

        for position, phase in enumerate(self.phases):
            if not isinstance(phase, PhaseConfig):
                raise TypeError(f"Phenotype '{self.name}': phases must be PhaseConfig. "
                                f"Got {type(phase).__name__} at position {position}.")
            # the phase index is used to look up the phase in the phenotype's phase list
            if phase.index != position:
                raise ValueError(f"Phenotype '{self.name}': phase '{phase.name}' is at position {position} but has "
                                 f"index {phase.index}. Phase indices must match their position in `phases`.")

        if self.senescent_phase is not None:
            if not isinstance(self.senescent_phase, PhaseConfig):
                raise TypeError(f"Phenotype '{self.name}': `senescent_phase` must be a PhaseConfig or None. "
                                f"Got {type(self.senescent_phase).__name__}.")
            if self.senescent_phase.index < len(self.phases):
                raise ValueError(f"Phenotype '{self.name}': senescent phase index {self.senescent_phase.index} "
                                 f"collides with the regular phase indices [0, {len(self.phases) - 1}].")

    @classmethod
    def from_dict(cls, data: dict):
        """
        Builds the config from a dict (e.g., loaded from JSON). Missing keys take the dataclass defaults.
        `"senescent_phase": false` or `null` disables the senescent phase; omitting it uses the default one.
        """
        cls._reject_unknown_keys(data)
        kwargs = {key: data[key] for key in ("name", "dt", "time_unit", "space_unit") if key in data}

        kwargs["phases"] = tuple(PhaseConfig.from_dict(phase_data) for phase_data in data.get("phases", ()))

        if "senescent_phase" in data:
            senescent = data["senescent_phase"]
            kwargs["senescent_phase"] = None if senescent is None or senescent is False else PhaseConfig.from_dict(senescent)

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
