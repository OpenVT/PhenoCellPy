import json
from dataclasses import dataclass, field


@dataclass
class VolumeConfig:
    target_fluid_fraction: float | None = None
    nuclear_solid_target: float | None = None
    cytoplasm_solid_target: float | None = None
    relative_rupture_volume: float | None = None

    @classmethod
    def from_dict(cls, data: dict):
        return cls(**data)


@dataclass
class TimingConfig:
    phase_duration: float = 10.0
    fixed_duration: bool = False

    @classmethod
    def from_dict(cls, data: dict):
        return cls(**data)


@dataclass
class EventConfig:
    division_at_phase_exit: bool = False
    removal_at_phase_exit: bool = False

    @classmethod
    def from_dict(cls, data: dict):
        return cls(**data)


@dataclass
class PhaseConfig:
    name: str
    index: int

    timing: TimingConfig = field(default_factory=TimingConfig)
    volume: VolumeConfig = field(default_factory=VolumeConfig)
    events: EventConfig = field(default_factory=EventConfig)

    @classmethod
    def from_dict(cls, data: dict):
        return cls(
            name=data["name"],
            index=data["index"],
            timing=TimingConfig.from_dict(
                data.get("timing", {})
            ),
            volume=VolumeConfig.from_dict(
                data.get("volume", {})
            ),
            events=EventConfig.from_dict(
                data.get("events", {})
            ),
        )


@dataclass
class PhenotypeConfig:
    name: str
    dt: float
    time_unit: str = "min"
    space_unit: str = "micrometer"
    phases: list[PhaseConfig] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: dict):
        return cls(
            name=data["name"],
            dt=data["dt"],
            time_unit=data.get("time_unit", "min"),
            space_unit=data.get("space_unit", "micrometer"),
            phases=[
                PhaseConfig.from_dict(phase_data)
                for phase_data in data.get("phases", [])
            ],
        )

    @classmethod
    def from_json_file(cls, filename: str):
        with open(filename, "r") as f:
            data = json.load(f)

        return cls.from_dict(data)


# ---------------------------------------------------------
# Load the definition.
# Nothing is evaluated or simulated here.
# ---------------------------------------------------------

phenotype = PhenotypeConfig.from_json_file("ki67_basic.json")

print(phenotype)