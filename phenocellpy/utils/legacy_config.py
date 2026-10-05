"""
BSD 3-Clause License

Copyright (c) 2023, Juliano Ferrari Gianlupi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# Builders of configs from the arguments of the old (pre-config) `CellVolumes`, `Phase`, and `Phenotype`
# constructors. Each builder has the signature, and defaults, of the old constructor, and returns the config the new
# constructor takes. Phases and senescent phases that were passed as Phase objects are passed as the configs built by
# `build_phase_config_legacy`.

from phenocellpy.types import (PhenotypeConfig, PhaseConfig, TimingConfig, VolumeConfig, VolumeRatesConfig,
                               EventConfig, FunctionsConfig)
from typing import Union, List

def build_cell_volumes_config_legacy(
        target_fluid_fraction=None,
        nuclear_fluid=None,
        nuclear_solid=None,
        nuclear_solid_target=None,
        cytoplasm_fluid=None,
        cytoplasm_solid=None,
        cytoplasm_solid_target=None,
        target_cytoplasm_to_nuclear_ratio=None,
        calcified_fraction=None,
        relative_rupture_volume=None,
        time_unit="min",
        space_unit="micrometer",
) -> VolumeConfig:
    """
    Builds the config of the old `CellVolumes(...)` constructor call with the same arguments. `time_unit` and
    `space_unit` are ignored: they are set by the phase the volume belongs to.

    :rtype: VolumeConfig
    """
    return VolumeConfig(
        target_fluid_fraction=target_fluid_fraction,
        nuclear_fluid=nuclear_fluid,
        nuclear_solid=nuclear_solid,
        nuclear_solid_target=nuclear_solid_target,
        cytoplasm_fluid=cytoplasm_fluid,
        cytoplasm_solid=cytoplasm_solid,
        cytoplasm_solid_target=cytoplasm_solid_target,
        target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
        calcified_fraction=calcified_fraction,
        relative_rupture_volume=relative_rupture_volume,
    )


def build_phase_config_legacy(
        index: int = None,
        previous_phase_index: int = None,
        next_phase_index: int = None,
        dt: float = None,
        time_unit: str = "min",
        space_unit="micrometer",
        name: str = None,
        division_at_phase_exit: bool = False,
        removal_at_phase_exit: bool = False,
        fixed_duration: bool = False,
        phase_duration: float = 10,
        entry_function=None,
        entry_function_args: List = None,
        exit_function=None,
        exit_function_args: List = None,
        arrest_function=None,
        arrest_function_args: List = None,
        check_transition_to_next_phase_function=None,
        check_transition_to_next_phase_function_args: List = None,
        simulated_cell_volume: float = None,
        cytoplasm_volume_change_rate=None,
        nuclear_volume_change_rate=None,
        calcification_rate=None,
        target_fluid_fraction=None,
        nuclear_fluid=None,
        nuclear_solid=None,
        nuclear_solid_target=None,
        cytoplasm_fluid=None,
        cytoplasm_solid=None,
        cytoplasm_solid_target=None,
        target_cytoplasm_to_nuclear_ratio=None,
        calcified_fraction=None,
        fluid_change_rate=None,
        relative_rupture_volume=None,
        user_phase_time_step=None,
        user_phase_time_step_args=(None,),
) -> PhaseConfig:
    """
    Builds the config of the old `Phases.Phase(...)` constructor call with the same arguments. `dt`, `time_unit`,
    and `space_unit` are ignored: they are set by the phenotype the phase belongs to.

    :rtype: PhaseConfig
    """
    return PhaseConfig(
        # the old Phase constructor replaced a None name and index with these
        name="unnamed" if name is None else name,
        index=0 if index is None else index,
        previous_phase_index=previous_phase_index,
        next_phase_index=next_phase_index,
        timing=TimingConfig(
            phase_duration=phase_duration,
            fixed_duration=fixed_duration,
        ),
        volume=VolumeConfig(
            target_fluid_fraction=target_fluid_fraction,
            nuclear_fluid=nuclear_fluid,
            nuclear_solid=nuclear_solid,
            nuclear_solid_target=nuclear_solid_target,
            cytoplasm_fluid=cytoplasm_fluid,
            cytoplasm_solid=cytoplasm_solid,
            cytoplasm_solid_target=cytoplasm_solid_target,
            target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio,
            calcified_fraction=calcified_fraction,
            relative_rupture_volume=relative_rupture_volume,
            simulated_cell_volume=simulated_cell_volume,
            rates=VolumeRatesConfig(
                cytoplasm_volume_change_rate=cytoplasm_volume_change_rate,
                nuclear_volume_change_rate=nuclear_volume_change_rate,
                fluid_change_rate=fluid_change_rate,
                calcification_rate=calcification_rate,
            ),
        ),
        events=EventConfig(
            division_at_phase_exit=division_at_phase_exit,
            removal_at_phase_exit=removal_at_phase_exit,
        ),
        functions=FunctionsConfig(
            entry_function=entry_function,
            entry_function_args=entry_function_args,
            exit_function=exit_function,
            exit_function_args=exit_function_args,
            arrest_function=arrest_function,
            arrest_function_args=arrest_function_args,
            check_transition_to_next_phase_function=check_transition_to_next_phase_function,
            check_transition_to_next_phase_function_args=check_transition_to_next_phase_function_args,
            user_phase_time_step=user_phase_time_step,
            user_phase_time_step_args=user_phase_time_step_args,
        ),
    )


def build_phenotype_config_legacy(
        name: str = "unnamed",
        dt: float = 1,
        time_unit: str = "min",
        space_unit="micrometer",
        phases: List = None,
        senescent_phase: Union[PhaseConfig, bool, None] = None,
        starting_phase_index: int = 0,
        user_phenotype_time_step=None,
        user_phenotype_time_step_args=(None,),
) -> PhenotypeConfig:
    """
    Builds the config of the old `Phenotype(...)` constructor call with the same arguments. `phases` and
    `senescent_phase` take configs built with the phase builders of this module (e.g.,
    :func:`build_phase_config_legacy`) instead of Phase objects.

    :rtype: PhenotypeConfig
    """
    if phases is None:
        phases = [build_phase_config_legacy(previous_phase_index=0, next_phase_index=0, dt=dt, time_unit=time_unit,
                                            space_unit=space_unit)]
    return PhenotypeConfig(
        name=name,
        dt=dt,
        time_unit=time_unit,
        space_unit=space_unit,
        phases=tuple(phases),
        senescent_phase=senescent_phase,
        starting_phase_index=starting_phase_index,
        user_phenotype_time_step=user_phenotype_time_step,
        user_phenotype_time_step_args=user_phenotype_time_step_args,
    )
