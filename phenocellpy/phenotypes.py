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

from warnings import warn

import phenocellpy.phases as Phases

from copy import deepcopy

from phenocellpy.utils.argcheck import _check_arguments


# from numpy.random import randint


# todo:
#  - a JJ tyson cycle model (https://www.ebi.ac.uk/biomodels/BIOMD0000000003)
#  - add biomodels ontology anotation
#  - switch how defalts are handled to be like phases.NecrosisLysed
#  - add cell heterogeneity. have (yet another) argument to set randomization to True/False (and possibly which distri-
#    tion to use, log-normal, normal, etc). If true, rates, durations, and volumes are slightly shuffled. Also have
#    an argument to set the variance
#  - interface class
#  - have the time unit define some unit conversions
#  - have some pre-built secretions/absorption and have it drive phenotype changes
#  - pre-calculate the transition probability when using the stochastic transition (no need to calculate it every step,
#    as it is fixed)


class Phenotype:
    """

    Base class to define a cell phenotype.

    Defines a cell phenotype, a sequence of phases with different behaviors. E.g., a quiescent-proliferating cell cycle
    is a phenotype with two phases (quiescence, and growth/division); the necrotic phenotype starts with a osmotic swe-
    ling phase, followed by dissolution of the cell into its media after it bursts.

    This class has methods to time-step the phenotype model (which time-steps all submodels), to change the phenotype
    phase to an arbitrary phase of the phenotype cycle, to go to the next phase in the cycle, and to go to a
    non-changing senescent phase.

    Methods:
    --------

    time_step_phenotype()
        Time-steps the phenotype model. Returns a tuple (cell changed phases, cell died, cell divided). See
        :func:`time_step_cycle` for further information.

    go_to_next_phase()
        Goes to the next phase in the cycle

    set_phase(index)
        Sets the current cycle's phase to be phase of index `index`

    go_to_quiescence()
        Moves cycle to quiescent phase

    user_phenotype_time_step(*args)
        User-defined function to be executed with the time-step

    Attributes
    ----------

    name : str
        Name of the phenotype

    dt : float
        Time-step size (in units of `time_unit`). Must be >0.

    time_unit : str
        Time unit. TODO: Defines time conversions

    phases : list
        Ordered list of phases this cycle goes through. Must be a list of :class:`Phases.Phase` objects.

    starting_phase_index : int
        Index of which phase to start at. Currently there is no option to start at a random phase, but that capability
        will be implemented soon. When it is, to start at a random phase set to -1

    senescent_phase : :class:`Phases.Phase`, None, or False
        If false the cycle won't have a (stand-alone) senescent phase defined. If None the default
        :class:`Phases.SenescentPhase` will be used as the stand-alone senescent phase. If it is a :class:`Phases.Phase`
        object it will be used as the stand-alone senescent phase.

    current_phase : :class:`Phases.Phase`
        The current (active) phase of the cycle.

    time_in_phenotype : float
        Total time elapsed for the cycle

    user_phenotype_time_step_args: list or tuple
        args for `user_phenotype_time_step`
    """

    def __init__(self, name: str = "unnamed", 
                 dt: float = 1, 
                 time_unit: str = "min", 
                 space_unit="micrometer",
                 phases: list = None, 
                 senescent_phase: Phases.Phase or False = None, 
                 starting_phase_index: int = 0,
                 user_phenotype_time_step=None, 
                 user_phenotype_time_step_args=(None,)):
        """
        :param name: Name for the phenotype
        :type str
        :param dt: time-step duration in units of `time_unit`
        :type float
        :param time_unit: Time unit
        :type str
        :param space_unit: Space unit
        :type space_unit: str
        :param phases: The different phases of the phenotype
        :type list of Phases.Phase
        :param senescent_phase: Special outside-of-phenotype-order senescent phase
        :type Phases.Phase or False
        :param starting_phase_index: Which phase to start the phenotype model at
        :type int
        :param user_phenotype_time_step: User-defined function to be executed with the time-step
        :type function
        :param user_phenotype_time_step_args: args for `user_phenotype_time_step`
        :type list or tuple
        """
        # todo: add __init__ parameters for custom functions for each class
        # todo: add alias for self.current_phase.volume, i.e. property self.volume that fetches
        #  self.current_phase.volume. If read-only easy to do, not sure how to do it if I want to keep write abilities
        self.name = name

        self.time_unit = time_unit
        self.space_unit = space_unit

        if dt <= 0 or dt is None:
            raise ValueError(f"'dt' must be greater than 0. Got {dt}.")
        self.dt = dt
        if phases is None:
            self.phases = [Phases.Phase(previous_phase_index=0, next_phase_index=0, dt=self.dt, time_unit=time_unit,
                                        space_unit=space_unit)]
        else:
            self.phases = phases
        if senescent_phase is None:
            self.senescent_phase = Phases.SenescentPhase(dt=self.dt)
        elif senescent_phase is not None and not senescent_phase:
            self.senescent_phase = False
        elif not isinstance(senescent_phase, Phases.Phase):
            raise ValueError(f"`senescent_phase` must Phases.Phase object, False, or None. Got {senescent_phase}")
        else:
            self.senescent_phase = senescent_phase
        if starting_phase_index is None:
            starting_phase_index = 0
        elif starting_phase_index == -1:  # random option
            # todo: fix this. it won't actually work for many cells, as this randomization happens at class init, but
            #  the same class object is then copied to the cells
            warn("Randomization of the initial phase is currently disabled. Setting the initial phase to be "
                 "phase of index 0.")
            starting_phase_index = 0
            # starting_phase_index = randint(0, len(self.phases) + 1)

        self.user_phenotype_time_step = user_phenotype_time_step
        if self.user_phenotype_time_step is not None:
            if not (type(user_phenotype_time_step_args) == list or type(user_phenotype_time_step_args) == tuple):
                raise ValueError(
                    f"`user_phenotype_time_step` is defined but `user_pheno_time_step_args` is not list or "
                    f"tuple.\nGot {type(user_phenotype_time_step_args)} instead")
            self.user_pheno_time_step_args = user_phenotype_time_step_args

        self.current_phase = self.phases[starting_phase_index]
        self.time_in_phenotype = 0

    def time_step_phenotype(self):
        """
        Time-steps the phenotype.

        Increments :attr:`time_in_cycle` by :attr:`dt`. Calls :func:`current_phase.time_step_phase`. If the phase time-
        step determines the cycle moves to the next phase (i.e., returns `True` for `next_phase`), calls
        :func:go_to_next_phase. If the phase time-step determines the cell exits the cell cycle and goes to senescence
        (i.e., returns `True` for `senes`) calls :func:`go_to_senescence`. If :attr:`time_in_cycle` is 0 and
        :attr:current_phase has an entry function calls :func:current_phase.entry_function.

        :return: Flags (bool) for phase changing, cell death, and cell division
        :rtype: tuple of bool
        """
        # if we initialize a phenotype that has a initial phase with an entry function we should execute it immediately
        if not self.time_in_phenotype and self.current_phase.entry_function is not None:
            self.current_phase.entry_function(*self.current_phase.entry_function_args)

        if self.user_phenotype_time_step is not None:
            self.user_phenotype_time_step(*self.user_pheno_time_step_args)

        self.time_in_phenotype += self.dt

        go_next_phase, exit_phenotype, transition_to_index = self.current_phase.time_step_phase()

        if go_next_phase:
            if transition_to_index is not None:
                phase_idx = self.current_phase.index
                old_next_phase_idx = self.current_phase.next_phase_index
                self.current_phase.next_phase_index = transition_to_index
            else:
                phase_idx = None
            changed_phases, cell_removed, cell_divides = self.go_to_next_phase()
            if phase_idx is not None:
                self.phases[phase_idx].next_phase_index = old_next_phase_idx
            return changed_phases, cell_removed, cell_divides
        elif exit_phenotype:
            self.go_to_senescence()
            changed_phases, cell_removed, cell_divides = (True, False, False)
            return changed_phases, cell_removed, cell_divides
        changed_phases, cell_removed, cell_divides = (False, False, False)
        return changed_phases, cell_removed, cell_divides

    def go_to_next_phase(self):
        """

        Sets the cycle phase to be the phase of index :attr:`current_phase.next_phase_index`.

        Gets the flags for division and death upon phase exit. Sets the cycle phase to be
        :attr:`current_phase.next_phase_index` by calling :func:`set_phase`. Calls :func:`current_phase.entry_function`
        (after phase change). Returns that phase change has occurred, and the division and death upon phase exit flags.

        :return: Flags (bool) for phase changing, cell death, and cell division
        :rtype: tuple of bool
        """
        changed_phases = True
        divides = self.current_phase.division_at_phase_exit
        removal = self.current_phase.removal_at_phase_exit
        self.set_phase(self.current_phase.next_phase_index)

        return changed_phases, removal, divides

    def set_phase(self, idx):
        """

        Sets cycle phase to be phase of index :param:`idx`.

        This function moves the cycle to an arbitrary phase and coordinates their volume attributes. Saves current
        :attr:`current_phase.volume` and :attr:`current_phase.volume` to variables, sets :attr:`current_phase` to be
        `phases[idx]`, resets :attr:`current_phase.volume` and :attr:`current_phase.volume` to be the previously saved
        values, sets :attr:`current_phase.time_in_phase` to 0.

        :param idx: index of list :attr:`phases`, which phase to go to.
        :return: No return
        """

        # get the current cytoplasm, nuclear, calcified volumes
        cyto_solid = self.current_phase.volume.cytoplasm_solid
        cyto_fluid = self.current_phase.volume.cytoplasm_fluid

        nucl_solid = self.current_phase.volume.nuclear_solid
        nucl_fluid = self.current_phase.volume.nuclear_fluid

        calc_frac = self.current_phase.volume.calcified_fraction

        # get the target volumes

        target_cytoplasm_solid = self.current_phase.volume.cytoplasm_solid_target
        # target_cyto_fluid_frac = self.current_phase.volume.target_cytoplasm_fluid_fraction

        target_nuclear_solid = self.current_phase.volume.nuclear_solid_target
        # target_nucl_fluid_frac = self.current_phase.volume.target_nuclear_fluid_fraction

        target_fluid_fraction = self.current_phase.volume.target_fluid_fraction

        # set parameters of next phase

        self.phases[idx].volume.cytoplasm_solid = cyto_solid
        self.phases[idx].volume.cytoplasm_fluid = cyto_fluid

        self.phases[idx].volume.nuclear_solid = nucl_solid
        self.phases[idx].volume.nuclear_fluid = nucl_fluid

        self.phases[idx].volume.calcified_fraction = calc_frac

        self.phases[idx].volume.cytoplasm_solid_target = target_cytoplasm_solid

        self.phases[idx].volume.nuclear_solid_target = target_nuclear_solid

        self.phases[idx].volume.target_fluid_fraction = target_fluid_fraction

        # set phase

        self.current_phase = self.phases[idx]
        self.current_phase.time_in_phase = 0

        if self.current_phase.entry_function is not None:
            self.current_phase.entry_function(*self.current_phase.entry_function_args)

    def go_to_senescence(self):
        """

        Sets cycle phase to be the senescent phase.

        This function checks that :attr:`senescent_phase` is a :class:`Phases.Phase`, if that's the case it sets
        :attr:`current_phase` to be :attr:`senescent_phase`. It also sets (after phase change) the
        :attr:`current_phase.volume` and :attr:`current_phase.volume` to be the previous phase :attr:`volume` by first
        saving it to a temporary variable. It resets :attr:`current_phase.time_in_phase` to be 0.

        :return: No return
        """
        if not isinstance(self.senescent_phase, Phases.Phase):
            return

        # get the current cytoplasm, nuclear, calcified volumes
        cyto_solid = self.current_phase.volume.cytoplasm_solid
        cyto_fluid = self.current_phase.volume.cytoplasm_fluid

        nucl_solid = self.current_phase.volume.nuclear_solid
        nucl_fluid = self.current_phase.volume.nuclear_fluid

        calc_frac = self.current_phase.volume.calcified_fraction

        # get the target volumes

        # target_cytoplasm_solid = self.current_phase.volume.cytoplasm_solid_target
        # # target_cyto_fluid_frac = self.current_phase.volume.target_cytoplasm_fluid_fraction
        #
        # target_nuclear_solid = self.current_phase.volume.nuclear_solid_target
        # # target_nucl_fluid_frac = self.current_phase.volume.target_nuclear_fluid_fraction
        #
        # target_fluid_fraction = self.current_phase.volume.target_fluid_fraction

        # setting the senescent phase volume parameters. As the cell is now senescent it shouldn't want to change its
        # volume, so we set the targets to be the current measurements
        self.senescent_phase.volume.cytoplasm_solid = cyto_solid
        self.senescent_phase.volume.cytoplasm_fluid = cyto_fluid

        self.senescent_phase.volume.nuclear_solid = nucl_solid
        self.senescent_phase.volume.nuclear_fluid = nucl_fluid

        self.senescent_phase.volume.nuclear_solid_target = nucl_solid
        self.senescent_phase.volume.cytoplasm_solid_target = cyto_solid

        self.senescent_phase.volume.calcified_fraction = calc_frac

        self.senescent_phase.volume.target_fluid_fraction = (cyto_fluid + nucl_fluid) / (nucl_solid + nucl_fluid +
                                                                                         cyto_fluid + cyto_solid)

        # set the phase to be senescent
        self.current_phase = self.senescent_phase
        self.current_phase.time_in_phase = 0

    def copy(self):
        return deepcopy(self)

    def __str__(self):
        phases = ""
        for p in self.phases:
            phases += f"{p._short_str}, "
        if len(phases) > 2:
            phases = phases[:-2]
        return f"{self.name} cycle, phases: {phases}, at memory {self.__repr__().split(' ')[-1][:-1]}"


class SimpleLiveCycle(Phenotype):
    """
    Inherits :class:`Phenotype`. Simplest alive cycle, it has only one phase. When progressing to the "next" phase the
    cell should divide.
    """

    def __init__(self, name="Simple Live", dt=0.1, time_unit="min", space_unit="micrometer", senescent_phase=False,
                 division_at_phase_exits=(True,), removal_at_phase_exits=(False,),
                 fixed_durations=(False,), phase_durations: list = (60 / 0.0432,),
                 entry_functions=(None,), entry_functions_args=(None,), exit_functions=(None,),
                 exit_functions_args=(None, ), arrest_functions=(None, ), arrest_functions_args=(None, ),
                 check_transition_to_next_phase_functions=(None,),
                 check_transition_to_next_phase_functions_args: list = (None,),
                 simulated_cell_volume=None, cytoplasm_volume_change_rate=(None,),
                 nuclear_volume_change_rate=(None, ), calcification_rate=(None,), calcified_fraction=(0,),
                 target_fluid_fraction=(None, ), nuclear_fluid=(None, ), nuclear_solid=(None, ),
                 nuclear_solid_target=(None, ), cytoplasm_fluid=(None,), cytoplasm_solid=(None,),
                 cytoplasm_solid_target=(None, ), target_cytoplasm_to_nuclear_ratio=(None, ),
                 fluid_change_rate=(None, ),
                 user_phenotype_time_step=None, user_phenotype_time_step_args=None, user_phases_time_step=None,
                 user_phases_time_step_args=None):
    # def __init__(self, time_unit: str = "min", space_unit="micrometer", name: str = "Simple Live", dt=1,
    #              user_phenotype_time_step=None, user_phenotype_time_step_args=None, user_phases_time_step=None,
    #              user_phases_time_step_args=None, phase_durations=[60 / 0.0432], fixed_durations=[None],
    #              cytoplasm_volume_change_rate=(None,)):
        if user_phases_time_step is None:
            user_phases_time_step = [None]
            user_phases_time_step_args = [None]

        phases = [
            Phases.Phase(index=0, previous_phase_index=0, next_phase_index=0, dt=dt, time_unit=time_unit,
                         space_unit=space_unit, name="alive",
                         division_at_phase_exit=division_at_phase_exits[0], phase_duration=phase_durations[0],
                         user_phase_time_step=user_phases_time_step[0],
                         user_phase_time_step_args=user_phases_time_step_args[0], fixed_duration=fixed_durations[0],
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[0],
                         removal_at_phase_exit=removal_at_phase_exits, entry_function=entry_functions[0],
                         entry_function_args=entry_functions_args[0], exit_function=exit_functions[0],
                         exit_function_args=exit_functions_args[0], arrest_function=arrest_functions[0],
                         arrest_function_args=arrest_functions_args[0],
                         check_transition_to_next_phase_function=check_transition_to_next_phase_functions[0],
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_functions_args[0],
                         simulated_cell_volume=simulated_cell_volume,
                         nuclear_volume_change_rate=nuclear_volume_change_rate[0],
                         calcification_rate=calcification_rate[0],
                         calcified_fraction=calcified_fraction[0], target_fluid_fraction=target_fluid_fraction[0],
                         nuclear_fluid=nuclear_fluid[0], nuclear_solid=nuclear_solid[0],
                         nuclear_solid_target=nuclear_solid_target[0], cytoplasm_fluid=cytoplasm_fluid[0],
                         cytoplasm_solid=cytoplasm_solid[0], cytoplasm_solid_target=cytoplasm_solid_target[0],
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                         fluid_change_rate=fluid_change_rate[0])]
        super().__init__(name=name, dt=dt, time_unit=time_unit, space_unit=space_unit, phases=phases,
                         senescent_phase=False, user_phenotype_time_step=user_phenotype_time_step,
                         user_phenotype_time_step_args=user_phenotype_time_step_args)


class Ki67Basic(Phenotype):
    """

    Inherits :class:`Phenotype`. Simple proliferating-senescent phase. Cell divides upon leaving Ki67+

    This is a two phase cycle. Ki67- (defined in :class:`Phases.Ki67Negative`) is the senescent phase, Ki67-'s mean du-
    ration is 4.59h (stochastic transition to Ki67+ by default). Ki67+ (defined in :class:`Phases.Ki67Positive`) is the
    proliferating phase. It is responsible for doubling the volume of the cell at a rate of
    [increase in volume]/[Ki67+ duration]. Ki67+ duration is fixed (by default) at 15.5 hours. Once the cell exits Ki67+
    it divides.

    """

    def __init__(self, name="Ki67 Basic", dt=0.1, time_unit="min", space_unit="micrometer", senescent_phase=False,
                 division_at_phase_exits=(False, True), removal_at_phase_exits=(False, False),
                 fixed_durations=(False, True), phase_durations: list = (4.59 * 60, 15.5 * 60.0),
                 entry_functions=(None, None), entry_functions_args=(None, None), exit_functions=(None, None),
                 exit_functions_args=(None, None), arrest_functions=(None, None), arrest_functions_args=(None, None),
                 check_transition_to_next_phase_functions=(None, None),
                 check_transition_to_next_phase_functions_args: list = (None, None),
                 simulated_cell_volume=None, cytoplasm_volume_change_rate=(None, None),
                 nuclear_volume_change_rate=(None, None), calcification_rate=(None, None), calcified_fraction=(0, 0),
                 target_fluid_fraction=(None, None), nuclear_fluid=(None, None), nuclear_solid=(None, None),
                 nuclear_solid_target=(None, None), cytoplasm_fluid=(None, None), cytoplasm_solid=(None, None),
                 cytoplasm_solid_target=(None, None), target_cytoplasm_to_nuclear_ratio=(None, None),
                 fluid_change_rate=(None, None),
                 user_phenotype_time_step=None, user_phenotype_time_step_args=None, user_phases_time_step=None,
                 user_phases_time_step_args=None):
        if user_phases_time_step is None:
            user_phases_time_step = [None, None]
            user_phases_time_step_args = [None, None]
        _check_arguments(2, name, division_at_phase_exits, removal_at_phase_exits, fixed_durations, phase_durations,
                         entry_functions, entry_functions_args, exit_functions, exit_functions_args, arrest_functions,
                         arrest_functions_args, check_transition_to_next_phase_functions,
                         check_transition_to_next_phase_functions_args,
                         cytoplasm_volume_change_rate, nuclear_volume_change_rate, calcification_rate,
                         calcified_fraction, target_fluid_fraction, nuclear_fluid, nuclear_solid, nuclear_solid_target,
                         cytoplasm_fluid, cytoplasm_solid, cytoplasm_solid_target, target_cytoplasm_to_nuclear_ratio,
                         fluid_change_rate)

        Ki67_positive = Phases.Ki67Positive(index=1, previous_phase_index=0, next_phase_index=0, dt=dt,
                                            time_unit=time_unit, space_unit=space_unit,
                                            division_at_phase_exit=division_at_phase_exits[1],
                                            removal_at_phase_exit=removal_at_phase_exits[1],
                                            fixed_duration=fixed_durations[1], phase_duration=phase_durations[1],
                                            entry_function=entry_functions[1],
                                            entry_function_args=entry_functions_args[1],
                                            exit_function=exit_functions[1], exit_function_args=exit_functions_args[1],
                                            arrest_function=arrest_functions[1],
                                            arrest_function_args=arrest_functions_args[1],
                                            check_transition_to_next_phase_function=
                                            check_transition_to_next_phase_functions[1],
                                            check_transition_to_next_phase_function_args=
                                            check_transition_to_next_phase_functions_args[1],
                                            simulated_cell_volume=simulated_cell_volume,
                                            cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[1],
                                            nuclear_volume_change_rate=nuclear_volume_change_rate[1],
                                            calcification_rate=calcification_rate[1],
                                            target_fluid_fraction=target_fluid_fraction[1],
                                            nuclear_fluid=nuclear_fluid[1], nuclear_solid=nuclear_solid[1],
                                            nuclear_solid_target=nuclear_solid_target[1],
                                            cytoplasm_fluid=cytoplasm_fluid[1], cytoplasm_solid=cytoplasm_solid[1],
                                            cytoplasm_solid_target=cytoplasm_solid_target[1],
                                            target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[1],
                                            calcified_fraction=calcified_fraction[1],
                                            fluid_change_rate=fluid_change_rate[1],
                                            user_phase_time_step=user_phases_time_step[1],
                                            user_phase_time_step_args=user_phases_time_step_args[1])

        Ki67_negative = Phases.Ki67Negative(index=0, previous_phase_index=1, next_phase_index=1, dt=dt,
                                            time_unit=time_unit, space_unit=space_unit,
                                            division_at_phase_exit=division_at_phase_exits[0],
                                            removal_at_phase_exit=removal_at_phase_exits[0],
                                            fixed_duration=fixed_durations[0], phase_duration=phase_durations[0],
                                            entry_function=entry_functions[0],
                                            entry_function_args=entry_functions_args[0],
                                            exit_function=exit_functions[0], exit_function_args=exit_functions_args[0],
                                            arrest_function=arrest_functions[0],
                                            arrest_function_args=arrest_functions_args[0],
                                            check_transition_to_next_phase_function=
                                            check_transition_to_next_phase_functions[0],
                                            check_transition_to_next_phase_function_args=
                                            check_transition_to_next_phase_functions_args[0],
                                            simulated_cell_volume=simulated_cell_volume,
                                            cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[0],
                                            nuclear_volume_change_rate=nuclear_volume_change_rate[0],
                                            calcification_rate=calcification_rate[0],
                                            target_fluid_fraction=target_fluid_fraction[0],
                                            nuclear_fluid=nuclear_fluid[0], nuclear_solid=nuclear_solid[0],
                                            nuclear_solid_target=nuclear_solid_target[0],
                                            cytoplasm_fluid=cytoplasm_fluid[0], cytoplasm_solid=cytoplasm_solid[0],
                                            cytoplasm_solid_target=cytoplasm_solid_target[0],
                                            target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                                            calcified_fraction=calcified_fraction[0],
                                            fluid_change_rate=fluid_change_rate[0],
                                            user_phase_time_step=user_phases_time_step[0],
                                            user_phase_time_step_args=user_phases_time_step_args[0])

        phases = [Ki67_negative, Ki67_positive]

        super().__init__(name=name, dt=dt, time_unit=time_unit, space_unit=space_unit, phases=phases,
                         senescent_phase=senescent_phase, user_phenotype_time_step=user_phenotype_time_step,
                         user_phenotype_time_step_args=user_phenotype_time_step_args)


class Ki67Advanced(Phenotype):
    """

    Inherits :class:`Phenotype`. Simple senescent-proliferating (mitosis)-rest cycle. Cell divides upon leaving Ki67+ pre

    This is a three phase cycle. Ki67- (defined in :class:`Phases.Ki67Negative`) is the senescent phase, Ki67-'s mean
    duration is 3.62h (stochastic transition to Ki67+ pre-mitotic by default). Ki67+ pre-mitotic (defined in
    :class:`Phases.Ki67PositivePreMitotic`) is the proliferating phase. It is responsible for doubling the volume of the
    cell at a rate of [increase in volume]/[Ki67+ pre duration]. Ki67+ pre-mitotic duration is fixed (by default) at 13
    hours. Once the cell exits Ki67+ pre-mitotic it divides and enters Ki67+ post-mitotic. Ki67+ post-mitotic (defined
    in :class:`Phases.Ki67PositivePostMitotic`) is a rest phase, Ki67+ post-mitotic duration is fixed (by default)
    at 2.5 hours. Afterwards the cycle loops back to :class:`Phases.Ki67Negative`.

    """

    def __init__(self, name="Ki67 Advanced", dt=0.1, time_unit="min", space_unit="micrometer", senescent_phase=False,
                 division_at_phase_exits=(False, True, False), removal_at_phase_exits=(False, False, False),
                 fixed_durations=(False, True, True), phase_durations: list = (3.62 * 60, 13.0 * 60.0, 2.5 * 60),
                 entry_functions=(None, None, None), entry_functions_args=(None, None, None),
                 exit_functions=(None, False, None), exit_functions_args=(None, None, None),
                 arrest_functions=(None, None, None), arrest_functions_args=(None, None, None),
                 check_transition_to_next_phase_functions=(None, None, None),
                 check_transition_to_next_phase_functions_args: list = (None, None, None), simulated_cell_volume=None,
                 cytoplasm_volume_change_rate=(None, None, None),
                 nuclear_volume_change_rate=(None, None, None), calcification_rate=(None, None, None),
                 calcified_fraction=(0, 0, 0),
                 target_fluid_fraction=(None, None, None), nuclear_fluid=(None, None, None),
                 nuclear_solid=(None, None, None),
                 nuclear_solid_target=(None, None, None), cytoplasm_fluid=(None, None, None),
                 cytoplasm_solid=(None, None, None),
                 cytoplasm_solid_target=(None, None, None), target_cytoplasm_to_nuclear_ratio=(None, None, None),
                 fluid_change_rate=(None, None, None),
                 user_phenotype_time_step=None, user_phenotype_time_step_args=None, user_phases_time_step=None,
                 user_phases_time_step_args=None):
        if user_phases_time_step is None:
            user_phases_time_step = [None, None, None]
            user_phases_time_step_args = [None, None, None]
        _check_arguments(3, name, division_at_phase_exits, removal_at_phase_exits, fixed_durations, phase_durations,
                         entry_functions, entry_functions_args, exit_functions, exit_functions_args, arrest_functions,
                         arrest_functions_args, check_transition_to_next_phase_functions,
                         check_transition_to_next_phase_functions_args,
                         cytoplasm_volume_change_rate, nuclear_volume_change_rate, calcification_rate,
                         calcified_fraction, target_fluid_fraction, nuclear_fluid, nuclear_solid, nuclear_solid_target,
                         cytoplasm_fluid, cytoplasm_solid, cytoplasm_solid_target, target_cytoplasm_to_nuclear_ratio,
                         fluid_change_rate)

        Ki67_negative = Phases.Ki67Negative(index=0, previous_phase_index=2, next_phase_index=1, dt=dt,
                                            time_unit=time_unit, space_unit=space_unit,
                                            division_at_phase_exit=division_at_phase_exits[0],
                                            removal_at_phase_exit=removal_at_phase_exits[0],
                                            fixed_duration=fixed_durations[0], phase_duration=phase_durations[0],
                                            entry_function=entry_functions[0],
                                            entry_function_args=entry_functions_args[0],
                                            exit_function=exit_functions[0], exit_function_args=exit_functions_args[0],
                                            arrest_function=arrest_functions[0],
                                            arrest_function_args=arrest_functions_args[0],
                                            check_transition_to_next_phase_function=
                                            check_transition_to_next_phase_functions[0],
                                            check_transition_to_next_phase_function_args=
                                            check_transition_to_next_phase_functions_args[0],
                                            simulated_cell_volume=simulated_cell_volume,
                                            cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[0],
                                            nuclear_volume_change_rate=nuclear_volume_change_rate[0],
                                            calcification_rate=calcification_rate[0],
                                            target_fluid_fraction=target_fluid_fraction[0],
                                            nuclear_fluid=nuclear_fluid[0], nuclear_solid=nuclear_solid[0],
                                            nuclear_solid_target=nuclear_solid_target[0],
                                            cytoplasm_fluid=cytoplasm_fluid[0], cytoplasm_solid=cytoplasm_solid[0],
                                            cytoplasm_solid_target=cytoplasm_solid_target[0],
                                            target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                                            calcified_fraction=calcified_fraction[0],
                                            fluid_change_rate=fluid_change_rate[0],
                                            user_phase_time_step=user_phases_time_step[0],
                                            user_phase_time_step_args=user_phases_time_step_args[0])

        Ki67_positive_pre = Phases.Ki67PositivePreMitotic(index=1, previous_phase_index=0, next_phase_index=2, dt=dt,
                                                          time_unit=time_unit, space_unit=space_unit,
                                                          division_at_phase_exit=division_at_phase_exits[1],
                                                          removal_at_phase_exit=removal_at_phase_exits[1],
                                                          fixed_duration=fixed_durations[1],
                                                          phase_duration=phase_durations[1],
                                                          entry_function=entry_functions[1],
                                                          entry_function_args=entry_functions_args[1],
                                                          exit_function=exit_functions[1],
                                                          exit_function_args=exit_functions_args[1],
                                                          arrest_function=arrest_functions[1],
                                                          arrest_function_args=arrest_functions_args[1],
                                                          check_transition_to_next_phase_function=
                                                          check_transition_to_next_phase_functions[1],
                                                          check_transition_to_next_phase_function_args=
                                                          check_transition_to_next_phase_functions_args[
                                                              1], simulated_cell_volume=simulated_cell_volume,
                                                          cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[
                                                              1],
                                                          nuclear_volume_change_rate=nuclear_volume_change_rate[1],
                                                          calcification_rate=calcification_rate[1],
                                                          target_fluid_fraction=target_fluid_fraction[1],
                                                          nuclear_fluid=nuclear_fluid[1],
                                                          nuclear_solid=nuclear_solid[1],
                                                          nuclear_solid_target=nuclear_solid_target[1],
                                                          cytoplasm_fluid=cytoplasm_fluid[1],
                                                          cytoplasm_solid=cytoplasm_solid[1],
                                                          cytoplasm_solid_target=cytoplasm_solid_target[1],
                                                          target_cytoplasm_to_nuclear_ratio=
                                                          target_cytoplasm_to_nuclear_ratio[1],
                                                          calcified_fraction=calcified_fraction[1],
                                                          fluid_change_rate=fluid_change_rate[1],
                                                          user_phase_time_step=user_phases_time_step[1],
                                                          user_phase_time_step_args=user_phases_time_step_args[1])

        Ki67_positive_post = Phases.Ki67PositivePostMitotic(index=2, previous_phase_index=1, next_phase_index=0, dt=dt,
                                                            time_unit=time_unit, space_unit=space_unit,
                                                            division_at_phase_exit=division_at_phase_exits[2],
                                                            removal_at_phase_exit=removal_at_phase_exits[2],
                                                            fixed_duration=fixed_durations[2],
                                                            phase_duration=phase_durations[2],
                                                            entry_function=entry_functions[2],
                                                            entry_function_args=entry_functions_args[2],
                                                            exit_function=exit_functions[2],
                                                            exit_function_args=exit_functions_args[2],
                                                            arrest_function=arrest_functions[2],
                                                            arrest_function_args=arrest_functions_args[2],
                                                            check_transition_to_next_phase_function=
                                                            check_transition_to_next_phase_functions[2],
                                                            check_transition_to_next_phase_function_args=
                                                            check_transition_to_next_phase_functions_args[2],
                                                            simulated_cell_volume=simulated_cell_volume,
                                                            cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[
                                                                2],
                                                            nuclear_volume_change_rate=nuclear_volume_change_rate[2],
                                                            calcification_rate=calcification_rate[2],
                                                            target_fluid_fraction=target_fluid_fraction[2],
                                                            nuclear_fluid=nuclear_fluid[2],
                                                            nuclear_solid=nuclear_solid[2],
                                                            nuclear_solid_target=nuclear_solid_target[2],
                                                            cytoplasm_fluid=cytoplasm_fluid[2],
                                                            cytoplasm_solid=cytoplasm_solid[2],
                                                            cytoplasm_solid_target=cytoplasm_solid_target[2],
                                                            target_cytoplasm_to_nuclear_ratio=
                                                            target_cytoplasm_to_nuclear_ratio[2],
                                                            calcified_fraction=calcified_fraction[2],
                                                            fluid_change_rate=fluid_change_rate[2],
                                                            user_phase_time_step=user_phases_time_step[2],
                                                            user_phase_time_step_args=user_phases_time_step_args[2])
        phases = [Ki67_negative, Ki67_positive_pre, Ki67_positive_post]
        super().__init__(name=name, dt=dt, time_unit=time_unit, space_unit=space_unit, phases=phases,
                         senescent_phase=senescent_phase, user_phenotype_time_step=user_phenotype_time_step,
                         user_phenotype_time_step_args=user_phenotype_time_step_args)


class FlowCytometryBasic(Phenotype):
    """
    Inherits :class:`Phenotype`. Basic flow cytometry model.

    Three-phase live cell cycle consists of G0/G1 -> S -> G2/M -> (back to) G0/G1. Reference phases durations from
    https://www.ncbi.nlm.nih.gov/books/NBK9876/. G0/G1 phase is defined in :class:`Phases.G0G1` is more representative
    of the quiescent phase than the first growth phase, the cell transitions stochastically from this phase, its
    (default) expected duration is 5.15h, transition to the next phase is stochastic. S phase is defined in
    :class:`Phases.S` is the phase responsible for doubling the cell volume, its  (default) expected duration is 8h,
    transition to the next phase is stochastic. G2/M is defined in :class:`Phases.G2M`. It is the pre-mitotic rest
    phase. The cell divides when exiting this phase. Its (default) expected duration is 5h, transition to the next phase
    is stochastic
    """

    def __init__(self, name="Flow Cytometry Basic", dt=0.1, time_unit="min", space_unit="micrometer",
                 senescent_phase=False,
                 division_at_phase_exits=(False, False, True), removal_at_phase_exits=(False, False, False),
                 fixed_durations=(False, False, False), phase_durations: list = (5.15 * 60, 8 * 60.0, 5 * 60),
                 entry_functions=(None, None, None), entry_functions_args=(None, None, None),
                 exit_functions=(None, None, None), exit_functions_args=(None, None, None),
                 arrest_functions=(None, None, None), arrest_functions_args=(None, None, None),
                 check_transition_to_next_phase_functions=(None, None, None),
                 check_transition_to_next_phase_functions_args: list = (None, None, None), simulated_cell_volume=None,
                 cytoplasm_volume_change_rate=(None, None, None),
                 nuclear_volume_change_rate=(None, None, None), calcification_rate=(None, None, None),
                 calcified_fraction=(0, 0, 0),
                 target_fluid_fraction=(None, None, None), nuclear_fluid=(None, None, None),
                 nuclear_solid=(None, None, None),
                 nuclear_solid_target=(None, None, None), cytoplasm_fluid=(None, None, None),
                 cytoplasm_solid=(None, None, None),
                 cytoplasm_solid_target=(None, None, None), target_cytoplasm_to_nuclear_ratio=(None, None, None),
                 fluid_change_rate=(None, None, None),
                 user_phenotype_time_step=None, user_phenotype_time_step_args=None, user_phases_time_step=None,
                 user_phases_time_step_args=None):
        if user_phases_time_step is None:
            user_phases_time_step = [None, None, None]
            user_phases_time_step_args = [None, None, None]
        _check_arguments(3, name, division_at_phase_exits, removal_at_phase_exits, fixed_durations, phase_durations,
                         entry_functions, entry_functions_args, exit_functions, exit_functions_args, arrest_functions,
                         arrest_functions_args, check_transition_to_next_phase_functions,
                         check_transition_to_next_phase_functions_args,
                         cytoplasm_volume_change_rate, nuclear_volume_change_rate, calcification_rate,
                         calcified_fraction, target_fluid_fraction, nuclear_fluid, nuclear_solid, nuclear_solid_target,
                         cytoplasm_fluid, cytoplasm_solid, cytoplasm_solid_target, target_cytoplasm_to_nuclear_ratio,
                         fluid_change_rate)

        G0G1 = Phases.G0G1(dt=dt, time_unit=time_unit, space_unit=space_unit,
                           division_at_phase_exit=division_at_phase_exits[0],
                           removal_at_phase_exit=removal_at_phase_exits[0], fixed_duration=fixed_durations[0],
                           phase_duration=phase_durations[0], entry_function=entry_functions[0],
                           entry_function_args=entry_functions_args[0], exit_function=exit_functions[0],
                           exit_function_args=exit_functions_args[0], arrest_function=arrest_functions[0],
                           arrest_function_args=arrest_functions_args[0],
                           check_transition_to_next_phase_function=check_transition_to_next_phase_functions[0],
                           check_transition_to_next_phase_function_args=check_transition_to_next_phase_functions_args[
                               0],
                           simulated_cell_volume=simulated_cell_volume,
                           cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[0],
                           nuclear_volume_change_rate=nuclear_volume_change_rate[0],
                           calcification_rate=calcification_rate[0], target_fluid_fraction=target_fluid_fraction[0],
                           nuclear_fluid=nuclear_fluid[0], nuclear_solid=nuclear_solid[0],
                           nuclear_solid_target=nuclear_solid_target[0], cytoplasm_fluid=cytoplasm_fluid[0],
                           cytoplasm_solid=cytoplasm_solid[0], cytoplasm_solid_target=cytoplasm_solid_target[0],
                           target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                           calcified_fraction=calcified_fraction[0], fluid_change_rate=fluid_change_rate[0],
                           user_phase_time_step=user_phases_time_step[0],
                           user_phase_time_step_args=user_phases_time_step_args[0])

        S = Phases.S(dt=dt, time_unit=time_unit, space_unit=space_unit,
                     division_at_phase_exit=division_at_phase_exits[1],
                     removal_at_phase_exit=removal_at_phase_exits[1], fixed_duration=fixed_durations[1],
                     phase_duration=phase_durations[1], entry_function=entry_functions[1],
                     entry_function_args=entry_functions_args[1], exit_function=exit_functions[1],
                     exit_function_args=exit_functions_args[1], arrest_function=arrest_functions[1],
                     arrest_function_args=arrest_functions_args[1],
                     check_transition_to_next_phase_function=check_transition_to_next_phase_functions[1],
                     check_transition_to_next_phase_function_args=check_transition_to_next_phase_functions_args[1],
                     simulated_cell_volume=simulated_cell_volume,
                     cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[1],
                     nuclear_volume_change_rate=nuclear_volume_change_rate[1],
                     calcification_rate=calcification_rate[1], target_fluid_fraction=target_fluid_fraction[1],
                     nuclear_fluid=nuclear_fluid[1], nuclear_solid=nuclear_solid[1],
                     nuclear_solid_target=nuclear_solid_target[1], cytoplasm_fluid=cytoplasm_fluid[1],
                     cytoplasm_solid=cytoplasm_solid[1], cytoplasm_solid_target=cytoplasm_solid_target[1],
                     target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[1],
                     calcified_fraction=calcified_fraction[1], fluid_change_rate=fluid_change_rate[1],
                     user_phase_time_step=user_phases_time_step[1],
                     user_phase_time_step_args=user_phases_time_step_args[1])

        G2M = Phases.G2M(dt=dt, time_unit=time_unit, space_unit=space_unit,
                         division_at_phase_exit=division_at_phase_exits[2],
                         removal_at_phase_exit=removal_at_phase_exits[2], fixed_duration=fixed_durations[2],
                         phase_duration=phase_durations[2], entry_function=entry_functions[2],
                         entry_function_args=entry_functions_args[2], exit_function=exit_functions[2],
                         exit_function_args=exit_functions_args[2], arrest_function=arrest_functions[2],
                         arrest_function_args=arrest_functions_args[2],
                         check_transition_to_next_phase_function=check_transition_to_next_phase_functions[2],
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_functions_args[2],
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[2],
                         nuclear_volume_change_rate=nuclear_volume_change_rate[2],
                         calcification_rate=calcification_rate[2], target_fluid_fraction=target_fluid_fraction[2],
                         nuclear_fluid=nuclear_fluid[2], nuclear_solid=nuclear_solid[2],
                         nuclear_solid_target=nuclear_solid_target[2], cytoplasm_fluid=cytoplasm_fluid[2],
                         cytoplasm_solid=cytoplasm_solid[2], cytoplasm_solid_target=cytoplasm_solid_target[2],
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[2],
                         calcified_fraction=calcified_fraction[2], fluid_change_rate=fluid_change_rate[2],
                         user_phase_time_step=user_phases_time_step[2],
                         user_phase_time_step_args=user_phases_time_step_args[2])

        phases = [G0G1, S, G2M]

        super().__init__(name=name, dt=dt, time_unit=time_unit, space_unit=space_unit, phases=phases,
                         senescent_phase=senescent_phase, user_phenotype_time_step=user_phenotype_time_step,
                         user_phenotype_time_step_args=user_phenotype_time_step_args)


class FlowCytometryAdvanced(Phenotype):
    """
    Inherits :class:`Phenotype`. Flow cytometry model.

    Four-phase live cell cycle consists of G0/G1 -> S -> G2 -> M -> (back to) G0/G1. Reference phases durations from
    https://www.ncbi.nlm.nih.gov/books/NBK9876/. G0/G1 phase is defined in :class:`Phases.G0G1` is more representative
    of the quiescent phase than the first growth phase, the cell transitions stochastically from this phase, its
    (default) expected duration is 4.98h, transition to the next phase is stochastic. S phase is defined in
    :class:`Phases.S` is the phase responsible for doubling the cell volume, its  (default) expected duration is 8h,
    transition to the next phase is stochastic. The (default) cell volume growth rate is
    [total volume growth]/[phase duration]. G2 is defined in :class:`Phases.G0G1` (using different parameters than this
    phenotype G0/G1 phase). It is the pre-mitotic rest phase. Its (default) expected duration is 4h, transition to the
    next phase is stochastic. M is defined in :class:`Phases.G2M`, it is the mitotic phase, the cell divides when exi-
    ting this phase. Its expected duration is 1h, transition from this phase is stochastic.
    """

    def __init__(self, name="Flow Cytometry Advanced", dt=0.1, time_unit="min", space_unit="micrometer",
                 senescent_phase=False,
                 division_at_phase_exits=(False, False, False, True),
                 removal_at_phase_exits=(False, False, False, False), fixed_durations=(False, False, False, False),
                 phase_durations: list = (4.98 * 60, 8 * 60.0, 4 * 60, 1 * 60),
                 entry_functions=(None, None, None, None), entry_functions_args=(None, None, None, None),
                 exit_functions=(None, None, None, None), exit_functions_args=(None, None, None),
                 arrest_functions=(None, None, None, None), arrest_functions_args=(None, None, None, None),
                 check_transition_to_next_phase_functions=(None, None, None, None),
                 check_transition_to_next_phase_functions_args: list = (None, None, None, None),
                 simulated_cell_volume=None,
                 cytoplasm_volume_change_rate=(None, None, None, None),
                 nuclear_volume_change_rate=(None, None, None, None), calcification_rate=(None, None, None, None),
                 calcified_fraction=(0, 0, 0, 0),
                 target_fluid_fraction=(None, None, None, None), nuclear_fluid=(None, None, None, None),
                 nuclear_solid=(None, None, None, None),
                 nuclear_solid_target=(None, None, None, None), cytoplasm_fluid=(None, None, None, None),
                 cytoplasm_solid=(None, None, None, None),
                 cytoplasm_solid_target=(None, None, None, None),
                 target_cytoplasm_to_nuclear_ratio=(None, None, None, None),
                 fluid_change_rate=(None, None, None, None),
                 user_phenotype_time_step=None, user_phenotype_time_step_args=None, user_phases_time_step=None,
                 user_phases_time_step_args=None):
        if user_phases_time_step is None:
            user_phases_time_step = len(phase_durations) * [None]
            user_phases_time_step_args = len(phase_durations) * [None]
        _check_arguments(4, name, division_at_phase_exits, removal_at_phase_exits, fixed_durations, phase_durations,
                         entry_functions, entry_functions_args, exit_functions, exit_functions_args, arrest_functions,
                         arrest_functions_args, check_transition_to_next_phase_functions,
                         check_transition_to_next_phase_functions_args,
                         cytoplasm_volume_change_rate, nuclear_volume_change_rate, calcification_rate,
                         calcified_fraction, target_fluid_fraction, nuclear_fluid, nuclear_solid, nuclear_solid_target,
                         cytoplasm_fluid, cytoplasm_solid, cytoplasm_solid_target, target_cytoplasm_to_nuclear_ratio,
                         fluid_change_rate)

        G0G1 = Phases.G0G1(dt=dt, time_unit=time_unit, space_unit=space_unit,
                           division_at_phase_exit=division_at_phase_exits[0],
                           removal_at_phase_exit=removal_at_phase_exits[0], fixed_duration=fixed_durations[0],
                           phase_duration=phase_durations[0], entry_function=entry_functions[0],
                           entry_function_args=entry_functions_args[0], exit_function=exit_functions[0],
                           exit_function_args=exit_functions_args[0], arrest_function=arrest_functions[0],
                           arrest_function_args=arrest_functions_args[0],
                           check_transition_to_next_phase_function=check_transition_to_next_phase_functions[0],
                           check_transition_to_next_phase_function_args=check_transition_to_next_phase_functions_args[
                               0],
                           simulated_cell_volume=simulated_cell_volume,
                           cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[0],
                           nuclear_volume_change_rate=nuclear_volume_change_rate[0],
                           calcification_rate=calcification_rate[0], target_fluid_fraction=target_fluid_fraction[0],
                           nuclear_fluid=nuclear_fluid[0], nuclear_solid=nuclear_solid[0],
                           nuclear_solid_target=nuclear_solid_target[0], cytoplasm_fluid=cytoplasm_fluid[0],
                           cytoplasm_solid=cytoplasm_solid[0], cytoplasm_solid_target=cytoplasm_solid_target[0],
                           target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                           calcified_fraction=calcified_fraction[0], fluid_change_rate=fluid_change_rate[0],
                           user_phase_time_step=user_phases_time_step[0],
                           user_phase_time_step_args=user_phases_time_step_args[0])

        S = Phases.S(dt=dt, time_unit=time_unit, space_unit=space_unit,
                     division_at_phase_exit=division_at_phase_exits[1],
                     removal_at_phase_exit=removal_at_phase_exits[1], fixed_duration=fixed_durations[1],
                     phase_duration=phase_durations[1], entry_function=entry_functions[1],
                     entry_function_args=entry_functions_args[1], exit_function=exit_functions[1],
                     exit_function_args=exit_functions_args[1], arrest_function=arrest_functions[1],
                     arrest_function_args=arrest_functions_args[1],
                     check_transition_to_next_phase_function=check_transition_to_next_phase_functions[1],
                     check_transition_to_next_phase_function_args=check_transition_to_next_phase_functions_args[1],
                     simulated_cell_volume=simulated_cell_volume,
                     cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[1],
                     nuclear_volume_change_rate=nuclear_volume_change_rate[1],
                     calcification_rate=calcification_rate[1], target_fluid_fraction=target_fluid_fraction[1],
                     nuclear_fluid=nuclear_fluid[1], nuclear_solid=nuclear_solid[1],
                     nuclear_solid_target=nuclear_solid_target[1], cytoplasm_fluid=cytoplasm_fluid[1],
                     cytoplasm_solid=cytoplasm_solid[1], cytoplasm_solid_target=cytoplasm_solid_target[1],
                     target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[1],
                     fluid_change_rate=fluid_change_rate[1],
                     user_phase_time_step=user_phases_time_step[1],
                     user_phase_time_step_args=user_phases_time_step_args[1])

        G2 = Phases.G0G1(index=2, previous_phase_index=1, next_phase_index=3, dt=dt, time_unit=time_unit,
                         space_unit=space_unit, name="G2",
                         division_at_phase_exit=division_at_phase_exits[2],
                         removal_at_phase_exit=removal_at_phase_exits[2], fixed_duration=fixed_durations[2],
                         phase_duration=phase_durations[2], entry_function=entry_functions[2],
                         entry_function_args=entry_functions_args[2], exit_function=exit_functions[2],
                         exit_function_args=exit_functions_args[2], arrest_function=arrest_functions[2],
                         arrest_function_args=arrest_functions_args[2],
                         check_transition_to_next_phase_function=check_transition_to_next_phase_functions[2],
                         check_transition_to_next_phase_function_args=check_transition_to_next_phase_functions_args[2],
                         simulated_cell_volume=simulated_cell_volume,
                         cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[2],
                         nuclear_volume_change_rate=nuclear_volume_change_rate[2],
                         calcification_rate=calcification_rate[2], target_fluid_fraction=target_fluid_fraction[2],
                         nuclear_fluid=nuclear_fluid[2], nuclear_solid=nuclear_solid[2],
                         nuclear_solid_target=nuclear_solid_target[2], cytoplasm_fluid=cytoplasm_fluid[2],
                         cytoplasm_solid=cytoplasm_solid[2], cytoplasm_solid_target=cytoplasm_solid_target[2],
                         target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[2],
                         fluid_change_rate=fluid_change_rate[2],
                         user_phase_time_step=user_phases_time_step[2],
                         user_phase_time_step_args=user_phases_time_step_args[2])

        M = Phases.G2M(index=3, previous_phase_index=2, next_phase_index=0, dt=dt, time_unit=time_unit,
                       space_unit=space_unit, name="M",
                       division_at_phase_exit=division_at_phase_exits[3],
                       removal_at_phase_exit=removal_at_phase_exits[3], fixed_duration=fixed_durations[3],
                       phase_duration=phase_durations[3], entry_function=entry_functions[3],
                       entry_function_args=entry_functions_args[3], exit_function=exit_functions[3],
                       exit_function_args=exit_functions_args[3], arrest_function=arrest_functions[3],
                       arrest_function_args=arrest_functions_args[3],
                       check_transition_to_next_phase_function=check_transition_to_next_phase_functions[3],
                       check_transition_to_next_phase_function_args=check_transition_to_next_phase_functions_args[3],
                       simulated_cell_volume=simulated_cell_volume,
                       cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[3],
                       nuclear_volume_change_rate=nuclear_volume_change_rate[3],
                       calcification_rate=calcification_rate[3], target_fluid_fraction=target_fluid_fraction[3],
                       nuclear_fluid=nuclear_fluid[3], nuclear_solid=nuclear_solid[3],
                       nuclear_solid_target=nuclear_solid_target[3], cytoplasm_fluid=cytoplasm_fluid[3],
                       cytoplasm_solid=cytoplasm_solid[3], cytoplasm_solid_target=cytoplasm_solid_target[3],
                       target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[3],
                       fluid_change_rate=fluid_change_rate[3],
                       user_phase_time_step=user_phases_time_step[3],
                       user_phase_time_step_args=user_phases_time_step_args[3])

        phases = [G0G1, S, G2, M]

        super().__init__(name=name, dt=dt, time_unit=time_unit, phases=phases, space_unit=space_unit,
                         senescent_phase=senescent_phase, user_phenotype_time_step=user_phenotype_time_step,
                         user_phenotype_time_step_args=user_phenotype_time_step_args)


class ApoptosisStandard(Phenotype):
    """
    Inherits :class:`Phenotype`. The standard apoptotic model.

    A single phase phenotype, tje apoptotic phase, defined in :class:Phases.Apoptosis`. The phase has fixed length,
    and reduces the cell volume, the cytoplasm diminishes at a rate of 1/60 1/min, the nucleus at a rate of 0.35/60
    1/min. The (already dead) cell should be removed from the simulation once this phase is concluded. By default, this
    phase duration is fixed at 8.6h

    """

    def __init__(self, name="Standard apoptosis model", dt=0.1, time_unit="min", space_unit="micrometer",
                 senescent_phase=False,
                 division_at_phase_exits=(False,), removal_at_phase_exits=(True,), fixed_durations=(True,),
                 phase_durations=(8.6 * 60,), entry_functions=(None,), entry_functions_args=(None,),
                 exit_functions=(None,), exit_functions_args=(None,), arrest_functions=(None,),
                 arrest_functions_args=(None,), check_transition_to_next_phase_functions=(None,),
                 check_transition_to_next_phase_functions_args=(None,), simulated_cell_volume=None,
                 cytoplasm_volume_change_rate=(1 / 60,),
                 nuclear_volume_change_rate=(0.35 / 60,), calcification_rate=(0,),
                 calcified_fraction=(0,),
                 target_fluid_fraction=(None,), nuclear_fluid=(None,), nuclear_solid=(None,),
                 nuclear_solid_target=(None,), cytoplasm_fluid=(None,), cytoplasm_solid=(None,),
                 cytoplasm_solid_target=(None,), target_cytoplasm_to_nuclear_ratio=(None,),
                 fluid_change_rate=(None,),
                 user_phenotype_time_step=None, user_phenotype_time_step_args=None, user_phases_time_step=None,
                 user_phases_time_step_args=None):
        if user_phases_time_step is None:
            user_phases_time_step = len(phase_durations) * [None]
            user_phases_time_step_args = len(phase_durations) * [None]
        _check_arguments(1, name, division_at_phase_exits, removal_at_phase_exits, fixed_durations, phase_durations,
                         entry_functions, entry_functions_args, exit_functions, exit_functions_args, arrest_functions,
                         arrest_functions_args, check_transition_to_next_phase_functions,
                         check_transition_to_next_phase_functions_args,
                         cytoplasm_volume_change_rate, nuclear_volume_change_rate, calcification_rate,
                         calcified_fraction, target_fluid_fraction, nuclear_fluid, nuclear_solid, nuclear_solid_target,
                         cytoplasm_fluid, cytoplasm_solid, cytoplasm_solid_target, target_cytoplasm_to_nuclear_ratio,
                         fluid_change_rate)

        apopto = Phases.Apoptosis(index=0, previous_phase_index=0, next_phase_index=0, dt=dt, time_unit=time_unit,
                                  space_unit=space_unit,
                                  name="Apoptosis", division_at_phase_exit=division_at_phase_exits[0],
                                  removal_at_phase_exit=removal_at_phase_exits[0], fixed_duration=fixed_durations[0],
                                  phase_duration=phase_durations[0], entry_function=entry_functions[0],
                                  entry_function_args=entry_functions_args[0], exit_function=exit_functions[0],
                                  exit_function_args=exit_functions_args[0], arrest_function=arrest_functions[0],
                                  arrest_function_args=arrest_functions_args[0],
                                  check_transition_to_next_phase_function=check_transition_to_next_phase_functions[0],
                                  check_transition_to_next_phase_function_args=
                                  check_transition_to_next_phase_functions_args[0],
                                  simulated_cell_volume=simulated_cell_volume,
                                  cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[0],
                                  nuclear_volume_change_rate=nuclear_volume_change_rate[0],
                                  calcification_rate=calcification_rate[0],
                                  target_fluid_fraction=target_fluid_fraction[0], nuclear_fluid=nuclear_fluid[0],
                                  nuclear_solid=nuclear_solid[0], nuclear_solid_target=nuclear_solid_target[0],
                                  cytoplasm_fluid=cytoplasm_fluid[0], cytoplasm_solid=cytoplasm_solid[0],
                                  cytoplasm_solid_target=cytoplasm_solid_target[0],
                                  target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                                  calcified_fraction=calcified_fraction[0], fluid_change_rate=fluid_change_rate[0],
                                  user_phase_time_step=user_phases_time_step[0],
                                  user_phase_time_step_args=user_phases_time_step_args[0])

        # a phase to help lyse the simulated cell, shouldn't do anything
        # debris = Phases.Phase(index=1, previous_phase_index=0, next_phase_index=1, dt=dt, time_unit=time_unit,
        #                       name="Debris", division_at_phase_exit=False, removal_at_phase_exit=True,
        #                       fixed_duration=True, phase_duration=1e6)

        phases = [apopto]

        super().__init__(name=name, dt=dt, time_unit=time_unit, space_unit=space_unit,
                         phases=phases, senescent_phase=senescent_phase,
                         user_phenotype_time_step=user_phenotype_time_step,
                         user_phenotype_time_step_args=user_phenotype_time_step_args)


class NecrosisStandard(Phenotype):
    """
    Inherits :class:`Phenotype`. Standard Necrosis model

    Two-phase phenotype. The first phase represents the osmotic swell of a necrotic cell, defined in
    :class:`Phases.NecrosisSwell` it doesn't have a set or expected phase duration, transition to the next phase happens
    when the cell reaches its rupture volume (twice the original volume, by default). Osmotic swelling rates and calci-
    fication rate in :class:`Phases.NecrosisSwell`. Second phase is the dissolution of the ruptured cell into its media,
    defined in :class:`Phases.NecrosisLysed`, dissolving rates can be found there. As a safeguard, this phase has a
    fixed duration of 60 days, after which, if the cell hasn't dissolved naturally, it should be removed from the si-
    mulation.

    """

    def __init__(self, name="Standard necrosis model", dt=0.1, time_unit="min", space_unit="micrometer",
                 senescent_phase=False,
                 division_at_phase_exits=(False, False), removal_at_phase_exits=(False, True),
                 fixed_durations=(False, True),
                 phase_durations=(None, None), entry_functions=(None, None), entry_functions_args=(None, None),
                 exit_functions=(None, None), exit_functions_args=(None, None), arrest_functions=(None, None),
                 arrest_functions_args=(None, None), check_transition_to_next_phase_functions=(None, None),
                 check_transition_to_next_phase_functions_args=(None, None), simulated_cell_volume=None,
                 cytoplasm_volume_change_rate=(None, None),
                 nuclear_volume_change_rate=(None, None),
                 calcification_rate=(None, None),
                 calcified_fraction=(0, 0),
                 target_fluid_fraction=(None, None), nuclear_fluid=(None, None), nuclear_solid=(None, None),
                 nuclear_solid_target=(None, None), cytoplasm_fluid=(None, None), cytoplasm_solid=(None, None),
                 cytoplasm_solid_target=(None, None), target_cytoplasm_to_nuclear_ratio=(None, None),
                 fluid_change_rate=(None, None),
                 user_phenotype_time_step=None, user_phenotype_time_step_args=None, user_phases_time_step=None,
                 user_phases_time_step_args=None):
        if user_phases_time_step is None:
            user_phases_time_step = len(phase_durations) * [None]
            user_phases_time_step_args = len(phase_durations) * [None]
        _check_arguments(2, name, division_at_phase_exits, removal_at_phase_exits, fixed_durations, phase_durations,
                         entry_functions, entry_functions_args, exit_functions, exit_functions_args, arrest_functions,
                         arrest_functions_args, check_transition_to_next_phase_functions,
                         check_transition_to_next_phase_functions_args,
                         cytoplasm_volume_change_rate, nuclear_volume_change_rate, calcification_rate,
                         calcified_fraction, target_fluid_fraction, nuclear_fluid, nuclear_solid, nuclear_solid_target,
                         cytoplasm_fluid, cytoplasm_solid, cytoplasm_solid_target, target_cytoplasm_to_nuclear_ratio,
                         fluid_change_rate)

        necro_swell = Phases.NecrosisSwell(index=0, previous_phase_index=0, next_phase_index=1, dt=dt,
                                           time_unit=time_unit, space_unit=space_unit,
                                           division_at_phase_exit=division_at_phase_exits[0],
                                           removal_at_phase_exit=removal_at_phase_exits[0],
                                           fixed_duration=fixed_durations[0], phase_duration=phase_durations[0],
                                           entry_function=entry_functions[0],
                                           entry_function_args=entry_functions_args[0], exit_function=exit_functions[0],
                                           exit_function_args=exit_functions_args[0],
                                           arrest_function=arrest_functions[0],
                                           arrest_function_args=arrest_functions_args[0],
                                           check_transition_to_next_phase_function=
                                           check_transition_to_next_phase_functions[0],
                                           check_transition_to_next_phase_function_args=
                                           check_transition_to_next_phase_functions_args[0],
                                           simulated_cell_volume=simulated_cell_volume,
                                           cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[0],
                                           nuclear_volume_change_rate=nuclear_volume_change_rate[0],
                                           calcification_rate=calcification_rate[0],
                                           target_fluid_fraction=target_fluid_fraction[0],
                                           nuclear_fluid=nuclear_fluid[0], nuclear_solid=nuclear_solid[0],
                                           nuclear_solid_target=nuclear_solid_target[0],
                                           cytoplasm_fluid=cytoplasm_fluid[0], cytoplasm_solid=cytoplasm_solid[0],
                                           cytoplasm_solid_target=cytoplasm_solid_target[0],
                                           target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[0],
                                           calcified_fraction=calcified_fraction[0],
                                           fluid_change_rate=fluid_change_rate[0],
                                           user_phase_time_step=user_phases_time_step[0],
                                           user_phase_time_step_args=user_phases_time_step_args[0])

        necro_lysed = Phases.NecrosisLysed(index=1, previous_phase_index=0, next_phase_index=1, dt=dt,
                                           time_unit=time_unit, space_unit=space_unit,
                                           division_at_phase_exit=division_at_phase_exits[1],
                                           removal_at_phase_exit=removal_at_phase_exits[1],
                                           fixed_duration=fixed_durations[1], phase_duration=phase_durations[1],
                                           entry_function=entry_functions[1],
                                           entry_function_args=entry_functions_args[1], exit_function=exit_functions[1],
                                           exit_function_args=exit_functions_args[1],
                                           arrest_function=arrest_functions[1],
                                           arrest_function_args=arrest_functions_args[1],
                                           check_transition_to_next_phase_function=
                                           check_transition_to_next_phase_functions[1],
                                           check_transition_to_next_phase_function_args=
                                           check_transition_to_next_phase_functions_args[1],
                                           simulated_cell_volume=simulated_cell_volume,
                                           cytoplasm_volume_change_rate=cytoplasm_volume_change_rate[1],
                                           nuclear_volume_change_rate=nuclear_volume_change_rate[1],
                                           calcification_rate=calcification_rate[1],
                                           target_fluid_fraction=target_fluid_fraction[1],
                                           nuclear_fluid=nuclear_fluid[1], nuclear_solid=nuclear_solid[1],
                                           nuclear_solid_target=nuclear_solid_target[1],
                                           cytoplasm_fluid=cytoplasm_fluid[1], cytoplasm_solid=cytoplasm_solid[1],
                                           cytoplasm_solid_target=cytoplasm_solid_target[1],
                                           target_cytoplasm_to_nuclear_ratio=target_cytoplasm_to_nuclear_ratio[1],
                                           calcified_fraction=calcified_fraction[1],
                                           fluid_change_rate=fluid_change_rate[1],
                                           user_phase_time_step=user_phases_time_step[1],
                                           user_phase_time_step_args=user_phases_time_step_args[1])

        phases = [necro_swell, necro_lysed]

        super().__init__(name=name, dt=dt, time_unit=time_unit, space_unit=space_unit, phases=phases,
                         senescent_phase=senescent_phase, user_phenotype_time_step=user_phenotype_time_step,
                         user_phenotype_time_step_args=user_phenotype_time_step_args)

        return


cycle_names = ["Simple Live", "Ki67 Basic", "Ki67 Advanced", "Flow Cytometry Basic", "Flow Cytometry Advanced",
               "Standard apoptosis model", "Standard necrosis model"]


def get_phenotype_by_name(name):
    """
    Fetches a (uninitialized) phenotype class

    :param name: Name of the phenotype model being fetched
    :type name: str
    :return: A phenotype model
    :rtype: :class:`Phenotype`
    """
    if name not in cycle_names:
        raise ValueError(f"{name} is not a pre-defined cycle")

    if name == "Simple Live":
        return SimpleLiveCycle
    elif name == "Ki67 Basic":
        return Ki67Basic
    elif name == "Ki67 Advanced":
        return Ki67Advanced
    elif name == "Flow Cytometry Basic":
        return FlowCytometryBasic
    elif name == "Flow Cytometry Advanced":
        return FlowCytometryAdvanced
    elif name == "Standard apoptosis model":
        return ApoptosisStandard
    elif name == "Standard necrosis model":
        return NecrosisStandard

    return Phenotype


if __name__ == "__main__":
    import numpy as np

    dt = 1
    print(cycle_names)

    test = Ki67Basic(dt=dt)
    cells = [type('', (), {})() for _ in range(2)]
    for c in cells:
        c.p = test.copy()

    for c in cells:
        print(c.p)
        [print(ph) for ph in c.p.phases]

    custom_p0 = Phases.Phase(index=0, previous_phase_index=-1, next_phase_index=1, dt=dt,
                             time_unit="min", space_unit="micrometer", name="custom_p0",
                             division_at_phase_exit=False, removal_at_phase_exit=False, fixed_duration=True,
                             phase_duration=20, entry_function=None, exit_function=None,
                             arrest_function=None,
                             check_transition_to_next_phase_function=None,
                             simulated_cell_volume=1,
                             cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None,
                             calcification_rate=None, target_fluid_fraction=None, nuclear_fluid=None,
                             nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                             cytoplasm_solid=None, cytoplasm_solid_target=None,
                             target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None,
                             fluid_change_rate=None, relative_rupture_volume=None,
                             user_phase_time_step=None, user_phase_time_step_args=(None,))


    def grow_phase_transition(*args):
        volume = args[0]
        doubling_volume = 0.9 * args[1]
        time_phase = args[2]
        phase_duration = args[3]
        return volume >= doubling_volume and time_phase > phase_duration


    custom_p1 = Phases.Ki67Positive(index=1, previous_phase_index=0, next_phase_index=2, dt=dt,
                                    time_unit="min", space_unit="micrometer", name="custom_p1",
                                    division_at_phase_exit=False, removal_at_phase_exit=False, fixed_duration=True,
                                    phase_duration=120, entry_function=None,
                                    entry_function_args=[None],
                                    exit_function=False, arrest_function=None,
                                    check_transition_to_next_phase_function=grow_phase_transition,
                                    check_transition_to_next_phase_function_args=[0, 9, 0, 9],
                                    simulated_cell_volume=1,
                                    cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None,
                                    calcification_rate=None, target_fluid_fraction=None, nuclear_fluid=None,
                                    nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                                    cytoplasm_solid=None, cytoplasm_solid_target=None,
                                    target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None,
                                    fluid_change_rate=None, relative_rupture_volume=None,
                                    user_phase_time_step=None, user_phase_time_step_args=(None,))

    stable_phase_1 = Phases.Phase(index=2, previous_phase_index=1, next_phase_index=3, dt=dt,
                                  time_unit="min", space_unit="micrometer", name="stable1",
                                  division_at_phase_exit=False, removal_at_phase_exit=False, fixed_duration=True,
                                  phase_duration=30, entry_function=None, exit_function=None,
                                  arrest_function=None,
                                  check_transition_to_next_phase_function=None,
                                  simulated_cell_volume=1,
                                  cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None,
                                  calcification_rate=None, target_fluid_fraction=None, nuclear_fluid=None,
                                  nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                                  cytoplasm_solid=None, cytoplasm_solid_target=None,
                                  target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None,
                                  fluid_change_rate=None, relative_rupture_volume=None,
                                  user_phase_time_step=None, user_phase_time_step_args=(None,))


    def shrink_phase_transition(*args):
        dt = args[0]
        phase_duration = args[1]
        total = args[2]
        total_target = args[3]
        time_check = np.random.uniform() < (1 - np.exp(-dt / phase_duration))
        volume_check = total <= 1.1 * total_target
        return time_check and volume_check


    shrink_phase = Phases.Ki67PositivePostMitotic(index=3, previous_phase_index=2, next_phase_index=0, dt=dt,
                                                  time_unit="min", space_unit="micrometer", name="shrink",
                                                  division_at_phase_exit=False, removal_at_phase_exit=False,
                                                  fixed_duration=False,
                                                  phase_duration=60, entry_function=None,
                                                  entry_function_args=[None], exit_function=None, arrest_function=None,
                                                  check_transition_to_next_phase_function=shrink_phase_transition,
                                                  check_transition_to_next_phase_function_args=[0, 1, 99, 0],
                                                  simulated_cell_volume=1,
                                                  cytoplasm_volume_change_rate=None, nuclear_volume_change_rate=None,
                                                  calcification_rate=None, target_fluid_fraction=None,
                                                  nuclear_fluid=None,
                                                  nuclear_solid=None, nuclear_solid_target=None, cytoplasm_fluid=None,
                                                  cytoplasm_solid=None, cytoplasm_solid_target=None,
                                                  target_cytoplasm_to_nuclear_ratio=None, calcified_fraction=None,
                                                  fluid_change_rate=None, relative_rupture_volume=None,
                                                  user_phase_time_step=None, user_phase_time_step_args=(None,))
    shrink_phase.check_transition_to_next_phase_function = shrink_phase_transition
    # shrink_phase.entry_function = Phases.Phase._halve_target_volume

    custom_pheno = Phenotype(name="oscillate volume with rests", dt=dt, time_unit="min",
                             space_unit="micrometer",
                             phases=[custom_p0, custom_p1, stable_phase_1,
                                     shrink_phase],
                             senescent_phase=False, starting_phase_index=0,
                             user_phenotype_time_step=None,
                             user_phenotype_time_step_args=[None, ])

    for i in range(10000):
        print(f"t={i}")
        changed_phase, died, divides = test.time_step_phenotype()
        print(changed_phase, died, divides)
        if custom_pheno.current_phase.name == "custom_p1":
            custom_pheno.current_phase.check_transition_to_next_phase_function_args = \
                [custom_pheno.current_phase.volume.total, custom_pheno.current_phase.volume.total_target,
                 custom_pheno.current_phase.time_in_phase, custom_pheno.current_phase.phase_duration]
        elif custom_pheno.current_phase.name == "shrink":
            custom_pheno.current_phase.check_transition_to_next_phase_function_args = \
                [dt, custom_pheno.current_phase.phase_duration,
                 custom_pheno.current_phase.volume.total, custom_pheno.current_phase.volume.total_target]

        changed_phase, died, divides = custom_pheno.time_step_phenotype()
        print(changed_phase, died, divides, custom_pheno.current_phase.volume.total)
