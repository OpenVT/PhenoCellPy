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

from numpy import exp
from numpy import log
from numpy.random import uniform

from phenocellpy.cell_volume import CellVolumes
from phenocellpy.types import PhaseConfig, TimingConfig, VolumeConfig, VolumeRatesConfig, EventConfig, FunctionsConfig

from copy import deepcopy



class Phase:
    """
    Base class to define phases of a cell phenotype.

    This base class is inherited by all phases. It defines some methods to time-step the Phase model, to transition to
    the next phase, and what should be done by the phenotype model when entering/exiting a phase. Initializes the cell
    volume model.

    Methods:
    -------

    time_step_phase()
        Time steps the phase. Returns a tuple (did the cell transition to the next phase, did the cell enter
        senescence). See `time_step_phase`'s documentation for further explanation.

    transition_to_next_phase(*args)
        One of the default transition functions (`_transition_to_next_phase_deterministic`,
        `_transition_to_next_phase_stochastic`) or a user defined function. If user defined, it will be called with
        `check_transition_to_next_phase_function_args` as args. Must return a bool denoting if the transition occurs or not.

    _transition_to_next_phase_deterministic()
        Default deterministic transition function. Returns `time_in_phase > phase_duration`

    _transition_to_next_phase_stochastic()
        Default stochastic transition function. Probability of transition depends on `dt` and `phase_duration`

    entry_function(*args)
        Optional function to be executed upon entering this phase. Some pre-build Phases have their own entry function
        already defined. It gets called using attribute `entry_function_args`. Must have no return

    exit_function(*args)
        Optional function to be executed just before exiting this phase. Some pre-build Phases have their own exit
        function already defined. It gets called using attribute `exit_function_args`. Must have no return

    arrest_function(*args)
        Optional function that returns true if the cell should exit the cell cycle and enter senescence

    user_phase_time_step(*args)
        User-defined function to be executed with the time-step


    update_volume()
        Function to update the volume of the cell. Calls the cell volume submodel (class:CellVolume) `update_volume`
        function

    _double_target_volume()
        Function that doubles all the target volumes of subclass class:CellVolume. Is used by some inheriting Phases as
        the entry function

    Parameters:
    -----------

    See :meth:`__init__`. The phase parameters are described in :class:`phenocellpy.types.PhaseConfig` and
    the configs it groups (:class:`TimingConfig`, :class:`VolumeConfig`, :class:`VolumeRatesConfig`,
    :class:`EventConfig`, :class:`FunctionsConfig`)

    Attributes
    ----------

    time_in_phase : float
        Time spent in this phase

    volume : class:cell_volume.CellVolumes
        Cell volume submodel

    """

    def __init__(self, config: PhaseConfig = None, dt: float = None, time_unit: str = "min",
                 space_unit: str = "micrometer"):
        """
        :param config: Parameters of the phase. If None, :meth:`default_config` is used. Parameters left as None take
            the value from :meth:`default_config`
        :type config: PhaseConfig
        :param dt: Time step duration in units of `time_unit`. `dt > 0`
        :type dt: float
        :param time_unit: What time units are used by the model (e.g., minutes, hours, days)
        :type time_unit: str
        :param space_unit: What space units are used by the model
        :type space_unit: str
        """
        defaults = self.default_config()
        if config is None:
            config = defaults
        elif not isinstance(config, PhaseConfig):
            raise TypeError(f"`config` must be a PhaseConfig. Got {type(config).__name__}.")
        else:
            config = config.fill_none_from(defaults)

        # assignments only run the single-field checks, the functions might have changed since the config was built
        config.functions.validate()

        self.config = config

        self.index = config.index

        self.previous_phase_index = config.previous_phase_index

        self.next_phase_index = config.next_phase_index

        self.time_unit = time_unit
        self.space_unit = space_unit

        if dt is None or dt <= 0:
            raise ValueError(f"'dt' must be greater than 0. Got {dt}.")
        self.dt = dt

        self.name = config.name

        self.division_at_phase_exit = config.events.division_at_phase_exit  # bool flagging for division
        self.removal_at_phase_exit = config.events.removal_at_phase_exit  # bool flagging for removal (e.g., death)

        self.fixed_duration = config.timing.fixed_duration

        self.phase_duration = config.timing.phase_duration

        self.time_in_phase = 0

        functions = config.functions

        self.entry_function = functions.entry_function  # function to be executed upon entering this phase
        self.entry_function_args = functions.entry_function_args

        self.exit_function = functions.exit_function  # function to be executed just before exiting this phase
        self.exit_function_args = functions.exit_function_args

        self.arrest_function = functions.arrest_function  # function determining if cell will exit cell cycle and become senescent
        self.arrest_function_args = functions.arrest_function_args

        if functions.check_transition_to_next_phase_function is None:
            self.check_transition_to_next_phase_function_args = [None]
            if self.fixed_duration:
                self.check_transition_to_next_phase_function = self._check_transition_to_next_phase_deterministic
            else:
                self.check_transition_to_next_phase_function = self._check_transition_to_next_phase_stochastic
        else:
            self.check_transition_to_next_phase_function_args = functions.check_transition_to_next_phase_function_args
            self.check_transition_to_next_phase_function = functions.check_transition_to_next_phase_function

        if config.volume.simulated_cell_volume is None:
            self.simulated_cell_volume = 1
        else:
            self.simulated_cell_volume = config.volume.simulated_cell_volume

        rates = config.volume.rates

        # the default rates are reference values for MCF-7, in 1/min
        if rates.cytoplasm_volume_change_rate is None:
            self.cytoplasm_volume_change_rate = 0.27 / 60.0
        else:
            self.cytoplasm_volume_change_rate = rates.cytoplasm_volume_change_rate

        if rates.nuclear_volume_change_rate is None:
            self.nuclear_volume_change_rate = 0.33 / 60.0
        else:
            self.nuclear_volume_change_rate = rates.nuclear_volume_change_rate
        if rates.calcification_rate is None:
            self.calcification_rate = 0
        else:
            self.calcification_rate = rates.calcification_rate

        if rates.fluid_change_rate is None:
            self.fluid_change_rate = 3.0 / 60.0
        else:
            self.fluid_change_rate = rates.fluid_change_rate

        self.user_phase_time_step = functions.user_phase_time_step

        self.user_phase_time_step_args = functions.user_phase_time_step_args

        self.volume = CellVolumes(config.volume)

    @classmethod
    def default_config(cls):
        """
        Default parameters of this phase class. Parameters a config leaves as None take their value from here.

        :rtype: PhaseConfig
        """
        return PhaseConfig(functions=FunctionsConfig(user_phase_time_step_args=(None,)))

    def update_volume(self):
        """
        Calls the cell volume submodel :function:`CellVolumes.update_volume`. Passes the current phase volume change
        rates as well as the timestep to it.

        :return: No return
        """
        self.volume.update_volume(
            self.dt, self.fluid_change_rate, self.nuclear_volume_change_rate, self.cytoplasm_volume_change_rate,
            self.calcification_rate
        )

    def _check_transition_to_next_phase_stochastic(self, *none):
        """
        Default stochastic phase transition function.

        Calculates a Poisson probability based on dt and self.phase_duration (p=1-exp(-dt/phase_duration), rolls a
        random number, and returns random number < probability.

        :param none: Not used. Placeholder in case of user defined function with args
        :return: bool. random number < probability of transition
        """

        # the approximation 1-exp(-x) ~ x can be used. That approximation has a difference of 0.005 at x=0.1, which I'd
        # find acceptable. TODO: implement a check on self.dt / self.phase_duration, if it is < .1 use the approximation

        prob = float(1 - exp(-self.dt / self.phase_duration))
        return uniform() < prob

    def _check_transition_to_next_phase_deterministic(self, *none):
        """
        Default deterministic phase transition function.

        If the time spent in this phase is greater than the phase duration, go to the next phase.

        :param none: Not used. Placeholder in case of user defined function with args
        :return:
        """
        return self.time_in_phase > self.phase_duration

    def time_step_phase(self):
        """

        Time steps the phase.

        This function increments the `time_in_phase` by `dt`. Updates the cell volume. Checks if the cell arrests its
        cycle (i.e., leaves the cycle; goes to senescence). If the cell doesn't senesce, this function checks if the
        cell
        should transition to the next phase.

        :return: tuple. First element of tuple: bool denoting if the cell moves to the next phase. Second element:
        denotes if the cell leaves the cell cycle and enters senescence.
        """
        self.time_in_phase += self.dt

        self.update_volume()

        transition_to_index = None

        if self.user_phase_time_step is not None:
            self.user_phase_time_step(*self.user_phase_time_step_args)

        if self.arrest_function is not None:
            exit_phenotype = self.arrest_function(*self.exit_function_args)
            go_to_next_phase_in_phenotype = False
            return go_to_next_phase_in_phenotype, exit_phenotype, transition_to_index
        else:
            exit_phenotype = False

        go_to_next_phase_in_phenotype = self.check_transition_to_next_phase_function(
            *self.check_transition_to_next_phase_function_args
        )

        if hasattr(go_to_next_phase_in_phenotype, "len") and len(go_to_next_phase_in_phenotype) > 1:
            transition_to_index = go_to_next_phase_in_phenotype[1]
            go_to_next_phase_in_phenotype = go_to_next_phase_in_phenotype[0]

        if go_to_next_phase_in_phenotype and self.exit_function is not None:
            self.exit_function(*self.exit_function_args)
            return go_to_next_phase_in_phenotype, exit_phenotype, transition_to_index
        return go_to_next_phase_in_phenotype, exit_phenotype, transition_to_index

    def _double_target_volume(self, *none):
        """

        Doubles the cell volume submodel (:class:`phenocellpy.cell_volume`) target volumes. Used by several cell cycle
        models to double the cell volume before mitosis

        :param none: Not used. This is a custom entry function, therefore it has to have args
        :return: No return
        """
        self.volume.nuclear_solid_target *= 2
        self.volume.cytoplasm_solid_target *= 2

    def _halve_target_volume(self, *none):
        """

        Halves the cell volume submodel (:class:`phenocellpy.cell_volume`) target volumes. Used by several cell cycle
        models to halve the cell volume after mitosis

        :param none: Not used. This is a custom entry function, therefore it has to have args
        :return: No return
        """
        self.volume.cytoplasm_solid_target /= 2
        self.volume.nuclear_solid_target /= 2  # 540

    def copy(self):
        return deepcopy(self)

    def __str__(self):
        return f"{self.name} phase, at memory {self.__repr__().split(' ')[-1][:-1]}"

    @property
    def _short_str(self):
        return f"{self.name} phase"


class SenescentPhase(Phase):
    """
    Default Senescent Phase. Inherits :class:`Phase`

    This senescent phase class is meant to be "outside" whatever phenotype progression is being used.

    """

    def __init__(self, config: PhaseConfig = None, dt: float = None, time_unit: str = "min",
                 space_unit: str = "micrometer"):
        super().__init__(config, dt=dt, time_unit=time_unit, space_unit=space_unit)
        return

    @classmethod
    def default_config(cls):
        return PhaseConfig(
            name="senescent", index=9999, next_phase_index=9999,
            timing=TimingConfig(phase_duration=60 * 24 * 60, fixed_duration=True),
            volume=VolumeConfig(rates=VolumeRatesConfig(cytoplasm_volume_change_rate=0, nuclear_volume_change_rate=0,
                                                    calcification_rate=0)),
        )


class Ki67Negative(Phase):
    """
    Inherits :class:`Phase`. Defines Ki 67- quiescent phase.

    This is a quiescent phenotype for cells that are replicating. Ki67 is a protein marker associated with
    proliferation.
    Transition to the next phase is set to be stochastic (the phase does not use a fixed duration) by default. Default
    expected phase duration is 4.59h, the phase transition rate is, therefore, dt/4.59 1/h.
    This phase does not calcify the cell. The parameters for this phase are based on the MCF-10A cell line
    https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-10a-cell-line
    https://www.ebi.ac.uk/ols/ontologies/bto/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FBTO_0001939
    """

    def __init__(self, config: PhaseConfig = None, dt: float = 0.1, time_unit: str = "min",
                 space_unit: str = "micrometer"):
        super().__init__(config, dt=dt, time_unit=time_unit, space_unit=space_unit)

    @classmethod
    def default_config(cls):
        return PhaseConfig(
            name="Ki 67-", index=0, previous_phase_index=1, next_phase_index=1,
            timing=TimingConfig(phase_duration=4.59 * 60, fixed_duration=False),
        )


class Ki67Positive(Phase):
    """

    Inherits :class:`Phase`. Defines Ki 67+ proliferating phase.

    This is a proliferating phenotype for cells that are replicating. Ki67 is a protein marker associated with
    proliferation. Transition to the next phase is set to be deterministic (the phase does use a fixed duration) by
    default.
    Default phase duration is 15.5h. By default, if no user defined custom entry function is defined (i.e.,
    `entry_function=None`), this phase will set its entry function to be :class:`Phase._double_target_volume`.
    By default, will set the volume change rates to be [change in volume]/[phase duration]. This phase does not calcify
    the cell.

    The parameters for this phase are based on the MCF-10A cell line
    https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-10a-cell-line
    https://www.ebi.ac.uk/ols/ontologies/bto/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FBTO_0001939

    """

    def __init__(self, config: PhaseConfig = None, dt: float = 0.1, time_unit: str = "min",
                 space_unit: str = "micrometer"):
        super().__init__(config, dt=dt, time_unit=time_unit, space_unit=space_unit)

        functions = self.config.functions
        if functions.entry_function is None:
            self.entry_function = self._double_target_volume
            self.entry_function_args = [None]

        if functions.exit_function is False:
            self.exit_function = None
            self.exit_function_args = [None]
        elif functions.exit_function is None:
            self.exit_function = self._halve_target_volume
            self.exit_function_args = [None]
        #TODO: fully implement the idea of doubling duration to derive rates
        #if cytoplasm_volume_change_rate is None and cytoplasm_doubling_duration is not None: #propose adding an argument to allow different phase times #that diverge from the total phase time. Do we want to build in safeties, i.e cytoplasm_doubling_duration>phasetime?
        #    cytoplasm_volume_change_rate = -np.log(0.05)/cytoplasm_doubling_duration
        #if nuclear_volume_change_rate is None and nuclear_doubling_duration is not None: #propose adding an argument to allow different phase times #that diverge from the total phase time. Do we want to build in safeties, i.e cytoplasm_doubling_duration>phasetime?
        #    nuclear_volume_change_rate = -np.log(0.05)/nuclear_doubling_duration
        #if fluid_change_rate is None and cytoplasm_doubling_duration is not None and nuclear_doubling_duration is not None:
        #    if cytoplasm_doubling_duration < nuclear_doubling_duration: #with custom rates we follow the assumption that the fluid intake rate is an order of magnitude faster then the fastest of the two rates
        #        fluid_change_rate = -np.log(0.05)/(cytoplasm_doubling_duration/10)
        #    else:
        #        fluid_change_rate = -np.log(0.05)/(nuclear_doubling_duration/10)
        #elif fluid_change_rate is None and cytoplasm_doubling_duration is not None:
        #    fluid_change_rate = -np.log(0.05)/(cytoplasm_doubling_duration/10) #with custom rates we follow the assumption that the fluid intake rate is an order of magnitude faster
        #elif fluid_change_rate is None and nuclear_doubling_duration is not None:
        #    fluid_change_rate = -np.log(0.05)/(cytoplasm_doubling_duration/10) #with custom rates we follow the assumption that the fluid intake rate is an order of magnitude faster

    @classmethod
    def default_config(cls):
        return PhaseConfig(
            name="Ki 67+", index=1, previous_phase_index=0, next_phase_index=0,
            timing=TimingConfig(phase_duration=15.5 * 60.0, fixed_duration=True),
            events=EventConfig(division_at_phase_exit=True, removal_at_phase_exit=False),
        )


class Ki67PositivePreMitotic(Ki67Positive):
    """

    Inherits :class:`Ki67Positive`. Defines Ki 67+ pre-mitotic proliferating phase. Only difference to
    :class:`Ki67Positive` is the phase length.

    This is a proliferating phenotype for cells that are replicating. Ki67 is a protein marker associated with
    proliferation. Transition to the next phase is set to be deterministic (the phase does use a fixed duration) by
    default.
    Default phase duration is 13h. By default, if no user defined custom entry function is defined (i.e.,
    `entry_function=None`), this phase will set its entry function to be :class:`Phase._double_target_volume`

    The parameters for this phase are based on the MCF-10A cell line
    https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-10a-cell-line
    https://www.ebi.ac.uk/ols/ontologies/bto/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FBTO_0001939

    """

    @classmethod
    def default_config(cls):
        return PhaseConfig(
            name="Ki 67+ pre-mitotic", index=1, previous_phase_index=0, next_phase_index=2,
            timing=TimingConfig(phase_duration=13.0 * 60.0, fixed_duration=True),
            events=EventConfig(division_at_phase_exit=True, removal_at_phase_exit=False),
        )


class Ki67PositivePostMitotic(Phase):
    """
    Inherits :class:`Phase`. Defines Ki 67+ post-mitotic phase, it represents the cell's reorganization.

    This is a rest phenotype for cells that are replicating. Ki67 is a protein marker associated with proliferation.
    Transition to the next phase is set to be deterministic (the phase does use a fixed duration) by default. Default
    phase duration is 2.5h. By default, if no user defined custom entry function is defined (i.e.,
    `entry_function=None`), this phase will set its entry function to be
    :class:`Ki67PositivePostMitotic._standard_Ki67_positive_postmit_entry_function`, which calls
    :class:`Phase._halve_target_volume`.

    The parameters for this phase are based on the MCF-10A cell line
    https://www.sciencedirect.com/topics/medicine-and-dentistry/mcf-10a-cell-line
    https://www.ebi.ac.uk/ols/ontologies/bto/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FBTO_0001939
    """

    def __init__(self, config: PhaseConfig = None, dt: float = 0.1, time_unit: str = "min",
                 space_unit: str = "micrometer"):
        super().__init__(config, dt=dt, time_unit=time_unit, space_unit=space_unit)

        if self.config.functions.entry_function is None:
            self.entry_function = self._standard_Ki67_positive_postmit_entry_function
            self.entry_function_args = [None]

    @classmethod
    def default_config(cls):
        return PhaseConfig(
            name="Ki 67+ post-mitotic", index=2, previous_phase_index=1, next_phase_index=0,
            timing=TimingConfig(phase_duration=2.5 * 60.0, fixed_duration=True),
            events=EventConfig(division_at_phase_exit=True, removal_at_phase_exit=False),
        )

    def _standard_Ki67_positive_postmit_entry_function(self, *args):
        """
        Calls :class:`Phase._halve_target_volume`.

        :param args: Not used
        :return:
        """
        self._halve_target_volume(*args)


class G0G1(Phase):
    """
    Inherits :class:`Phase`. Defines G0/G1 phase, it more representative of the quiescent phase than the first growth
    phase.

    This is a quiescent phenotype for cells that are replicating. Transition to the next phase is set to be stochastic
    (the phase does not use a fixed duration) by default. Default
    expected phase duration is 5.15h, the phase transition rate is, therefore, dt/5.15 1/h.
    This phase does not calcify the cell. Reference phase duration from https://www.ncbi.nlm.nih.gov/books/NBK9876/
    """

    def __init__(self, config: PhaseConfig = None, dt: float = 0.1, time_unit: str = "min",
                 space_unit: str = "micrometer"):
        super().__init__(config, dt=dt, time_unit=time_unit, space_unit=space_unit)

    @classmethod
    def default_config(cls):
        return PhaseConfig(
            name="G0/G1", index=0, previous_phase_index=2, next_phase_index=1,
            timing=TimingConfig(phase_duration=5.15 * 60.0, fixed_duration=False),
        )


class S(Phase):
    """
    Inherits :class:`Phase`. Defines S phase, it more representative of the growth phase than the inter-growth rest.

    This is a growth phenotype for cells that are replicating. Transition to the next phase is set to be stochastic
    (the phase does not use a fixed duration) by default. Default expected phase duration is 8h, the phase transition
    rate is, therefore, dt/8 1/h. By default, will set the volume change rates to be
    [change in volume]/[phase duration]. This phase does not calcify the cell. Reference phase duration from
    https://www.ncbi.nlm.nih.gov/books/NBK9876/
    """

    def __init__(self, config: PhaseConfig = None, dt: float = 0.1, time_unit: str = "min",
                 space_unit: str = "micrometer"):
        super().__init__(config, dt=dt, time_unit=time_unit, space_unit=space_unit)

        if self.config.functions.entry_function is None:
            self.entry_function = self._double_target_volume
            self.entry_function_args = [None]

        # without a set rate, the rates are derived from the initial volumes
        volume = self.config.volume
        rates = volume.rates
        phase_duration = self.phase_duration

        if rates.cytoplasm_volume_change_rate is None and volume.cytoplasm_fluid is not None and \
                volume.cytoplasm_solid is not None:
            self.cytoplasm_volume_change_rate = (volume.cytoplasm_fluid + volume.cytoplasm_solid) / (phase_duration / dt)

        elif rates.cytoplasm_volume_change_rate is None and volume.cytoplasm_fluid is not None:
            self.cytoplasm_volume_change_rate = volume.cytoplasm_fluid / (phase_duration / dt)

        elif rates.cytoplasm_volume_change_rate is None and volume.cytoplasm_solid is not None:
            self.cytoplasm_volume_change_rate = volume.cytoplasm_solid / (phase_duration / dt)

        elif rates.cytoplasm_volume_change_rate is None:
            self.cytoplasm_volume_change_rate = 1

        if rates.nuclear_volume_change_rate is None and volume.cytoplasm_fluid is not None and \
                volume.cytoplasm_solid is not None:
            self.nuclear_volume_change_rate = (volume.nuclear_fluid + volume.nuclear_solid) / (phase_duration / dt)

        elif rates.nuclear_volume_change_rate is None and volume.cytoplasm_fluid is not None:
            self.nuclear_volume_change_rate = volume.nuclear_fluid / (phase_duration / dt)

        elif rates.nuclear_volume_change_rate is None and volume.cytoplasm_solid is not None:
            self.nuclear_volume_change_rate = volume.nuclear_solid / (phase_duration / dt)

        elif rates.nuclear_volume_change_rate is None:
            self.nuclear_volume_change_rate = 1

        if rates.fluid_change_rate is None and volume.cytoplasm_fluid is not None and volume.nuclear_fluid is not None:
            self.fluid_change_rate = (volume.cytoplasm_fluid + volume.nuclear_fluid) / (phase_duration / dt)
        elif rates.fluid_change_rate is None and volume.cytoplasm_fluid is not None:
            self.fluid_change_rate = volume.cytoplasm_fluid / (phase_duration / dt)
        elif rates.fluid_change_rate is None and volume.nuclear_fluid is not None:
            self.fluid_change_rate = volume.nuclear_fluid / (phase_duration / dt)
        else:
            self.fluid_change_rate = 1

    @classmethod
    def default_config(cls):
        return PhaseConfig(
            name="S", index=1, previous_phase_index=0, next_phase_index=2,
            timing=TimingConfig(phase_duration=8 * 60.0, fixed_duration=False),
        )


class G2M(Phase):
    """
    Inherits :class:`Phase`. Defines G2M phase, it more representative of the mitosis phase than the growth.

    This is a growth phenotype for cells that are replicating. Transition to the next phase is set to be stochastic
    (the phase does not use a fixed duration) by default. Default expected phase duration is 5h, the phase transition
    rate is, therefore, dt/5 1/h. By default, will set the volume change rates to be
    [change in volume]/[phase duration]. This phase does not calcify the cell. Reference phase duration from
    https://www.ncbi.nlm.nih.gov/books/NBK9876/
    """

    def __init__(self, config: PhaseConfig = None, dt: float = 0.1, time_unit: str = "min",
                 space_unit: str = "micrometer"):
        super().__init__(config, dt=dt, time_unit=time_unit, space_unit=space_unit)

        if self.config.functions.exit_function is None:
            self.exit_function = self._halve_target_volume
            self.exit_function_args = [None]

    @classmethod
    def default_config(cls):
        return PhaseConfig(
            name="G2/M", index=2, previous_phase_index=1, next_phase_index=0,
            timing=TimingConfig(phase_duration=5 * 60.0, fixed_duration=False),
            events=EventConfig(division_at_phase_exit=True, removal_at_phase_exit=False),
        )


class Apoptosis(Phase):
    """
    Inherits :class:`Phase`. Defines apoptotic phenotype phase.

    This is a dead cell phenotype phase, the cell will shrink itself and should be removed from the simulation when this
    phase ends. Transition to the next phase is set to be deterministic (the phase does use a fixed duration) by
    default. Default phase duration is 8.6h. By default, if no custom user defined entry function is used (i.e.,
    `entry_function=None`), entry function is set to :class:`Apoptosis._standard_apoptosis_entry`.
    :class:`Apoptosis._standard_apoptosis_entry` sets all the cell target volumes from :class:`phenocellpy.cell_volume`
    to 0. The default mass change rates are `cytoplasm_volume_change_rate = 1/60` [volume/min],
    `nuclear_volume_change_rate = 0.35 / 60` [volume/min], `fluid_change_rate = 3 / 60`. This phase does not calcify
    the cell.
    """

    def __init__(self, config: PhaseConfig = None, dt: float = 0.1, time_unit: str = "min",
                 space_unit: str = "micrometer"):
        super().__init__(config, dt=dt, time_unit=time_unit, space_unit=space_unit)

        if self.config.functions.entry_function is None:
            self.entry_function = self._standard_apoptosis_entry
            self.entry_function_args = [None]

        self.relative_rupture_volume = self.config.volume.relative_rupture_volume

    @classmethod
    def default_config(cls):
        return PhaseConfig(
            name="Apoptosis", index=0, previous_phase_index=0, next_phase_index=0,
            timing=TimingConfig(phase_duration=8.6 * 60.0, fixed_duration=True),
            events=EventConfig(division_at_phase_exit=False, removal_at_phase_exit=True),
            volume=VolumeConfig(relative_rupture_volume=2,
                                rates=VolumeRatesConfig(cytoplasm_volume_change_rate=1 / 60,
                                                        nuclear_volume_change_rate=0.35 / 60,
                                                        fluid_change_rate=3 / 60, calcification_rate=0)),
        )

    def _standard_apoptosis_entry(self, *none):
        """
        Zeroes all the cell's target volumes. Keeps the nuclear to cytoplasm ratio the same.

        :param none: Not used. This is a custom entry function, therefore it has to have args
        :return:
        """

        # shrink cell
        self.volume.target_fluid_fraction = 0
        self.volume.cytoplasm_solid_target = 0
        self.volume.nuclear_solid_target = 0


class NecrosisSwell(Phase):
    """
    Inherits :class:`Phase`. Swelling part of the necrosis process.

    Represents the osmotic swell a necrotic cell goes through. By default, this phase uses a custom transition function
    (i.e., `check_transition_to_next_phase_functions=None`), it can be overwritten by a user defined one. The custom
    transition
    function is :class:`NecrosisSwell._necrosis_transition_function`, it returns true when the cell becomes bigger than
    its rupture volume. The default relative rupture volume is 2, i.e., the cell ruptures after doubling in volume.
    By default, if no custom user defined entry function is used (i.e., `entry_function=None`), entry function is set
    to :class:`NecrosisSwell._standard_necrosis_entry_function`. It zeroes the solid target volumes and the target
    cytoplasm to nuclear ratio, and sets the target fluid fraction to 1. This causes the cell to increase its volume.
    The default volume change rates are `cytoplasm_volume_change_rate = 0.0032 / 60.0`,
    `nuclear_volume_change_rate = 0.013 / 60.0`, `fluid_change_rate = 0.67 / 60.0`,
    `calcification_rate = 0.0042 / 60.0`. This phase does calcify the cell.
    """

    def __init__(self, config: PhaseConfig = None, dt: float = 0.1, time_unit: str = "min",
                 space_unit: str = "micrometer"):
        super().__init__(config, dt=dt, time_unit=time_unit, space_unit=space_unit)

        functions = self.config.functions
        if functions.entry_function is None:
            self.entry_function = self._standard_necrosis_entry_function
            self.entry_function_args = [None]

        if functions.check_transition_to_next_phase_function is None:
            self.check_transition_to_next_phase_function = self._necrosis_transition_function
            self.check_transition_to_next_phase_function_args = [None]

    @classmethod
    def default_config(cls):
        # this phase uses by default a custom transition check, so the duration is here to avoid issues
        return PhaseConfig(
            name="Necrotic (swelling)", index=0, previous_phase_index=0, next_phase_index=1,
            timing=TimingConfig(phase_duration=9e99, fixed_duration=False),
            volume=VolumeConfig(relative_rupture_volume=2,
                                rates=VolumeRatesConfig(cytoplasm_volume_change_rate=0.0032 / 60.0,
                                                        nuclear_volume_change_rate=0.013 / 60.0,
                                                        fluid_change_rate=0.67 / 60.0,
                                                        calcification_rate=0.0042 / 60.0)),
        )

    def _standard_necrosis_entry_function(self, *none):
        """
        Responsible for causing the osmotic swell.

        Zeroes the solid target volumes, sets the target cytoplasm to nuclear
        ratio to 0, sets the target fluid fraction to 1, sets the rupture volume to be double the current total volume.
        :param none: Not used. This is a custom entry function, therefore it has to have args
        :return: No return
        """

        # the cell wants to degrade the solids and swell by osmosis
        self.volume.target_fluid_fraction = 1
        self.volume.nuclear_solid_target = 0
        self.volume.cytoplasm_solid_target = 0

        self.volume.target_cytoplasm_to_nuclear_ratio = 0

        # set rupture volume

        self.volume.rupture_volume = self.volume.relative_rupture_volume * self.volume.total

    def _necrosis_transition_function(self, *none):
        """
        Custom phase transition function. The simulated cell should only change phase once it bursts (i.e., when its
        volume is above the rupture volume), it cares not how long or how little time it takes to reach that state.
        :param none: Not used. This is a custom transition function, therefore it has to have args
        :return: Flag for phase transition
        :rtype: bool
        """
        return self.volume.total > self.volume.rupture_volume


class NecrosisLysed(Phase):
    """
    Inherits :class:`Phase`. Ruptured necrotic cell

    Represents the already ruptured necrotic cell. The modeler should find a way to represent this fragmentary state,
    either create several tiny cells from the cell that ruptured, or create a disconnected cell.  By default, the
    transition to the next phase is deterministic and the phase lasts 60 days. The simulated cell should shrink and
    disappear before then, this is a safeguard to remove the cell "by hand" if it hasn't. By default, if no custom user
    defined entry function is used (i.e., `entry_function=None`), entry function is set to
    :class:`NecrosisLysed._standard_lysis_entry_function`. It zeroes all target volumes from
    :class:`phenocellpy.cell_volume`. The default volume change rates are:
    `cytoplasm_volume_change_rate = 0.0032 / 60.0`, `nuclear_volume_change_rate = 0.013 / 60.0`,
    `fluid_change_rate = 0.050 / 60.0`, `calcification_rate = 0.0042 / 60.0`. This phase calcifies the cell.

    """

    def __init__(self, config: PhaseConfig = None, dt: float = 0.1, time_unit: str = "min",
                 space_unit: str = "micrometer"):
        super().__init__(config, dt=dt, time_unit=time_unit, space_unit=space_unit)

        if self.config.functions.entry_function is None:
            self.entry_function = self._standard_lysis_entry_function
            self.entry_function_args = [None]

    @classmethod
    def default_config(cls):
        # 60 days, the cell should disappear naturally before then, but if it hasn't we do it
        return PhaseConfig(
            name="Necrotic (lysed)", index=1, previous_phase_index=0, next_phase_index=-1,
            timing=TimingConfig(phase_duration=60 * 60 * 24, fixed_duration=True),
            events=EventConfig(division_at_phase_exit=False, removal_at_phase_exit=True),
            volume=VolumeConfig(relative_rupture_volume=9e99,
                                rates=VolumeRatesConfig(cytoplasm_volume_change_rate=0.0032 / 60.0,
                                                        nuclear_volume_change_rate=0.013 / 60.0,
                                                        fluid_change_rate=0.050 / 60.0,
                                                        calcification_rate=0.0042 / 60.0)),
        )

    def _standard_lysis_entry_function(self, *none):
        """
        Zeroes all the cell's target volumes. Also zeroes the nuclear to cytoplasm ratio.
        :param none: Not used. This is a custom entry function, therefore it has to have args
        :return:
        """
        self.volume.target_fluid_fraction = 0
        self.volume.nuclear_solid_target = 0
        self.volume.cytoplasm_solid_target = 0

        self.volume.target_cytoplasm_to_nuclear_ratio = 0

        # set rupture volume

        self.volume.rupture_volume = self.volume.relative_rupture_volume * self.volume.total

def get_phase_class(name):
    """
    Fetches a (uninitialized) phase class by its class name (e.g., "Ki67Positive")

    :param name: Name of the phase class
    :type name: str
    :return: A phase class
    :rtype: type
    """
    phase_class = globals().get(name)
    if not (isinstance(phase_class, type) and issubclass(phase_class, Phase)):
        raise ValueError(f"{name!r} is not a phase class of phenocellpy.phases")
    return phase_class


def main():
    return


if __name__ == '__main__':
    dt = 1
    phase = Phase(dt=dt)
    sen = SenescentPhase(dt=dt)
    cells = [type('', (), {})() for _ in range(10)]
    for c in cells:
        c.p = sen.copy()

    for c in cells:
        print(c.p.volume)

    phase = Phase(dt=dt)
    sen = SenescentPhase(dt=dt)
    ki67n = Ki67Negative(dt=dt)
    test_ki67p = Ki67Positive(dt=dt)
    ki67ppre = Ki67PositivePreMitotic(dt=dt)
    ki67ppos = Ki67PositivePostMitotic(dt=dt)
    g0g1 = G0G1(dt=dt)
    s = S(dt=dt)
    g2m = G2M(dt=dt)
    ap = Apoptosis(dt=dt)
    necsw = NecrosisSwell(dt=dt)
    necLys = NecrosisLysed(dt=dt)

    def grow_phase_transition(*args):
        return args[0] >= args[1] and args[2] > args[4]

    def double_target_volumes(self, *none):
        self.volume.nuclear_solid_target *= 2
        self.volume.cytoplasm_solid_target *= 2

    custom = Phase(
        PhaseConfig(
            name="custom", index=1, previous_phase_index=0, next_phase_index=2,
            timing=TimingConfig(phase_duration=120, fixed_duration=True),
            volume=VolumeConfig(simulated_cell_volume=1),
            functions=FunctionsConfig(entry_function=double_target_volumes, entry_function_args=[None],
                                      check_transition_to_next_phase_function=grow_phase_transition,
                                      check_transition_to_next_phase_function_args=[0, 9, 0, 9],
                                      user_phase_time_step_args=(None, )),
        ),
        dt=dt, time_unit="min", space_unit="micrometer"
    )

    # print(test_ki67p.index)
    for _ in range(1000):
        print(phase.name, phase.time_step_phase(), phase.volume.total)
        print(sen.name, sen.time_step_phase(), sen.volume.total)
        print(ki67n.name, ki67n.time_step_phase(), ki67n.volume.total)
        print(test_ki67p.name, test_ki67p.time_step_phase(), test_ki67p.volume.total)
        print(ki67ppre.name, ki67ppre.time_step_phase(), ki67ppre.volume.total)
        print(ki67ppos.name, ki67ppos.time_step_phase(), ki67ppos.volume.total)
        print(g0g1.name, g0g1.time_step_phase(), g0g1.volume.total)
        print(s.name, s.time_step_phase(), s.volume.total)
        print(g2m.name, g2m.time_step_phase(), g2m.volume.total)
        print(ap.name, ap.time_step_phase(), ap.volume.total)
        print(necsw.name, necsw.time_step_phase(), necsw.volume.total)
        print(necLys.name, necLys.time_step_phase(), necLys.volume.total)
        print(custom.name, custom.time_step_phase(), custom.volume.total)
