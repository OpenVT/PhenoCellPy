try:
    from warnings import deprecated
except ImportError:  # Python < 3.13. Avoids depending on typing_extensions
    def deprecated(message):
        def decorator(func):
            return func
        return decorator

@deprecated("This function may become deprecated as we build a proper type system")
def _check_arguments(number_phases, phase_names, division_at_phase_exits, removal_at_phase_exits, fixed_durations,
                     phase_durations, entry_functions, entry_functions_args, exit_functions, exit_functions_args,
                     arrest_functions, arrest_functions_args, transitions_to_next_phase, transitions_to_next_phase_args,
                     cytoplasm_volume_change_rate, nuclear_volume_change_rate, calcification_rate, calcified_fraction,
                     target_fluid_fraction, nuclear_fluid, nuclear_solid, nuclear_solid_target, cytoplasm_fluid,
                     cytoplasm_solid, cytoplasm_solid_target, target_cytoplasm_to_nuclear_ratio, fluid_change_rate):
    """
    Checks that the numbers of parameters passed to Phenotype classes matches how many phases the phenotype has. E.g.,
    the Ki67Basic Phenotype has 2 phases, therefore it should receive 2 phase names, 2 flags for division at phase
    exit, etc.

    :param number_phases: How many phases compose the phenotype
    :type int
    :param phase_names: Names of the phases
    :type list
    :param division_at_phase_exits: Flags for division at phase exit
    :type list
    :param removal_at_phase_exits: Flags for removal at phase exits
    :type list
    :param fixed_durations: Flags for fixed duration phase
    :type list
    :param phase_durations: Time lengths of the phases
    :type list
    :param entry_functions: Functions that are executed upon entering the phase
    :type list
    :param entry_functions_args: List of lists of arguments for the entry functions
    :type list
    :param exit_functions: Functions that are executed as the phsae is exited
    :type list
    :param exit_functions_args: List of lists of arguments for the exit functions
    :type list
    :param arrest_functions: Functions defining exit from the phenotype and entrance to senescence
    :type list
    :param arrest_functions_args: List of lists of arguments for the arrest functions
    :type list
    :param transitions_to_next_phase: Functions for phase transition
    :type list
    :param transitions_to_next_phase_args: List of lists of arguments for the transition functions
    :type list
    :param cytoplasm_volume_change_rate: Change rates for cytoplasmic mass
    :type list
    :param nuclear_volume_change_rate: Change rates for nuclear mass
    :type list
    :param calcification_rate: Rates of calcification
    :type list
    :param calcified_fraction: Initially calcified fractions
    :type list
    :param target_fluid_fraction: Target fluid fractions for the phases
    :type list
    :param nuclear_fluid: Fluid nuclear volumes for the phases
    :type list
    :param nuclear_solid: Solid nuclear volumes for the phases
    :type list
    :param nuclear_solid_target: Target solid volumes for the phases
    :type list
    :param cytoplasm_fluid: Fluid cytoplasmic volume for each phase
    :type list
    :param cytoplasm_solid: Solid cytoplasmic volumes for each phase
    :type list
    :param cytoplasm_solid_target: Target solid cytoplasmic volumes for each phase
    :type list
    :param target_cytoplasm_to_nuclear_ratio: Target nuclear volume in relation to cytoplasmic volume
    :type list
    :param fluid_change_rate: Change rate for the fluid parts of the cell
    :type list
    :return: None
    """
    if len(division_at_phase_exits) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(division_at_phase_exits)} division flags defined")
    elif type(division_at_phase_exits) != list and type(division_at_phase_exits) != tuple:
        raise TypeError(f"`division_at_phase_exits` must be a list or tuple, got {type(division_at_phase_exits)}")

    if len(removal_at_phase_exits) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(removal_at_phase_exits)} removal flags defined")
    elif type(removal_at_phase_exits) != list and type(removal_at_phase_exits) != tuple:
        raise TypeError(f"`removal_at_phase_exits` must be a list or tuple, got {type(removal_at_phase_exits)}")

    if type(fixed_durations) != list and type(fixed_durations) != tuple:
        raise TypeError(f"`fixed_durations` must be a list or tuple, got {type(fixed_durations)}")
    elif len(fixed_durations) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(fixed_durations)} fixed duration flags defined")

    if len(phase_durations) != number_phases:
        raise ValueError(f"{phase_names} has {number_phases} phases, {len(phase_durations)} durations defined")
    elif type(phase_durations) != list and type(phase_durations) != tuple:
        raise TypeError(f"`phase_durations` must be a list or tuple, got {type(phase_durations)}")

    if len(entry_functions) != number_phases:
        raise ValueError(f"{phase_names} has {number_phases} phases, {len(entry_functions)} entry functions defined")
    elif type(entry_functions) != list and type(entry_functions) != tuple:
        raise TypeError(f"`entry_functions` must be a list or tuple, got {type(entry_functions)}")

    if len(entry_functions_args) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(entry_functions_args)} entry functions args defined")
    elif type(entry_functions_args) != list and type(entry_functions_args) != tuple:
        raise TypeError(f"`entry_functions_args` must be a list or tuple, got {type(entry_functions_args)}")
    #
    if len(exit_functions) != number_phases:
        raise ValueError(f"{phase_names} has {number_phases} phases, {len(exit_functions)} exit functions defined")
    elif type(exit_functions) != list and type(exit_functions) != tuple:
        raise TypeError(f"`exit_functions` must be a list or tuple, got {type(exit_functions)}")

    if len(exit_functions_args) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(exit_functions_args)} entry functions args defined")
    elif type(exit_functions_args) != list and type(exit_functions_args) != tuple:
        raise TypeError(f"`entry_functions_args` must be a list or tuple, got {type(exit_functions_args)}")
    #
    if len(arrest_functions) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(arrest_functions)} arrest functions defined")
    elif type(arrest_functions) != list and type(arrest_functions) != tuple:
        raise TypeError(f"`arrest_functions` must be a list or tuple, got {type(exit_functions)}")

    if len(arrest_functions_args) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(arrest_functions_args)} arrest functions args defined")
    elif type(arrest_functions_args) != list and type(arrest_functions_args) != tuple:
        raise TypeError(f"`arrest_functions_args` must be a list or tuple, got {type(arrest_functions_args)}")
    #
    if len(transitions_to_next_phase) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(transitions_to_next_phase)} transition functions defined")
    elif type(transitions_to_next_phase) != list and type(transitions_to_next_phase) != tuple:
        raise TypeError(
            f"`check_transition_to_next_phase_function` must be a list or tuple, got {type(transitions_to_next_phase)}")

    if len(transitions_to_next_phase_args) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(transitions_to_next_phase_args)} transition functions args "
            f"defined")
    elif type(transitions_to_next_phase_args) != list and type(transitions_to_next_phase_args) != tuple:
        raise TypeError(
            f"`check_transition_to_next_phase_functions_args` must be a list or tuple, got {type(arrest_functions_args)}")

    #
    if len(nuclear_volume_change_rate) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(nuclear_volume_change_rate)} nuclear volume change rates "
            f"defined")
    elif type(nuclear_volume_change_rate) != list and type(nuclear_volume_change_rate) != tuple:
        raise TypeError(
            f"`nuclear_volume_change_rate` must be a list or tuple, got {type(nuclear_volume_change_rate)}")

    #
    if len(calcification_rate) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(calcification_rate)} calcification rates "
            f"defined")
    elif type(calcification_rate) != list and type(calcification_rate) != tuple:
        raise TypeError(
            f"`calcification_rate` must be a list or tuple, got {type(calcification_rate)}")

    #
    if len(calcified_fraction) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(calcified_fraction)} calcified fractions defined")
    elif type(calcified_fraction) != list and type(calcified_fraction) != tuple:
        raise TypeError(
            f"`calcified_fraction` must be a list or tuple, got {type(calcified_fraction)}")
        #
    if len(cytoplasm_volume_change_rate) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(cytoplasm_volume_change_rate)} cytoplasm volume change "
            f"rates defined")
    elif type(cytoplasm_volume_change_rate) != list and type(cytoplasm_volume_change_rate) != tuple:
        raise TypeError(
            f"`calcified_fraction` must be a list or tuple, got {type(cytoplasm_volume_change_rate)}")

    #
    if len(target_fluid_fraction) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(target_fluid_fraction)} target fluid fractions defined")
    elif not (type(target_fluid_fraction) == list or type(target_fluid_fraction) == tuple):
        raise TypeError(
            f"`target_fluid_fraction` must be a list or tuple, got {type(target_fluid_fraction)}")
    #
    if len(nuclear_fluid) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(nuclear_fluid)} nuclear fluid volumes defined")
    elif not (type(nuclear_fluid) == list or type(nuclear_fluid) == tuple):
        raise TypeError(
            f"`nuclear_fluid` must be a list or tuple, got {type(nuclear_fluid)}")
    #
    if len(nuclear_solid) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(nuclear_solid)} nuclear solid volume defined")
    elif not (type(nuclear_solid) == list or type(nuclear_solid) == tuple):
        raise TypeError(
            f"`nuclear_solid` must be a list or tuple, got {type(nuclear_solid)}")
    #
    if len(nuclear_solid_target) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(nuclear_solid_target)} target target solid volumes defined")
    elif not (type(nuclear_solid_target) == list or type(nuclear_solid_target) == tuple):
        raise TypeError(
            f"`nuclear_solid_target` must be a list or tuple, got {type(nuclear_solid_target)}")

    #
    if len(cytoplasm_fluid) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(cytoplasm_fluid)} cytoplasm fluid volumes defined")
    elif not (type(cytoplasm_fluid) == list or type(cytoplasm_fluid) == tuple):
        raise TypeError(
            f"`cytoplasm_fluid` must be a list or tuple, got {type(cytoplasm_fluid)}")
    #
    if len(cytoplasm_solid) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(cytoplasm_solid)} cytoplasm solid volumes defined")
    elif not (type(cytoplasm_solid) == list or type(cytoplasm_solid) == tuple):
        raise TypeError(
            f"`cytoplasm_solid` must be a list or tuple, got {type(cytoplasm_solid)}")
    #
    if len(cytoplasm_solid_target) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(cytoplasm_solid_target)} cytoplasm target solid volumes "
            f"defined")
    elif not (type(cytoplasm_solid_target) == list or type(cytoplasm_solid_target) == tuple):
        raise TypeError(
            f"`cytoplasm_solid_target` must be a list or tuple, got {type(cytoplasm_solid_target)}")
    #
    if len(target_cytoplasm_to_nuclear_ratio) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(target_cytoplasm_to_nuclear_ratio)} target cytoplasm to "
            f"nuclear ratios defined")
    elif not (type(target_cytoplasm_to_nuclear_ratio) == list or type(target_cytoplasm_to_nuclear_ratio) == tuple):
        raise TypeError(
            f"`target_cytoplasm_to_nuclear_ratio` must be a list or tuple, got "
            f"{type(target_cytoplasm_to_nuclear_ratio)}")
    #
    if len(fluid_change_rate) != number_phases:
        raise ValueError(
            f"{phase_names} has {number_phases} phases, {len(fluid_change_rate)} fluid change rates defined")
    elif not (type(fluid_change_rate) == list or type(fluid_change_rate) == tuple):
        raise TypeError(
            f"`fluid_change_rate` must be a list or tuple, got {type(fluid_change_rate)}")
