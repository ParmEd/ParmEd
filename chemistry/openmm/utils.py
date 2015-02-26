"""
This package contains some useful functionality for common tasks in OpenMM
"""
from chemistry import unit as u

def energy_decomposition(structure, context):
    """
    This computes the energy of every force group in the given structure and
    computes the energy for each force group for the given Context.  Note, the
    context must have positions already assigned.

    Parameters
    ----------
    structure : Structure
        This should be the Structure object from which the System object in the
        Context was created
    context : mm.Context
        The OpenMM context set up for computing forces and energies

    Returns
    -------
    dict {str:float}
        A dictionary mapping the name of the force group (taken from the
        attribute names of the format XXX_FORCE_GROUP in the structure object)
        with the energy of that group in kcal/mol
    """
    all_names = dict()
    force_group_names = dict()
    energy_components = dict()
    kcal = u.kilocalories_per_mole

    for attr in dir(structure):
        if attr.endswith('_FORCE_GROUP'):
            val = getattr(structure, attr)
            all_names[val] = attr.replace('_FORCE_GROUP', '').lower()

    for force in context.getSystem().getForces():
        gp = force.getForceGroup()
        force_group_names[gp] = all_names[gp]

    for grp, name in force_group_names.iteritems():
        state = context.getState(getEnergy=True, groups=1<<grp)
        energy_components[name] = state.getPotentialEnergy().value_in_unit(kcal)

    e = context.getState(getEnergy=True).getPotentialEnergy()
    energy_components['total'] = e.value_in_unit(kcal)

    return energy_components
