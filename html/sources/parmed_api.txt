The ParmEd API
==============

In some cases, ``parmed.py`` interpreter is insufficient to do what you need it
to do, or perhaps you want to use some of ParmEd's functionality in your own
Python script or program.

In that case, the ParmEd API is here! Every :class:`Action
<ParmedTools.ParmedActions.Action>` subclass is available for import from
:mod:`ParmedTools`. For a full listing of available actions, you are forward to
the auto-generated list in :doc:`the parmed.py page <parmed>`.

Importing the actions
---------------------

All actions can be imported from ``ParmedTools`` with the same camelcase shown
in the ``parmed.py`` page. Examples are shown below::

    >>> from ParmedTools import change, addLJType, changeRadii, tiMerge
    >>> from ParmedTools import setBond, deleteBond, addDihedral, addPDB

Using the actions
-----------------

Each of the actions described on the previous page can be invoked either on a
parm list (a few, like :doc:`interpolate`, actually *require* a parm list with
multiple topologies), or a raw :class:`AmberParm <chemistry.amber.Amberparm>`
object. Note that some actions do not support all :class:`AmberParm
<chemistry.amber.AmberParm>` subclasses (like those for the CHARMM or Amoeba
force fields).

The first step is to instantiate the action with the :class:`AmberParm
<chemistry.amber.AmberParm>` instance as the first argument. The arguments that
each action requires in the ParmEd interpreter can then be given as separate
arguments to the constructor or all as a single string in the second argument.
Keyword arguments can optionally be given as keyword arguments to the action
constructor.

Examples of valid syntax for :doc:`addLJType` are shown below::

    from ParmedTools import addLJType
    from chemistry.amber import AmberParm

    parm = AmberParm('trx.prmtop')

    # All arguments as one string, a la parmed.py
    action = addLJType(parm, "@1 radius 1.5 epsilon 0.5")

    # Equivalent; all arguments separate
    action = addLJType(parm, "@1", "radius", 1.5, "epsilon", 0.5)

    # Also equivalent; keyword arguments given as keywords
    action = addLJType(parm, "@1", radius=1.5, epsilon=0.5)

Note that simply *instantiating* the object does not do anything. You need to
explicitly *execute* it as well, as shown below (continuting from the example
above)::

    action.execute()

The informational message printed to the ParmEd output are available by casting
the action to a string::

    print str(action)

    # Equivalent:
    print '%s' % action
