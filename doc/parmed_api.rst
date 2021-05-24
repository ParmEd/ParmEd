The API for ParmEd Command-line Utilities
=========================================

In some cases, ``parmed`` interpreter is insufficient to do what you need it
to do, or perhaps you want to use some of ``parmed``'s functionality in your own
Python script or program.

In that case, the ParmEd API is here! Every :class:`Action
<parmed.tools.actions.Action>` subclass is available for import from
:mod:`parmed.tools`. For a full listing of available actions, you are forward to
the auto-generated list in :doc:`the parmed page <parmed>`.

Importing the actions
---------------------

All actions can be imported from ``parmed.tools`` with the same camelCase shown
in the ``parmed`` page. Examples are shown below::

    >>> from parmed.tools import change, addLJType, changeRadii, tiMerge
    >>> from parmed.tools import setBond, deleteBond, addDihedral, addPDB

Using the actions
-----------------

Each of the actions described on the previous page can be invoked either on a
parm list (a few, like :doc:`interpolate`, actually *require* a parm list with
multiple topologies), or a raw :class:`AmberParm <parmed.amber.Amberparm>`
object. Note that some actions do not support all :class:`AmberParm
<parmed.amber.AmberParm>` subclasses (like those for the CHARMM or Amoeba
force fields).

The first step is to instantiate the action with the :class:`AmberParm
<parmed.amber.AmberParm>` instance as the first argument. The arguments that
each action requires in the ParmEd interpreter must then be given as separate
arguments to the constructor. Keyword arguments can optionally be given as
keyword arguments to the action constructor.

Examples of valid syntax for :doc:`addLJType` are shown below::

    from parmed.tools import addLJType
    from parmed.amber import AmberParm

    parm = AmberParm('trx.prmtop')

    # All arguments separate
    action = addLJType(parm, "@1", "radius", 1.5, "epsilon", 0.5)

    # Also equivalent; keyword arguments given as keywords
    action = addLJType(parm, "@1", radius=1.5, epsilon=0.5)

Note that simply *instantiating* the object does not do anything. You need to
explicitly *execute* it as well, as shown below (continuting from the example
above)::

    action.execute()

The informational message printed to the ParmEd output are available by casting
the action to a string::

    print(str(action))

    # Equivalent:
    print('%s' % action)

*Note*: A backwards-incompatible change was introduced after version 2.7.3 in
which including all arguments as a single string in the first argument was
supported. However, it was impossible to maintain this behavior *and* support
file name paths with whitespace in them. As a result, this
backwards-incompatible change was made.
