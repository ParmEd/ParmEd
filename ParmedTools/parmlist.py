"""
List of topology file objects that can be edited in ParmEd. They can be indexed
via either the name of the topology file or by the order in which they were
loaded.
"""

from chemistry.amber import AmberFormat, AmberParm, AmoebaParm, ChamberParm
from ParmedTools.exceptions import DuplicateParm, ParmIndexError

class ParmList(object):
    """
    List of AmberParm objects index-able via either prmtop name or prmtop index
    (based on order added)
    """
    def __init__(self):
        """ Must be instantiated by itself """
        self._parm_names = []
        self._parm_instances = []
        self.parm = None # The active parm instance
        self.prev = None
        self.current_index = 0

    def set_new_active(self, idx):
        """ Sets a new active parm """
        self.current_index = self.index(idx)
        self.parm = self[self.current_index]

    def add_parm(self, parm, rst7=None):
        """ Add a parm to the list """
        # Make sure this parm is not part of the list already
        if str(parm) in self._parm_names:
            raise DuplicateParm('%s already in ParmList' % parm)
        # Convert a string to an AmberParm or add an AmberParm directly
        if not isinstance(parm, AmberFormat):
            parm = AmberFormat(parm)
            # From the parm data, we should be able to tell whether it was a
            # chamber topology or regular topology. Take the proper view and
            # add it to the list
            if 'CTITLE' in parm.flag_list:
                parm = parm.view(ChamberParm)
            elif 'AMOEBA_FORCEFIELD' in parm.flag_list:
                parm = parm.view(AmoebaParm)
            else:
                parm = parm.view(AmberParm)
        # Otherwise, add in the new parm's name
        self._parm_names.append(str(parm))
        self._parm_instances.append(parm)
        # A newly added topology file is the currently active parm
        self.current_index = len(self._parm_instances) - 1
        self.parm = parm
        # If we have a restart file, load it
        if rst7 is not None:
            parm.load_rst7(rst7)

    def index(self, idx):
        """ See what the index of the requested parm is (can be int or str) """
        if isinstance(idx, int):
            if idx < 0 or idx >= len(self._parm_instances):
                raise ParmIndexError('index out of range for ParmList')
            return idx
        else:
            try:
                return self._parm_names.index(str(idx))
            except ValueError:
                raise ParmIndexError('%s prmtop not in ParmList' % idx)

    def __getitem__(self, idx):
        """ Allow an integer index or string index to identify a parm """
        return self._parm_instances[self.index(idx)]

    def __contains__(self, name):
        """ See if a given parm name is in the list """
        if isinstance(name, int):
            return name >= 0 and name < len(self)
        else:
            return str(name) in self._parm_names
   
    def __iter__(self):
        """ Iterate through the prmtop names """
        return iter(self._parm_instances)

    def __len__(self):
        return len(self._parm_instances)

    def empty(self):
        return len(self._parm_instances) == 0
