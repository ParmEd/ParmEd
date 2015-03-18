""" Exceptions used in Parmed """
from sys import stderr

class ParmError(Exception):
    """ Base parmed error """
    def __init__(self, msg='parmed error'):
        self.msg = msg
    def __str__(self):
        return self.msg

class ParmWarning(Warning, ParmError):
    """ Base parmed warning """
    def __init__(self, msg='parmed warning'):
        self.msg = msg
    def __str__(self):
        return self.msg

class SeriousParmWarning(ParmWarning):
    """ These warnings are more serious, and are fatal in strict operation """

class ChangeRadiiError(ParmError):
    pass

class WriteOFFError(ParmError):
    pass

class ParmedUtilsError(ParmError):
    pass

class ParmedChangeError(ParmError):
    pass

class ParmedAddLJTypeError(ParmError):
    pass

class ChangeLJPairError(ParmError):
    pass

class LJ_TypeError(ParmError):
    pass

class ParmedMoleculeError(ParmError):
    pass

class CoarseGrainError(ParmError):
    pass

class ChangeStateError(ParmError):
    pass

class SetParamError(ParmError):
    pass

class DeleteDihedralError(ParmError):
    pass

class ArgumentError(ParmError):
    pass

class NoArgument(ParmError):
    pass

class InterpreterError(ParmError):
    pass

class AmberIncompatibleWarning(ParmWarning):
    pass

class BadParmWarning(ParmWarning):
    pass

class FixableParmWarning(ParmWarning):
    pass

class NonfatalWarning(ParmWarning):
    pass

class NonUniversalWarning(ParmWarning):
    pass

class MissingDisulfide(ParmWarning):
    pass

class LongBondWarning(ParmWarning):
    pass

class NonexistentParm(ParmError):
    pass

class NonexistentParmWarning(SeriousParmWarning):
    pass

class DuplicateParm(ParmError):
    pass

class AmbiguousParmError(ParmError):
    pass

class IncompatibleParmsError(ParmError):
    pass

class AddPDBWarning(ParmWarning):
    pass

class AddPDBError(ParmError):
    pass

class LJ12_6_4Error(ParmError):
    pass

class DuplicateParamWarning(SeriousParmWarning):
    pass

class HMassRepartitionError(ParmError):
    pass

class SimulationError(ParmError):
    pass

class SimulationWarning(ParmWarning):
    pass

class UnhandledArgumentWarning(SeriousParmWarning):
    pass

class InputError(ParmError):
    pass

class ParmIndexError(ParmError, IndexError):
    pass

class FileExists(ParmError):
    pass

class FileDoesNotExist(ParmError, IOError):
    pass

class ChamberError(ParmError):
    pass

class TiMergeError(ParmError):
    pass

class WarningList(list):
    """ List of warnings """
   
    def __init__(self, empty_msg='No warnings found'):
        self._empty_msg = empty_msg
        list.__init__(self)

    def append(self, *args):
        raise NotImplemented('Use warn() to add to WarningList!')

    extend = append

    def warn(self, msg, exc_type=ParmWarning):
        """ Adds a warning to the list """
        list.append(self, (exc_type, msg))

    def dump(self, dest=stderr, ncols=80):
        """ Dump a list of all warnings to the destination """
        if len(self) == 0:
            dest.write(self._empty_msg + '\n')
            return

        dest.write('%d total warnings\n\n' % len(self))

        for w in self:
            words = ('%s: %s' % (w[0].__name__, w[1])).split()
            prstr = words[0] + ' '
            indent_chars = len(words[0]) + 1
            i = 1
            while i < len(words):
                if prstr and len(prstr) + len(words[i]) > ncols:
                    dest.write(prstr + '\n')
                    prstr = ' ' * indent_chars
                prstr += words[i] + ' '
                i += 1
            if prstr:
                dest.write(prstr + '\n')
            dest.write('\n')
#           dest.write('%s: %s\n' % (w[0].__name__, w[1]))
