""" Various utilities used by ParmEd that don't really fit elsewhere """

__all__ = ['six', 'io', 'timer', 'which']

def which(prog):
    """ Returns the full path of a program if it exists in PATH

    Parameters
    ----------
    prog : str
        The name of a program to try and locate in PATH

    Returns
    -------
    path : str or None
        The full path of the program. If it cannot be found, None
    """
    import os
    def is_exe(fpath):
        if os.path.isdir(fpath): return False
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    fpath, fprog = os.path.split(prog)
    if fpath:
        if is_exe(prog):
            return prog
        return None
    for fpath in os.environ['PATH'].split(os.pathsep):
        trial = os.path.join(fpath, prog)
        if is_exe(trial):
            return trial
    return None
