"""
Some helpful utilities for the chemistry package
"""
from compat24 import wraps
import warnings

def chemistry_deprecated(old_func, new_func):
    """ Decorator indicating that something is deprecated

    Parameters
    ----------
    old_func : str
        The name of the old function that used to be part of the API
    new_func : str
        The name of the new function to call instead

    Examples
    --------
    When new_func is a string, this can be used as an informative decorator.

    >>> @chemistry_deprecated("old", "new")
    ... def old(arg1, arg2, kwarg1=None, kwarg2=None):
    ...     pass
    ... 
    >>> old(1, 2, 3, 4)

    This can be used to adorn an existing function with no additional effects as
    shown below

    >>> def new(arg1, arg2, kwarg1=None, kwarg2=None):
    ...    pass
    >>> old = chemistry_deprecated('old', 'new')(new)
    """
    def outer(func):
        @wraps(func)
        def wrapped(*args, **kwargs):
            warnings.warn('%s is deprecated. Use the %s function instead' %
                          (func.__name__, new_func), DeprecationWarning)
            return func(*args, **kwargs)
        return wrapped
    return outer
