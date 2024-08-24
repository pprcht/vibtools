# pyirtools/__init__.py
  
from ._irtools import __version__
from .calculator import IRtoolsCalculator, matchscore
from .readers import *
from .printouts import *

# Dynamically add all functions from the readers module to __all__
import inspect
from . import readers
from . import printouts

__all__ = ["__version__", "IRtoolsCalculator", "matchscore"]

# Add all functions from readers module to __all__
__all__.extend(name for name, obj in inspect.getmembers(readers) if inspect.isfunction(obj))
__all__.extend(name for name, obj in inspect.getmembers(printouts) if inspect.isfunction(obj))

