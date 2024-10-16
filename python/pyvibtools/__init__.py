# pyvibtools/__init__.py
  
from ._vibtools import __version__
from .calculator import vibtoolsCalculator, matchscore
from .readers import *
from .printouts import *
from .multiplot import *

# Dynamically add all functions from the readers module to __all__
import inspect
from . import readers
from . import printouts
from . import multiplot

__all__ = ["__version__", "vibtoolsCalculator", "matchscore"]

# Add all functions from readers module to __all__
__all__.extend(name for name, obj in inspect.getmembers(readers) if inspect.isfunction(obj))
__all__.extend(name for name, obj in inspect.getmembers(printouts) if inspect.isfunction(obj))
__all__.extend(name for name, obj in inspect.getmembers(multiplot) if inspect.isfunction(obj))
