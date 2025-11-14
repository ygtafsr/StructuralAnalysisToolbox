
from __future__ import annotations
from . import mapdl

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from structuralanalysistoolbox import model
    

_model : model = None

map ={
    "add_MPC" : mapdl._MPC184
}
