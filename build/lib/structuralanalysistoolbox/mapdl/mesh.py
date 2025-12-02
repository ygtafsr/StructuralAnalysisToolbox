
from __future__ import annotations
from dataclasses import dataclass
from functools import cached_property

from structuralanalysistoolbox.mapdl import element

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from structuralanalysistoolbox import model
    from ansys.mapdl import reader as pymapdl_reader
    import pyvista as pv

class Node:
    def __init__(self, model : model.Model):

        self.id : int
        self.dof : int
        self.mapdl = model.mapdl

    def _create_free_node(self):
        """finds max node id, creates a new one with the next number."""
    
class Mesh():
    """
    Mesh Class to represent Finite Elements and Nodes.
    """
    def __init__(self, arv : pymapdl_reader.Archive):
        self.archieve = arv
        self.grid : pv.UnstructuredGrid = arv.grid
        self.nset_dict = arv.node_components
        self.elset_dict = arv.element_components

    def element_types(self):
        etypes = []
        for item in self.archieve.ekey:
            pyvista_etype_idx = item[0]
            pyvista_etype = item[1]
            etype = element.VALID_ETYPES[pyvista_etype]()
            etype.midx = pyvista_etype_idx
            etypes.append(etype)
        return etypes




    

