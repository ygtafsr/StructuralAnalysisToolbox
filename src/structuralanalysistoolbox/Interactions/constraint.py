
from ansys.mapdl.core import Mapdl

# parameters are components retrieved through 'mapdl.components[comp_name]'
# components type must be 'node'

def create_RBE2(mapdl_instance : Mapdl, dependent_nodes : str, independent_node : str) -> None:
    
    # Check current processor ???
    # Check whether components are exist ???
    mapdl_instance.finish()
    mapdl_instance.prep7()

    # get the components content
    _dependent_nodes = mapdl_instance.components[dependent_nodes].items
    _independent_node = mapdl_instance.components[independent_node].items[0]

    _current_max_type_number = int(mapdl_instance.get('max_type_number','ETYP','','NUM','MAX'))
 
    # create the MPC184 elements
    mapdl_instance.et(_current_max_type_number+1,'MPC184',1,1)
    mapdl_instance.type(_current_max_type_number+1)

    for node in _dependent_nodes:
        mapdl_instance.e(node, _independent_node)

    mapdl_instance.allsel()