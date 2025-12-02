
import json
from enum import Enum
from dataclasses import asdict 
from ..config.settings import DEFAULT_MAT_DATABASE_PATH

class filter(Enum):
    byName = 1

def _serialize(database): 
    """transform database object to json object """
    with open(DEFAULT_MAT_DATABASE_PATH, mode='w', encoding='utf-8') as file:
        json.dump(database, file, indent=2)

def _deserialize():
    """transform json object to database object (not yet to a material object)"""
    with open(DEFAULT_MAT_DATABASE_PATH, mode="r", encoding="utf-8") as file:
        _database = json.load(file)
    return _database

_database = _deserialize()
_materials = _database["materials"]
    
def delete(mat_name : str):
    del _materials[mat_name] # delete from memory
    _serialize(_database) # write all database from start

def list_materials():
    """List materials by name."""
    mat = ""
    for material in list(_materials.keys()):
        mat += material + "\n"
    return print(mat)

def _create(self, material):
    """converts a material object to a dict, then saves to the database file"""
    _materials[material.name] = _convert_to_dict(material)
    _serialize(_database) # write all database from start

def _read(mat_name : str):
    """Returns database object (data structure that represents material in the memory)"""
    if exist(mat_name):
        return _materials[mat_name]

def _update( material):
    delete(mat_name=material.name) # delete from memory
    _create(material=material) # create new object
    _serialize(_database) # write all database from start ????

def exist( mat_name : str) -> bool:
    if mat_name in _materials.keys(): return True
    else: return False

def _convert_to_dict(material) -> dict:
    """Converts material to a dict (database object)"""
    _library_material = {}
    # material labels
    _library_material["category"] = material.category
    _library_material["remarks"] = material.remarks
    _library_material["source"] = material.source
    _library_material["last update"] = material.last_update
    _library_material["material models"] = {}
    # iterate over material models in material object
    for mat_model_name, mat_model in material.material_models.items():
        # Create a dict for each material model (pyhsical, elasticity, plasticity etc.)
        _material_model_dict = {} 
        _material_model_dict["model name"] = type(mat_model).__name__
        # iterate over material label in each material Model
        for data_label, data_value in asdict(mat_model).items():
            if type(data_value) is list:
                _material_model_dict[data_label] = data_value
            else:
                _material_model_dict[data_label] = data_value
        _library_material["material models"][mat_model_name] = _material_model_dict
    return _library_material
   





