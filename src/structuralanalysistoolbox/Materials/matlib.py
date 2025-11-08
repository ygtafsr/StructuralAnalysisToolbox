

from pathlib import Path
import json
from enum import Enum
from dataclasses import asdict 

'''CRUD operations for materialdatabase.json'''
DB_PATH = Path(r"C:\Users\yigit\StructuralAnalysisToolbox\src\structuralanalysistoolbox\materials\materialdatabase.json")

class filter(Enum):
    byName = 1

class MatDatabase:
    
    def __init__(self):
        self._database = _deserialize()
        self._materials = self._database["materials"]

    def save(self, material):
        if self._exist(material.name):
            self._update(material)
        else:
            self._create(material) 
    
    def delete(self, mat_name : str):
        del self._materials[mat_name] # delete from memory
        _serialize(self._database) # write all database from start

    def list(self):
        """List materials by name."""
        mat = ""
        for material in list(self._materials.keys()):
            mat += material + "\n"
        print(mat)

    def _create(self, material):
        """converts a material object to a dict, then saves to the database file"""
        self._materials[material.name] = _convert_to_dict(material)
        _serialize(self._database) # write all database from start

    def _read(self, mat_name : str):
        """Returns database object (data structure that represents material in the memory)"""
        if self._exist(mat_name):
            return self._materials[mat_name]

    def _update(self, material):
        self.delete(mat_name=material.name) # delete from memory
        self._create(material=material) # create new object
        _serialize(self._database) # write all database from start

    def _exist(self, mat_name : str) -> bool:
        if mat_name in self._materials.keys(): return True
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
        # iterate over material data in each material Model
        for data_name, data_value in asdict(mat_model).items():
            _material_model_dict[data_name] = data_value
        _library_material["material models"][mat_model_name] = _material_model_dict
    return _library_material
   

def _serialize(database): 
    """transform database object to json object """
    with open(DB_PATH, mode='w', encoding='utf-8') as file:
        json.dump(database, file, indent=2)

def _deserialize():
    """transform json object to database object (not yet to a material object)"""
    with open(DB_PATH, mode="r", encoding="utf-8") as file:
        _database = json.load(file)
    return _database


