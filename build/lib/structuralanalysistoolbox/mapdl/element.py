
from dataclasses import dataclass, field, fields
from abc import ABC, abstractmethod
from typing import Literal


@dataclass
class EType(ABC):
    midx : int = 0  # Element type index (Same as in solver)
    ignore_at_execution = False

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        lines = [f"{cls}("]
        # include dataclass fields
        for f in fields(self):
            value = getattr(self, f.name)
            lines.append(f"  {f.name} = {value!r},")
        lines.append(")")
        return "\n".join(lines)

SURF185_KEYOPT_MAPPING = {
    2 : ("technology", ("FULL INTEGRATION", "REDUCED INTEGRATION", "ENHANCED STRAIN", "SIMPLIFIED ENHANCED STRAIN")),
    3 : ("layer", ("LAYERED", "HOMOGENEOUS")),
    6 : ("formulation", ("PURE DISPLACEMENT", "MIXED UP")),
    15 : ("pml_absorbing", ("NOT INCLUDED", "INCLUDED")),
    16 : ("steady_state", ("DISABLED", "ENABLED")),
    17 : ("extra_surface_output", ("BASIC", "EXTRA")),
}

@dataclass
class Solid185(EType):
    """8-Node Structural Solid"""

    name = "SOLID185"

    technology : Literal["FULL INTEGRATION",
                         "REDUCED INTEGRATION",
                         "ENHANCED STRAIN",
                         "SIMPLIFIED ENHANCED STRAIN"] = "FULL INTEGRATION"

    formulation : Literal["PURE DISPLACEMENT",
                          "MIXED UP"] = "PURE DISPLACEMENT"
    
    layering : Literal["LAYERED",
                       "HOMOGENEOUS"] = "HOMOGENEOUS"

    pml_absorbing : Literal["NOT INCLUDED",
                            "INCLUDED"] = "NOT INCLUDED"
    
    steady_state : Literal["DISABLED",
                           "ENABLED"] = "DISABLED"      
    
    extra_surface_output : Literal["BASIC",
                                   "EXTRA"] = "BASIC"

    def get_keyopts(self) -> tuple: # ((keyopt_number, keyopt_value), ...)
        """Returns a tuple of keyopt number and corresponding value based on the current attribute settings."""

        return ((2 ,SURF185_KEYOPT_MAPPING[2][1].index(self.technology)),
                (3, SURF185_KEYOPT_MAPPING[3][1].index(self.layering)),
                (6, SURF185_KEYOPT_MAPPING[6][1].index(self.formulation)),
                (15, SURF185_KEYOPT_MAPPING[15][1].index(self.pml_absorbing)),
                (16, SURF185_KEYOPT_MAPPING[16][1].index(self.steady_state)),
                (17, SURF185_KEYOPT_MAPPING[17][1].index(self.extra_surface_output)))

    def set_keyopt(self, keyopt_number: int, value: str):
        """Sets the attribute corresponding to the keyopt_number to the given value."""

        # raise error if keyopt_number is invalid
        if keyopt_number not in SURF185_KEYOPT_MAPPING:
            raise ValueError(f"Invalid keyopt number: {keyopt_number}")
        # raise error if value is invalid
        if value not in SURF185_KEYOPT_MAPPING[keyopt_number][1]:
            raise ValueError(f"Invalid value for keyopt {keyopt_number}: {value}")  
        # set the attribute
        setattr(self, SURF185_KEYOPT_MAPPING[keyopt_number][0], 
                      SURF185_KEYOPT_MAPPING[keyopt_number][1][value])

    def __repr__(self):
        return super().__repr__()

@dataclass
class Solid186:
    """20-Node Structural Solid"""

    name = "SOLID186"

    technology : Literal["FULL INTEGRATION",
                         "REDUCED INTEGRATION"] = "REDUCED INTEGRATION"
    
    formulation : Literal["PURE DISPLACEMENT",
                          "MIXED UP"] = "PURE DISPLACEMENT"
    
    layering : Literal["LAYERED",
                       "HOMOGENEOUS"] = "HOMOGENEOUS"
    
    pml_absorbing : Literal["INCLUDED",
                            "NOT INCLUDED"] = "NOT INCLUDED"
    
    steady_state : Literal["ENABLED",
                          "DISABLED"] = "DISABLED"
    
    extra_surface_output : Literal["BASIC",
                                   "EXTRA"] = "BASIC"
    

    def __repr__(self):
        return super().__repr__()
    
@dataclass
class Solid187:

    name = "SOLID187"
    formulation : Literal["PURE DISPLACEMENT",
                          "MIXED UP"] = "MIXED UP"
    
    pml_absorbing : Literal["INCLUDED",
                            "NOT INCLUDED"] = "NOT INCLUDED"
    
    steady_state : Literal["ENABLED",
                          "DISABLED"] = "DISABLED"
    
    extra_surface_output : Literal["BASIC",
                                   "EXTRA"] = "BASIC"

    def __repr__(self):
        return super().__repr__()

@dataclass
class Solid285:

    name = "SOLID285"
    formulation : Literal["PURE DISPLACEMENT",
                          "MIXED UP WITH CONSTANT PRESSURE",
                          "MIXED UP WITH LINEAR PRESSURE"] = "PURE DISPLACEMENT"
    
    pml_absorbing : Literal["INCLUDED",
                            "NOT INCLUDED"] = "NOT INCLUDED"
    
    steady_state : Literal["ENABLED",
                          "DISABLED"] = "DISABLED"
    
    extra_surface_output : Literal["BASIC",
                                   "EXTRA"] = "BASIC"
    def __repr__(self):
        return super().__repr__()

@dataclass
class Plane182:

    name = "PLANE182"

    technology : Literal["FULL INTEGRATION",
                         "REDUCED INTEGRATION",
                         "ENHANCED STRAIN",
                         "SIMPLIFIED ENHANCED STRAIN"] = "FULL INTEGRATION"
    
    shape : Literal["QUAD 4",
                    "TRIA 3"] = "QUAD 4"

    behaviour : Literal["PLANE STRESS",
                       "AXISYMMETRIC",
                        "PLANE STRAIN",
                        "PLANE STRESS WITH THICKNESS",
                        "GENERALIZED PLANE STRAIN",
                        "AXISYMMETRIC WITH TORSION"] = "PLANE STRESS"
    
    formulation : Literal["PURE DISPLACEMENT",
                          "MIXED UP"] = "PURE DISPLACEMENT"
    
    pml_absorbing : Literal["INCLUDED",
                            "NOT INCLUDED"] = "NOT INCLUDED"
    
    extra_surface_output : Literal["BASIC",
                                   "EXTRA"] = "BASIC"

    def __repr__(self):
        return super().__repr__()

@dataclass
class Plane183:

    name = "PLANE183"

    shape : Literal["QUAD 8",
                    "TRIA 6"] = "QUAD 8"
    
    behaviour : Literal["PLANE STRESS",
                       "AXISYMMETRIC",
                        "PLANE STRAIN",
                        "PLANE STRESS WITH THICKNESS",
                        "GENERALIZED PLANE STRAIN",
                        "AXISYMMETRIC WITH TORSION"] = "PLANE STRESS"
    
    formulation : Literal["PURE DISPLACEMENT",
                          "MIXED UP"] = "PURE DISPLACEMENT"
    
    pml_absorbing : Literal["INCLUDED",
                            "NOT INCLUDED"] = "NOT INCLUDED"
    
    extra_surface_output : Literal["BASIC",
                                   "EXTRA"] = "BASIC"

    def __repr__(self):
        return super().__repr__()
    
##################################
### Surface Effect Elements
##################################

class Surf153:
    """2-D Structural Surface Effect"""
    pass

SURF154_KEYOPT_MAPPING = {
    "pressure_cs": (2, ("ELEMENT CS", "LOCAL CS")),
    "midside_nodes": (4, ("INCLUDE", "EXCLUDE")),
    "normal_pressure_direction": (6, ("CALCULATED", "POSITIVE ONLY", "NEGATIVE ONLY")),
    "large_deflection_area": (7, ("USE NEW AREA", "ORIGINAL AREA")),
    "ocean_loading": (8, ("IGNORE", "APPLY")),
    "press_vector_orient": (11, ("PRJ AREA + INCLUDE TANGENT",
                                 "PRJ AREA + EXCLUDE TANGENT",
                                 "FULL AREA + INCLUDE TANGENT")),
    "effect_of_direction": (12, ("PRESSURE APPLIED", "PRESSURE SUPPRESSED"))
}

class Surf154(EType):
    """3-D Structural Surface Effect"""
    name = "SURF154"
    pressure_cs : Literal["ELEMENT CS", "LOCAL CS"] = "ELEMENT CS"
    midside_nodes : Literal["INCLUDE", "EXCLUDE"] = "INCLUDE"
    normal_pressure_direction : Literal["CALCULATED", "POSITIVE ONLY", "NEGATIVE ONLY"] | None = None
    large_deflection_area : Literal["USE NEW AREA", "ORIGINAL AREA"] = "USE NEW AREA"
    ocean_loading : Literal["IGNORE", "APPLY"] = "IGNORE"
    press_vector_orient : Literal["PRJ AREA + INCLUDE TANGENT",
                                  "PRJ AREA + EXCLUDE TANGENT",
                                  "FULL AREA + INCLUDE TANGENT"] | None = None
    effect_of_direction : Literal["PRESSURE APPLIED", "PRESSURE SUPPRESSED"] | None = None
    
    def get_keyopts(self) -> list:
        keyopts = []
        for key, value in SURF154_KEYOPT_MAPPING.items():
            value = getattr(self, key)
            if value is not None: 
                keyopt_number, options = SURF154_KEYOPT_MAPPING[key]
                keyopt_value = options.index(value)
                keyopts.append((keyopt_number, keyopt_value))
        return keyopts
                    
    def __repr__(self):
        return super().__repr__()

##################################
### MPC184
##################################

MPC184_KEYOPT_MAPPING = {
    1 : ("behaviour", ("RigidLink", "RigidBeam", "Slider", "Revolute", "Universal", "Slot", "Point", "Translational", "Cylindrical", "Planar", "Weld", "Orient", "Spherical", "General", "Screw")),
    2 : ("method", ("Direct Elimination", "Lagrange Multiplier", "Penalty Based Method")),
    3 : ("configuration", ("x-axis revolute", "y-axis revolute", "displacement and rotational", "only displacement")),
    5 : ("stress_stiffness", ("Include", "Exclude"))}

class MPC184(EType):
    
    name = "MPC184"
    behaviour : int = 0
    method : Literal["Direct Elimination", "Lagrange Multiplier", "Penalty Based Method"] = "Lagrange Multiplier"
    contribution : Literal["Include", "Exclude"] = "Include"

    def get_keyopts(self) -> tuple:
        return ((1, self.behaviour),
                (2, ["Direct Elimination", "Lagrange Multiplier", "Penalty Based Method"].index(self.method)),
                (3, ["Include", "Exclude"].index(self.contribution)))

    def __repr__(self):
        return super().__repr__()

##################################
### CONTACT and TARGET Elements
##################################

CONTA174_KEYOPT_MAPPING = {}

class Conta174(EType):
    
    name = "CONTA174"
    real_constants : dict
    friction : float    # TB material property
    dof : int = 0   # keyopt1
    contact_algorithm : Literal["Augmented Lagrange", 
                                "Penalty", 
                                "MPC",
                                "Lagrange+Penalty", 
                                "Pure Lagrange"] = "Augmented Lagrange" # keyopt2
    
    units_of_normal_contact_stiffness : Literal["FORCE/LENGTH^3", "FORCE/LENGTH"] = "FORCE/LENGTH^3" # keyopt3

    location_of_contact_detection_point : Literal["On Gauss Point", 
                                                  "On Nodal Point-Normal From Contact Surface",
                                                  "On Nodal Point-Normal To Target Surface",
                                                  "On Nodal Point-Normal From Contact Surface",
                                                  "On Nodal Point-Normal From Contact Surface"] = "On Gauss Point" # keyopt4
    
    # To define surface based constraints :: For keyopt2 = MPC(2) or Lagrange_Multiplier(3)
    surface_based_constraint : Literal["Force Distribution",
                                       "Rigid Surface",
                                       "Coupling"]
    
    # result section :: How motion of the result section is accounted for in large deflection analysis
    result_section : Literal["Moves with Contact Elements",
                             "Moves with Underlying Elements"] = "Moves with Contact Elements" # keyopt4
    
    automated_adjustment : Literal["Close_Gap",
                                   "Reduce_Penetration",
                                   "Close Gap/Reduce Penetration"]
    



    
    



class Targe170:
    pass


########################################################
### SECTION Definitons for the elements with two points
########################################################

@dataclass
class SectionType:
    """
    SECTYPE, SECID, Type, Subtype, Name, REFINEKEY
    Associates section type information with a section ID number.
    https://ansyshelp.ansys.com/public/account/secured?returnurl=///Views/Secured/corp/v242/en/ans_cmd/Hlp_C_SECTYPE.html
    """
    id : int
    type : Literal[""]

@dataclass
class SectionData:
    """
    SECDATA, VAL1, VAL2, VAL3, VAL4, VAL5, VAL6, VAL7, VAL8, VAL9, VAL10, VAL11, VAL12
    Describes the geometry of a section.
    https://ansyshelp.ansys.com/public/account/secured?returnurl=///Views/Secured/corp/v242/en/ans_cmd/Hlp_C_SECDATA.html
    """
    pass

"""####################################
### MPC184 ELEMENT CONFIGURATION
####################################

MPCtypes = {
            "RigidLink" : 0,
            "RigidBeam" : 1,
            "Slider" : 3,
            "Revolute" : 6,
            "Universal" : 7,
            "Slot" : 8,
            "Point" : 9,
            "Translational" : 10,
            "Cylindrical" : 11,
            "Planar" : 12,
            "Weld" : 13,
            "Orient" : 14,
            "Spherical" : 15,
            "General" : 16,
            "Screw" : 17
            }

behaviour = Literal["RigidLink",
                    "RigidBeam",
                    "Slider",
                    "Revolute",
                    "Universal",
                    "Slot",
                    "Point",
                    "Translational",
                    "Cylindrical",
                    "Planar",
                    "Weld",
                    "Orient",
                    "Spherical",
                    "General",
                    "Screw"]
                    
                
method = Literal["Direct Elimination",
                "Lagrange Multiplier",
                "Penalty Based Method"]
                

config = Literal["x-axis revolute",
                 "y-axis revolute",
                 "displacement and rotational",
                 "only displacement"]

stress_stiffness = Literal["Include",
                           "Exclude"]"""

"""@dataclass
class MPC184:
    
    keyopt1 : behaviour # Element behavior
    keyopt2 : method # Element constraint imposition method
    keyopt4 : None | config # Element configuration
    keyopt5 : None | stress_stiffness  # Geometric stress-stiffness contributions for rigid beams"""


VALID_ETYPES = {
    154: Surf154,
    185: Solid185,
    186: Solid186,
    187: Solid187,
    285: Solid285,
    182: Plane182,
    183: Plane183
}

name : Literal["VERTEX", # PyVista Cell Types
                   "LINE",
                   "TRIA",
                   "QUAD",
                   "TETRA",
                   "HEX",
                   "WEDGE",
                   "PYRAMID",
                   "SURFACE"]
_etype_map = {
    "Surf154" : "SURFACE",
    "Solid185" : "HEX",
    "Solid186" : "HEX",
    "Solid187" : "HEX",
    "Solid285" : ""
}



"""
# element type to VTK conversion function call map
# 0: skip
# 1: Point
# 2: Line (linear or quadratic)
# 3: Shell
# 4: 3D Solid (Hexahedral, wedge, pyramid, tetrahedral)
# 5: Tetrahedral
# 6: Line (always linear)

_etype_map = [
    0,
    2,  # LINK1
    3,  # PLANE2
    2,  # BEAM3
    2,  # BEAM4
    4,  # SOLID5
    0,  # UNUSED6
    1,  # COMBIN7
    2,  # LINK8
    0,  # INFIN9
    2,  # LINK10
    2,  # LINK11
    2,  # CONTAC12
    3,  # PLANE13
    2,  # COMBIN14
    0,  # FLUID15
    2,  # PIPE16
    2,  # PIPE17
    2,  # PIPE18
    0,  # SURF19
    2,  # PIPE20
    1,  # MASS21
    3,  # SURF22
    2,  # BEAM23
    2,  # BEAM24
    3,  # PLANE25
    0,  # UNUSED26
    0,  # MATRIX27
    3,  # SHELL28
    3,  # FLUID29
    4,  # FLUID30
    2,  # LINK31
    2,  # LINK32
    2,  # LINK33
    2,  # LINK34
    3,  # PLANE35
    0,  # SOURC36
    2,  # COMBIN37
    2,  # FLUID38
    2,  # COMBIN39
    2,  # COMBIN40
    3,  # SHELL41
    3,  # PLANE42
    3,  # SHELL43
    2,  # BEAM44
    4,  # SOLID45
    4,  # SOLID46
    0,  # INFIN47
    0,  # UNUSED48
    0,  # UNUSED49
    0,  # MATRIX50
    3,  # SHELL51
    0,  # CONTAC52
    3,  # PLANE53
    3,  # BEAM54
    3,  # PLANE55
    0,  # UNUSED56
    3,  # SHELL57
    0,  # UNUSED58
    2,  # PIPE59
    2,  # PIPE60
    2,  # SHELL61
    4,  # SOLID62
    3,  # SHELL63
    4,  # SOLID64
    4,  # SOLID65
    2,  # FLUID66
    3,  # PLANE67
    2,  # LINK68
    4,  # SOLID69
    4,  # SOLID70
    1,  # MASS71
    0,  # UNUSED72
    0,  # UNUSED73
    0,  # UNUSED74
    3,  # PLANE75
    0,  # UNUSED76
    3,  # PLANE77
    3,  # PLANE78
    3,  # FLUID79
    4,  # FLUID80
    3,  # FLUID81
    3,  # PLANE82
    3,  # PLANE83
    0,  # UNUSED84
    0,  # UNUSED85
    0,  # UNUSED86
    5,  # SOLID87
    3,  # VISCO88
    4,  # VISCO89
    4,  # SOLID90
    3,  # SHELL91
    5,  # SOLID92
    3,  # SHELL93
    0,  # CIRCU94
    4,  # SOLID95
    4,  # SOLID96
    4,  # SOLID97
    5,  # SOLID98
    3,  # SHELL99
    4,  # USER100
    4,  # USER101
    4,  # USER102
    4,  # USER103
    4,  # USER104
    4,  # USER105
    3,  # VISCO106
    4,  # VISCO107
    4,  # VISCO108
    0,  # TRANS109
    0,  # INFIN110
    0,  # INFIN111
    0,  # ROT112
    0,  # UNUSED113
    0,  # UNUSED114
    3,  # INTER115
    2,  # FLUID116
    4,  # EDGE117
    3,  # HF118
    5,  # HF119
    4,  # HF120
    3,  # PLANE121
    4,  # SOLID122
    5,  # SOLID123
    0,  # CIRCU124
    0,  # CIRCU125
    2,  # TRANS126
    0,  # UNUSED127
    0,  # UNUSED128
    2,  # FLUID129
    3,  # FLUID130
    3,  # SHELL131
    3,  # SHELL132
    0,  # UNUSED133
    0,  # UNUSED134
    0,  # TRANS135
    3,  # FLUID136
    0,  # FLUID137
    0,  # FLUID138
    0,  # FLUID139
    5,  # ROM140
    0,  # UNUSED141
    0,  # UNUSED142
    3,  # SHELL143
    0,  # ROM144
    0,  # UNUSED145
    0,  # UNUSED146
    0,  # UNUSED147
    0,  # UNUSED148
    0,  # UNUSED149
    0,  # UNUSED150
    2,  # SURF151
    3,  # SURF152
    2,  # SURF153
    3,  # SURF154
    3,  # SURF155
    2,  # SURF156
    3,  # SHELL157
    0,  # UNUSED158
    0,  # SURF159
    0,  # UNUSED160
    2,  # BEAM161
    0,  # UNUSED162
    3,  # SHELL163
    4,  # SOLID164
    0,  # UNUSED165
    0,  # UNUSED166
    0,  # UNUSED167
    5,  # SOLID168
    0,  # TARGE169
    3,  # TARGE170
    2,  # CONTA171
    2,  # CONTA172
    3,  # CONTA173
    3,  # CONTA174
    1,  # CONTA175
    2,  # CONTA176
    2,  # CONTA177
    2,  # CONTA178
    0,  # PRETS179
    2,  # LINK180
    3,  # SHELL181
    3,  # PLANE182
    3,  # PLANE183
    0,  # RBAR184
    4,  # SOLID185
    4,  # SOLID186
    5,  # SOLID187
    6,  # BEAM188
    2,  # BEAM189
    4,  # SOLSH190
    0,  # UNUSED191
    4,  # INTER192
    0,  # INTER193
    0,  # INTER194
    0,  # INTER195
    0,  # LAYER196
    0,  # LAYER197
    0,  # UNUSED198
    0,  # UNUSED199
    0,  # UNUSED200
    0,  # FOLLW201
    0,  # INTER202
    0,  # INTER203
    0,  # INTER204
    0,  # INTER205
    0,  # INTER206
    0,  # UNUSED207
    2,  # SHELL208
    2,  # SHELL209
    0,  # UNUSED210
    0,  # UNUSED211
    3,  # CPT212
    3,  # CPT213
    6,  # COMBI214
    4,  # CPT215
    6,  # CPT216
    6,  # CPT217
    3,  # FLUID218
    3,  # FLUID219
    4,  # FLUID220
    5,  # FLUID221
    3,  # PLANE222
    3,  # PLANE223
    0,  # UNUSED224
    0,  # SOLID225
    4,  # SOLID226
    5,  # SOLID227
    0,  # UNUSED228
    0,  # UNUSED229
    3,  # PLANE230
    4,  # SOLID231
    5,  # SOLID232
    3,  # PLANE233
    0,  # UNUSED234
    0,  # UNUSED235
    4,  # SOLID236
    5,  # SOLID237
    3,  # PLANE238
    4,  # SOLID239
    5,  # SOLID240
    0,  # HSFLD241
    0,  # HSFLD242
    0,  # UNUSED243
    0,  # UNUSED244
    0,  # UNUSED245
    0,  # UNUSED246
    0,  # UNUSED247
    0,  # UNUSED248
    0,  # UNUSED249
    2,  # COMBI250
    2,  # SURF251
    3,  # SURF252
    0,  # UNUSED253
    0,  # UNUSED254
    0,  # UNUSED255
    0,  # UNUSED256
    0,  # INFIN257
    0,  # INFIN258
    0,  # INFIN259
    0,  # UNUSED260
    4,  # UNUSED261
    0,  # GLINK262
    0,  # REINF263
    0,  # REINF264
    0,  # REINF265
    0,  # UNUSED266
    0,  # UNUSED267
    0,  # UNUSED268
    0,  # UNUSED269
    0,  # UNUSED270
    0,  # UNUSED271
    0,  # SOLID272
    0,  # SOLID273
    0,  # UNUSED274
    0,  # UNUSED275
    0,  # UNUSED276
    0,  # UNUSED277
    4,  # SOLID278
    4,  # SOLID279
    2,  # CABLE280
    3,  # SHELL281
    3,  # SHELL282
    3,  # SHELL283
    0,  # UNUSED284
    5,  # SOLID285
    0,  # UNUSED286
    0,  # UNUSED287
    2,  # PIPE288
    2,  # PIPE289
    2,  # ELBOW290
    5,  # SOLID291
    3,  # PLANE292
    3,  # PLANE293
    0,  # UNUSED294
    0,  # UNUSED295
    0,  # UNUSED296
    0,  # UNUSED297
    0,  # UNUSED298
    0,  # UNUSED299
    0,  # USER300
]
"""