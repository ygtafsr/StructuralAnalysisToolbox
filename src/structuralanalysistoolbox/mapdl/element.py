
from dataclasses import dataclass, fields
from enum import Enum
from typing import Literal, ClassVar



##########################
## ELEMENT TECHNOLOGIES  
##########################

class Tech185(Enum): 
    full_integration = 0
    reduced_integration = 1
    enhanced_strain = 2
    simplified_enhanced_strain = 3

class Tech186(Enum): 
    reduced_integration = 0
    full_integration = 1

class Tech182(Enum): 
    full_integration = 0
    reduced_integration = 1
    enhanced_strain = 2
    simplified_enhanced_strain = 3

##########################
## ELEMENT FORMULATIONS
##########################

class Form185(Enum): 
    pure_displacement = 0
    mixed_up = 1

class Form186(Enum): 
    pure_displacement = 0
    mixed_up = 1

class Form187(Enum): 
    pure_displacement = 0
    mixed_up_cons_press = 1
    mixed_up_lin_press = 2

class Form182(Enum): 
    pure_displacement = 0
    mixed_up = 1

class Form183(Enum): 
    pure_displacement = 0
    mixed_up = 1

class Form285(Enum):
    mixed_up = 0
    pure_displacement = 1

###############################

class Behv(Enum):
    plane_stress = 0
    axisymmetric = 1
    plane_strain = 2
    plane_stress_with_thickness = 3
    generalized_plane_strain = 5
    axisymmetric_with_torsion = 6
    
class Layer(Enum): 
    non_layered = 0
    layered = 1

class Shape(Enum):
    quad_8 = 0
    tria_6 = 1

class PML(Enum) : 
    not_include = 0
    include = 1

class Steady(Enum):
    disabled = 0
    enabled = 1

class SurfOut(Enum):
    basic = 0
    extra = 4

@dataclass
class EType:
    name : str = ""
    midx : int = 0

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        lines = [f"{cls}("]
        # include class-level constant
        lines.append(f"  name = '{self.name}',")
        # include dataclass fields
        for f in fields(self):
            value = getattr(self, f.name)
            lines.append(f"  {f.name} = {value!r},")
        lines.append(")")
        return "\n".join(lines)

@dataclass
class Solid185(EType):
    """8-Node Structural Solid"""
    name: ClassVar[str] = "SOLID185"
    technology : Tech185 = Tech185.full_integration
    formulation : Form185 = Form185.pure_displacement
    layering : Layer = Layer.non_layered
    pml_absorbing : PML =  PML.not_include
    steady_state : Steady = Steady.disabled
    extra_surface_output : SurfOut = SurfOut.basic

    def __repr__(self):
        return super().__repr__()

@dataclass
class Solid186:
    """20-Node Structural Solid"""
    name: ClassVar[str] = "SOLID186"
    technology : Tech186 = Tech186.reduced_integration
    formulation : Form186 = Form186.pure_displacement
    layering : Layer = Layer.non_layered
    pml_absorbing : PML = PML.not_include
    steady_state : Steady = Steady.disabled
    extra_surface_output : SurfOut = SurfOut.basic

    def __repr__(self):
        return super().__repr__()

@dataclass
class Solid187:

    name: ClassVar[str] = "SOLID187"
    formulation : Form187 = Form187.pure_displacement
    pml_absorbing : PML = PML.not_include
    steady_state : Steady = Steady.disabled
    extra_surface_output : SurfOut = SurfOut.basic

    def __repr__(self):
        return super().__repr__()

@dataclass
class Solid285:

    name: ClassVar[str] = "SOLID285"
    formulation : Form285 = Form285.mixed_up
    pml_absorbing : PML = PML.not_include
    steady_state : Steady = Steady.disabled
    extra_surface_output : SurfOut = SurfOut.basic

    def __repr__(self):
        return super().__repr__()

@dataclass
class Plane182:

    name: ClassVar[str] = "PLANE182"
    technology : Tech182
    behaviour : Behv
    formulation : Form182 = Form182.pure_displacement
    pml_absorbing : PML = PML.not_include
    extra_surface_output : SurfOut = SurfOut.basic

    def __repr__(self):
        return super().__repr__()

@dataclass
class Plane183:

    name: ClassVar[str] = "PLANE183"
    shape : Shape
    behaviour : Behv
    formulation : Form183 = Form183.pure_displacement
    pml_absorbing : PML = PML.not_include
    extra_surface_output : SurfOut = SurfOut.basic

    def __repr__(self):
        return super().__repr__()

##################################
### CONTACT and TARGET Elements
##################################

class Conta174:
    pass

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
    185: Solid185,
    186: Solid186,
    187: Solid187,
    285: Solid285,
    182: Plane182,
    183: Plane183
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