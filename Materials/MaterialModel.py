
# Material Model Base File

## MP: Defines a linear material property as a constant or a function of temperature.
# Material Property Label (lab)

## MPDATA: Defines property data associated with the temperature table

## MPTEMP: Defines a temperature table for material properties.

## MPLIST: Lists linear material properties.

# *****************************************************************************************

## TB (Create Material Data Table): Activates a data table for material properties or special element input.

## TBDATA: Defines data for the material data table.

## TBFIELD

## TBTEMP

# *****************************************************************************************

## MPWRITE: Writes linear material properties in the database to a file (if the LIB option is not specified)
#  or writes both linear and nonlinear material properties (if LIB is specified) from the database to a file.

## MPREAD: Reads a file containing material properties.

## MPLIB: Sets the default material library read and write paths.

