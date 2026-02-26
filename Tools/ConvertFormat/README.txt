# Converting a grid file from the FullSWOF_2D format to a GIS format (AscGrid)

xyz2asc.c is a C program converting a file in XYZ format (used by FullSWOF_2D)
into the AscGrid format (ASCII interchange format for Arc/Info Grid).
The AscGrid format should be suitable for your GIS software.
More details about the AscGrid format here: <https://www.loc.gov/preservation/digital/formats/fdd/fdd000421.shtml>

To get a description of the command-line inputs, simply launch the executable.

# Converting a grid file from a GIS format (AscGrid) to the FullSWOF_2D format

asc2xyz.c is a C program converting a file in AscGrid format (ASCII interchange format for Arc/Info Grid)
into the XYZ format used by FullSWOF_2D.
The AscGrid format should be suitable for your GIS software.
More details about the AscGrid format here: <https://www.loc.gov/preservation/digital/formats/fdd/fdd000421.shtml>

To get a description of the command-line inputs, simply launch the executable.
