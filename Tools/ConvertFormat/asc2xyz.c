/**
 * @file asc2xyz.c
 * @author Frederic Darboux <frederic.darboux@inra.fr> (2017)
 * @version 0.00.04
 * @date 2017-12-08
 *
 * @brief ArcGIS ASCII --> FullSWOF_2D XYZ
 * @details
 * Converts an ArcGIS ASCII file to a XYZ file for FullSWOF_2D.
 *
 * @copyright License Cecill-V2 \n
 * <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>
 *
 * (c) CNRS - Universite d'Orleans - BRGM - INRAE
 */
/*
 *
 * This file is part of FullSWOF_2D software.
 * <https://sourcesup.renater.fr/projects/fullswof-2d/>
 *
 * FullSWOF_2D = Full Shallow-Water equations for Overland Flow,
 * in two dimensions of space.
 * This software is a computer program whose purpose is to compute
 * solutions for 2D Shallow-Water equations.
 *
 * LICENSE
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * <http://www.cecill.info>.
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading, using, modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean that it is complicated to manipulate, and that also
 * therefore means that it is reserved for developers and experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and, more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *
 ******************************************************************************/


/*TODO: make header reading more flexible*/
/*TODO: read line by line to improve error checking*/

/*CHANGES:*/
/*2017-12-08 Version 0.00.04: include ctype.h  file to use function tolower()*/
/*2017-10-20 Version 0.00.03: Improved programming syntax*/
/*2017-08-02 Version 0.00.02: Improved check for header. Improved output header*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>

#define ASC_HEADER_NBLINE 6
#define TRUE 1
#define FALSE 0

#define LENGTH_OF_LINE 120

typedef float OUTDATA_TYPE;

static char *version = "asc2xyz - 0.00.04 - 2017-12-08";

int strcicmp(char const *a, char const *b);

int main(int argc, char *argv[]){
  char *filename_in, *filename_out;
  FILE *In, *Out;
  char ascheader[ASC_HEADER_NBLINE][LENGTH_OF_LINE], varName[LENGTH_OF_LINE];
  char nbcolFound, nbrowFound, cellsizeFound;
  float x, y, varValue;
  long int i, j;
  OUTDATA_TYPE *z;
  long int nbcol, nbrow, nbcells, startposition, ascHeaderNbLine;
  float cellsize;

  if(argc != 3){
    fprintf(stderr, "\t%s\n", version);
    fprintf(stderr, "\tUsage: %s input.asc output.xyz \n", argv[0]);
    exit(EXIT_FAILURE);
  }/*if*/

  filename_in  = malloc(strlen(argv[1])+1);
  if(!filename_in){
    fprintf(stderr,"asc2xyz: ERROR: data allocation failure for input file\n");
    exit(EXIT_FAILURE);
  }/*if*/
  
  filename_out = malloc(strlen(argv[2])+1);
  if(!filename_out){
    fprintf(stderr,"asc2xyz: ERROR: data allocation failure for output file\n");
    exit(EXIT_FAILURE);
  }/*if*/

  strcpy(filename_in,  argv[1]);
  strcpy(filename_out, argv[2]);

  In=fopen(filename_in, "r");
  if(In==NULL){
    fprintf(stderr,"\nasc2xyz: problem opening the input file\n");
    exit(EXIT_FAILURE);
  }/*if*/

  Out=fopen(filename_out, "w");
  if(Out==NULL){
    fprintf(stderr,"\nasc2xyz: problem opening the output file\n");
    exit(EXIT_FAILURE);
  }/*if*/

  /*read header of asc file*/
  for(i=0; i<ASC_HEADER_NBLINE; i++){
    if(fgets(ascheader[i], LENGTH_OF_LINE, In) == NULL){
      fprintf(stderr,"\nasc2xyz: problem reading the header\n");
      exit(EXIT_FAILURE);
    }
  }

  /*From header to parameters*/
  ascHeaderNbLine=0;
  nbcolFound=FALSE;
  nbrowFound=FALSE;
  cellsizeFound=FALSE;

  for(i=0; i<ASC_HEADER_NBLINE; i++){
    if(sscanf(ascheader[i], "%s %f\n", varName, &varValue) == 2){
      if(strcicmp(varName, "ncols")==0){
        nbcol=(int)varValue;
        ascHeaderNbLine++;
        nbcolFound=TRUE;
      }else if(strcicmp(varName, "nrows")==0){
        nbrow=(int)varValue;
        ascHeaderNbLine++;
        nbrowFound=TRUE;
      }else if(strcicmp(varName, "xllcorner")==0){
        ascHeaderNbLine++;
      }else if(strcicmp(varName, "yllcorner")==0){
        ascHeaderNbLine++;
      }else if(strcicmp(varName, "cellsize")==0){
        cellsize=(float)varValue;
        ascHeaderNbLine++;
        cellsizeFound=TRUE;
      }else if(strcicmp(varName, "NODATA_value")==0){
        ascHeaderNbLine++;
      }else{
        fprintf(stderr,"\nasc2xyz: ERROR: Unknown parameter '%s' in the header\n", varName);
        exit(EXIT_FAILURE);
      }/*else*/
    }/*if*/
    else{/*a line of the header has an unknown name*/
	    fprintf(stderr,"\nasc2xyz: ERROR: problem while reading the AscGrid header\n");
      exit(EXIT_FAILURE);
    }/*else*/
  }/*for*/
  if(ascHeaderNbLine != ASC_HEADER_NBLINE){/*Note: not sure this statement is useful*/
    fprintf(stderr,"\nasc2xyz: ERROR: wrong number of lines in the AscGrid header\n");
    exit(EXIT_FAILURE);
  }/*if*/

  if(nbcolFound == FALSE){
    fprintf(stderr,"\nasc2xyz: ERROR: parameter 'ncols' not found in the AscGrid header\n");
    exit(EXIT_FAILURE);
  }/*if*/
  if(nbrowFound == FALSE){
    fprintf(stderr,"\nasc2xyz: ERROR: parameter 'nrows' not found in the AscGrid header\n");
    exit(EXIT_FAILURE);
  }/*if*/
  if(cellsizeFound == FALSE){/*Note: not sure this statement is useful*/
    fprintf(stderr,"\nasc2xyz: ERROR: parameter 'cellsize' not found in the AscGrid header\n");
    exit(EXIT_FAILURE);
  }/*if*/

  /*Check parameter values*/
  if(nbcol<=0){
    fprintf(stderr,"\nasc2xyz: ERROR: problem reading the header (ncols <= 0)\n");
    exit(EXIT_FAILURE);
  }/*if*/
  if(nbrow<=0){
    fprintf(stderr,"\nasc2xyz: ERROR: problem reading the header (nrows <= 0)\n");
    exit(EXIT_FAILURE);
  }/*if*/
  if(cellsize<=0){
    fprintf(stderr,"\nasc2xyz: ERROR: problem reading the header (cellsize <= 0)\n");
    exit(EXIT_FAILURE);
  }/*if*/

  /*Alloc memory*/
  nbcells=nbcol*nbrow;

  z=malloc((size_t)((nbcells)*sizeof(OUTDATA_TYPE)));
  if(!z){
    fprintf(stderr,"\nasc2xyz: ERROR: data allocation failure\n");
    exit(EXIT_FAILURE);
  }/*if*/

  /*read z and store them in memory*/
  for(i=0;i<nbcells;i++){
    if(fscanf(In, "%f", &z[i]) != 1){
      fprintf(stderr,"\nasc2xyz: ERROR: problem reading the value of z in the input file (position=%li)\n", i);
      exit(EXIT_FAILURE);
    }/*if*/
  }/*for*/
  /*write header in the xyz*/
  fprintf(Out, "#x[i] y[j] z[i][j]\n");
  for(i=0; i<ASC_HEADER_NBLINE; i++){
    fprintf(Out,"#ASCGRID %s",ascheader[i]);
  }/*for*/

  /*write x, y and z in xyz*/
  startposition=nbcol*(nbrow-1);
  for(i=0;i<nbcol;i++){
    x=(float)i*cellsize+cellsize/2;
    for(j=0;j<nbrow;j++){
      y=(float)j*cellsize+cellsize/2;
      fprintf(Out,"%.20g %.20g %.6g\n", x, y, z[startposition-nbcol*j+i]); /* %.20g and %.6g is a trick to avoid trailing zeros*/
    }/*for*/
  }/*for*/


  free(z);
  free(filename_in);
  free(filename_out);
  fclose(In);
  fclose(Out);

  return EXIT_SUCCESS;
}/*main*/

int strcicmp(char const *a, char const *b) /*From https://stackoverflow.com/questions/5820810/case-insensitive-string-comp-in-c*/
{
  for (;; a++, b++) {
    int d = tolower(*a) - tolower(*b);
    if (d != 0 || !*a)
      return d;
  }/*for*/
}
