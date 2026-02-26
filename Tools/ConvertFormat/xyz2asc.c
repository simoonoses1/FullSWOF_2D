/**
 * @author Frederic Darboux <frederic.darboux@inra.fr> (2016-2017)
 * @file xyz2asc.c
 * @version 0.00.07
 * @date 2018-02-27
 *
 * @brief FullSWOF_2D XYZ --> ArcGIS ASCII
 * @details
 * Converts a XYZ file for FullSWOF_2D to an ArcGIS ASCII file.
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

/*CHANGES:*/
/*Version 0.00.07: make the program tolerant to empty lines*/
/*Version 0.00.06: better procedure for file declaration*/
/*Version 0.00.05: improve help; improve programming syntax*/
/*Version 0.00.04: improve error messages*/
/*Version 0.00.03: include ctype.h  file to use function tolower()*/
/*Version 0.00.02: use of strcicmp() to make tests case-insensitive*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>

#define FALSE 0
#define TRUE  1
#define ASCHEADER        "#ASCGRID"
#define ASCHEADER_NBCHAR 9
#define ASCHEADER_NBLINE 6

#define LENGTH_OF_LINE 120

typedef float OUTDATA_TYPE;


static char *version = "xyz2asc - 0.00.07 - 2018-02-27";

int strcicmp(char const *a, char const *b);
void headerformat(void);


int main(int argc, char *argv[]){
  char *filename_in, *filename_out;
  FILE *In, *Out;
  char line[LENGTH_OF_LINE], dummy[LENGTH_OF_LINE], varName[LENGTH_OF_LINE];
  char isHeaderLine;
  float x, y, varValue;
  long int i, j, ascHeaderNbLine;
  OUTDATA_TYPE *z;
  long int nbcells, startposition;
  long int ncols, nrows, nbcoord;

  if(argc != 3){
    fprintf(stderr, "\t%s\n", version);
    fprintf(stderr, "\tUsage: %s input.xyz output.asc \n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Note that the input file should contain 3 columns only.\n");
    fprintf(stderr, "If your initial file has more than 3 columns, you may extract the relevant columns using:\n");
    fprintf(stderr, "\tawk '{print $1,$2,$4}' initial.xyz > input-col1-2-4.xyz\n");
    fprintf(stderr, "\n");
    headerformat();
    exit(EXIT_FAILURE);
  }/*if*/

  filename_in  = malloc(strlen(argv[1])+1);
  if(!filename_in){
    fprintf(stderr,"xyz2asc: ERROR: data allocation failure for input file\n");
    exit(EXIT_FAILURE);
  }/*if*/

  filename_out = malloc(strlen(argv[2])+1);
  if(!filename_out){
    fprintf(stderr,"xyz2asc: ERROR: data allocation failure for output file\n");
    exit(EXIT_FAILURE);
  }/*if*/

  strcpy(filename_in,  argv[1]);
  strcpy(filename_out, argv[2]);

  In=fopen(filename_in, "r");
  if(In==NULL){
    fprintf(stderr,"\nxyz2asc: problem opening the input file\n");
    exit(EXIT_FAILURE);
  }/*if*/

  Out=fopen(filename_out, "w");
  if(Out==NULL){
    fprintf(stderr,"\nxyz2asc: problem opening the output file\n");
    exit(EXIT_FAILURE);
  }/*if*/

  /*Read header*/
  ascHeaderNbLine=0;
  isHeaderLine=TRUE;
  while((TRUE == isHeaderLine) && (fgets(line, LENGTH_OF_LINE, In) != NULL)){
    if (line[0]=='#'){
      /*This is a comment line*/
      /*Is this comment line a standard AscGrid header line?*/
      sscanf(line, "%s %s\n", varName, dummy);
      if(strcicmp(varName, ASCHEADER)==0){
	/*save the parameter*/
        if(sscanf(line, "%s %s %f\n", dummy, varName, &varValue) == 3){
        if(strcicmp(varName, "ncols")==0){
	          ncols=(int)varValue;
            ascHeaderNbLine++;
            fprintf(Out,"%s", line+ASCHEADER_NBCHAR);
          }else if(strcicmp(varName, "nrows")==0){
	          nrows=(int)varValue;
            ascHeaderNbLine++;
            fprintf(Out,"%s", line+ASCHEADER_NBCHAR);
          }else if(strcicmp(varName, "xllcorner")==0){
            ascHeaderNbLine++;
            fprintf(Out,"%s", line+ASCHEADER_NBCHAR);
          }else if(strcicmp(varName, "yllcorner")==0){
            ascHeaderNbLine++;
            fprintf(Out,"%s", line+ASCHEADER_NBCHAR);
          }else if(strcicmp(varName, "cellsize")==0){
            ascHeaderNbLine++;
            fprintf(Out,"%s", line+ASCHEADER_NBCHAR);
          }else if(strcicmp(varName, "NODATA_value")==0){
            ascHeaderNbLine++;
            fprintf(Out,"%s", line+ASCHEADER_NBCHAR);
          }else{
	    fprintf(stderr,"\nxyz2asc: Unknown parameter '%s' in the header\n", varName);
	    headerformat();
            exit(EXIT_FAILURE);
          }/*else*/
        }/*if*/
        else{
        /*wrong format of header*/
	  fprintf(stderr,"\nxyz2asc: problem while reading the AscGrid header\n");
	  headerformat();
          exit(EXIT_FAILURE);
        }/*else*/
      }/*if*/
    }/*if*/
    else{/*this is not a header line*/
      isHeaderLine=FALSE;
    }/*else*/
  }/*while*/

  if(ascHeaderNbLine != ASCHEADER_NBLINE){
    fprintf(stderr,"\nxyz2asc: wrong number of lines in the AscGrid header\n");
    headerformat();
    exit(EXIT_FAILURE);
  }/*if*/

  /*Note: First data point is now in line[]*/

  /*Alloc memory*/
  nbcells=ncols*nrows;
  z = malloc((size_t)((nbcells)*sizeof(OUTDATA_TYPE)));
  if(!z){
    fprintf(stderr,"xyz2asc: ERROR: data allocation failure\n");
    exit(EXIT_FAILURE);
  }/*if*/

  /*read the data points*/
  i=0;
  do{
    nbcoord=sscanf(line, "%f %f %f\n", &x, &y, &z[i]);
    if(3==nbcoord){
      /*number of coordinates OK*/
      i++;
    }/*if*/
    else if(nbcoord<1){
      /*probably an empty line ==> do nothing for this line*/
    }/*else if*/
    else{
      /*wrong number of coordinates*/
      fprintf(stderr,"\nxyz2asc: problem while reading a data point (point %li)\n", i);
      exit(EXIT_FAILURE);
    }/*else*/
  }while(fgets(line, LENGTH_OF_LINE, In) != NULL);

  if(i != nbcells){
      fprintf(stderr,"\nxyz2asc: wrong number of data points (%li instead of %li)\n", i, nbcells);
      exit(EXIT_FAILURE);
  }/*if*/

  /*write the data points*/
  startposition=nrows-1;
  for(i=0;i<nrows;i++){
    for(j=0;j<ncols;j++){
      fprintf(Out,"%.6g ", z[startposition+nrows*j-i]); /* %.20g is a trick to avoid trailing zeros*/
    }/*for*/
    fprintf(Out,"\n");
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

void headerformat(void)
{
    fprintf(stderr, "The input file must contain a header looking like:\n");
    fprintf(stderr, "\t#ASCGRID ncols        10000\n");
    fprintf(stderr, "\t#ASCGRID nrows        8500\n");
    fprintf(stderr, "\t#ASCGRID xllcorner    120.2\n");
    fprintf(stderr, "\t#ASCGRID yllcorner    2541.3\n");
    fprintf(stderr, "\t#ASCGRID cellsize     10\n");
    fprintf(stderr, "\t#ASCGRID NODATA_value -9999\n\n");
}
