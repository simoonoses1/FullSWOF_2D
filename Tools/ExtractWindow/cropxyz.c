/**
 * @file cropxyz.c
 * @author Frederic Darboux <frederic.darboux@inra.fr> (2017)
 * @version 0.00.02
 * @date 2017-10-19
 *
 * @brief Crop a FullSWOF_2D XYZ file
 * @details
 * Crop a XYZ file for FullSWOF_2D, i.e. extract a subset of rectangular shape.
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


/*TODO: Process Asc header, if any*/

/*CHANGES:*/
/*2017-10-20 - V0.00.02: Improved programming syntax*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LENGTH_OF_LINE 120

static char *version = "cropxyz - 0.0.2 - 2017-10-20";

int strcicmp(char const *a, char const *b);

int main(int argc, char *argv[]){
  char *filename_in, *filename_out;
  FILE *In, *Out;
  char line[LENGTH_OF_LINE];
  long int lineNumber;
  double x, y, z;
  double xmin, ymin, xmax, ymax;

  if(argc != 7){
    fprintf(stderr, "\t%s\n", version);
    fprintf(stderr, "\tUsage: %s input.xyz xmin ymin xmax ymax output.xyz\n", argv[0]);
    exit(EXIT_FAILURE);
  }/*if*/

  filename_in  = malloc(strlen(argv[1])+1);
  if(!filename_in){
    fprintf(stderr,"cropxyz: ERROR: data allocation failure for input file\n");
    exit(EXIT_FAILURE);
  }/*if*/
  
  filename_out = malloc(strlen(argv[6])+1);
  if(!filename_out){
    fprintf(stderr,"cropxyz: ERROR: data allocation failure for output file\n");
    exit(EXIT_FAILURE);
  }/*if*/
  
  strcpy(filename_in, argv[1]);
  xmin=atof(argv[2]);
  ymin=atof(argv[3]);
  xmax=atof(argv[4]);
  ymax=atof(argv[5]);
  strcpy(filename_out, argv[6]);

  In=fopen(filename_in, "r");
  if(In==NULL){
    fprintf(stderr,"\ncropxyz: problem opening the input file\n");
    exit(EXIT_FAILURE);
  }/*if*/

  Out=fopen(filename_out, "w");
  if(Out==NULL){
    fprintf(stderr,"\ncropxyz: problem opening the output file\n");
    exit(EXIT_FAILURE);
  }/*if*/


  /*Process file*/
  lineNumber=0;
  while(fgets(line, LENGTH_OF_LINE, In) != NULL){
    lineNumber++;
    if (line[0]=='#'){
      /*This is a comment line, copy it*/
	    fprintf(Out,"%s", line);
    }/*if*/
    else{/*this is not a comment line*/
      if(sscanf(line, "%lf %lf %lf\n", &x, &y, &z) == 3){
        /*we got a coordinate point*/
        if((x>=xmin)&&(x<=xmax)&&(y>=ymin)&&(y<=ymax)){
          /*This point is in the requested range*/
          fprintf(Out,"%s", line);
        }/*if*/
      }else{
        fprintf(stderr,"\ncropxyz: line number %li: format problem\n", lineNumber);
        exit(EXIT_FAILURE);
      }/*else*/
    }/*else*/
  }/*while*/

  free(filename_in);
  free(filename_out);
  fclose(In);
  fclose(Out);

  return EXIT_SUCCESS;
}/*main*/
