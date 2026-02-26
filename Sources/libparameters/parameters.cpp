/**
 * @file parameters.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2011-2017)
 * @author Frederic Darboux <frederic.darboux@orleans.inra.fr> (2014)
 * @author: Marco Mancini <marco.mancini@univ-orleans.fr> (2019)
 * @author: Carine Lucas <carine.lucas@univ-orleans.fr> (2020-2023)
 * @version 1.10.00
 * @date 2023-03-17
 *
 * @brief Gets parameters
 * @details
 * Reads the parameters, checks their values.
 *
 * @copyright License Cecill-V2  \n
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

#include "parameters.hpp"


void Parameters::setparameters(const char * FILENAME){

  /**
   * @details
   * Gets all the parameters from the file FILENAME, check and affect them.
   * The values used by FullSWOF\_2D are saved in the file parameters.dat.
   * These values are also printed in the terminal when the code is run.
   * @param[in] FILENAME name of the paramters file.
   * @warning parameters.txt: ERROR: ***.
   * @warning parameters.txt: WARNING: ***.
   * @warning ERROR: the *** file *** does not exist in the directory Inputs.
   * @warning Impossible to open the *** file. Verify if the directory *** exists.
   * @note If a value cannot be affected correctly, the code will exit with failure termination code.\n
   * If parameters.dat cannot be opened, the code will exit with failure termination code.
   */

  // string path_input_directory ("./Inputs/");
  path_input_directory = "./Inputs/";
  string path_output_directory ("./Outputs");
  Parser fileParser(FILENAME);

  /*-----------------------------Nxcell:--------------------------------------------------*/
  Nxcell = atoi(fileParser.getValue("<Nxcell>").c_str());
  if (Nxcell < 1){
    cerr<< " parameters.txt: ERROR: the number of cells must be greater or equal to 1." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------Nycell:--------------------------------------------------*/
  Nycell = atoi(fileParser.getValue("<Nycell>").c_str());
  if (Nycell < 1){
    cerr<< " parameters.txt: ERROR: the number of cells must be greater or equal to 1." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------T:-------------------------------------------------------*/
  T = atof(fileParser.getValue("<T>").c_str());
  if (T < 0.0){
    cerr<< " parameters.txt: ERROR: the final time must be non-negative." << endl;
    exit(EXIT_FAILURE);
  }
  /*-----------------------------nbtimes:-------------------------------------------------*/
  nbtimes = atoi(fileParser.getValue("<nbtimes>").c_str());
  if (nbtimes < 0){
    cerr<< " parameters.txt: ERROR: the number of times saved must be non-negative." << endl;
    exit(EXIT_FAILURE);
  }
  if (1 == nbtimes){
    cerr<< " parameters.txt: WARNING: the initial and the final time are already saved. Using nbtimes=0." << endl;
    nbtimes = 0;
  }

  /*-----------------------------scheme_type:---------------------------------------------*/
  scheme_type = atoi(fileParser.getValue("<scheme_type>").c_str());
  if (scheme_type < 1 || scheme_type > 2){
    cerr << " parameters.txt: ERROR:  the choice "<< scheme_type <<" for the type of scheme does not exist." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------cfl_fix:-------------------------------------------------*/
  cfl_fix = atof(fileParser.getValue("<cflfix>").c_str());
  if (cfl_fix > 1.0){
    cerr<< " parameters.txt: WARNING: you are running the code with CFL = "<< cfl_fix<< endl;
  }
  if (cfl_fix < 0.0){
    cerr<< " parameters.txt: ERROR: the CFL must be non-negative." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------dt_fix:--------------------------------------------------*/
  dt_fix = atof(fileParser.getValue("<dtfix>").c_str());
  if ( 2 == scheme_type){
    if (dt_fix <= 0.0){
      cerr<< " parameters.txt: ERROR: the time step must be non-negative." << endl;
      exit(EXIT_FAILURE);
    }
  }

  /*-----------------------------L:-------------------------------------------------------*/
  L = atof(fileParser.getValue("<L>").c_str());
  if (L <= 0.0){
    cerr<< " parameters.txt: ERROR: the length of the domain must be (strictly) positive." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------l:-------------------------------------------------------*/
  l = atof(fileParser.getValue("<l>").c_str());
  if (l <= 0.0){
    cerr<< " parameters.txt: ERROR: the width of the domain must be (strictly) positive." << endl;
    exit(EXIT_FAILURE);
  }


  dx = L/Nxcell;
  dy = l/Nycell;



  /*-----------------------------Lbound:--------------------------------------------------*/
  //  type of left boundary condition

  string Ltype_bc=fileParser.getValue("<L_bc_init>").c_str();
  if (Ltype_bc.size()<1){
    cerr <<" parameters.txt: ERROR: the type of left boundary condition must be specified." << endl;
    exit(EXIT_FAILURE);
  }
  Lbound_type=atoi(Ltype_bc.c_str());

  switch (Lbound_type){
    case 1 :
    break;
    case 2 ://In case of constant coefficient
    parser_output=fileParser.getValue("<Lbound>").c_str();
    if (parser_output.size()<1){
      cerr <<" parameters.txt: ERROR: the left boundary condition must be specified." << endl;
      exit(EXIT_FAILURE);
    }
    Lbound[1] = atoi(parser_output.c_str());
    if (Lbound[1] < 1 || Lbound[1] > 5){
      cerr <<" parameters.txt: ERROR: the left boundary condition " << Lbound[1]<< " does not exist." << endl;
      exit(EXIT_FAILURE);
    }

    /*-----------------------------left_imp_discharge:--------------------------------------*/
    // left boundary condition : imposed h, q
    parser_Limp_q=fileParser.getValue("<left_imp_discharge>").c_str();
    left_imp_discharge[1] = atof(parser_Limp_q.c_str());

    /*-----------------------------left_imp_h:----------------------------------------------*/
    parser_Limp_h=fileParser.getValue("<left_imp_h>").c_str();
    left_imp_h[1] = atof(parser_Limp_h.c_str());
    break;
    default:
    cerr <<" parameters.txt: ERROR: the type of left boundary condition " << Lbound_type << " does not exist." << endl;
    exit(EXIT_FAILURE);
  }


  /*-----------------------------Rbound:--------------------------------------------------*/
  // type of right boundary condition
  string Rtype_bc=fileParser.getValue("<R_bc_init>").c_str();
  if (Rtype_bc.size()<1){
    cerr <<" parameters.txt: ERROR: the type of right boundary condition must be specified." << endl;
    exit(EXIT_FAILURE);
  }
  Rbound_type=atoi(Rtype_bc.c_str());
  switch (Rbound_type){
    case 1 :
    break;
    case 2 : //In case of constant coefficient
    parser_output=fileParser.getValue("<Rbound>").c_str();
    if (parser_output.size()<1){
      cerr <<" parameters.txt: ERROR: the right boundary condition must be specified." << endl;
      exit(EXIT_FAILURE);
    }
    Rbound[1] = atoi(parser_output.c_str());
    if (Rbound[1] < 1 || Rbound[1] > 5){
      cerr <<" parameters.txt: ERROR: the right boundary condition " << Rbound[1]<< " does not exist." << endl;
      exit(EXIT_FAILURE);
    }

    /*-----------------------------right_imp_discharge:-------------------------------------*/
    parser_Rimp_q=fileParser.getValue("<right_imp_discharge>").c_str();
    right_imp_discharge[1] = atof(parser_Rimp_q.c_str());
    /*-----------------------------right_imp_h:---------------------------------------------*/
    parser_Rimp_h=fileParser.getValue("<right_imp_h>").c_str();
    right_imp_h[1] = atof(parser_Rimp_h.c_str());
    break;

    default :
    cerr <<" parameters.txt: ERROR: the type of right boundary condition " << Rbound_type << " does not exist." << endl;
    exit(EXIT_FAILURE);
  }


  // In case of Left and Right boundary condition are constant coefficients
  if((2==Lbound_type) && (2==Rbound_type)){

    switch (Lbound[1]){    // LEFT BOUNDARY CONDITION

      case 1:  // imposed water height
      if (parser_Limp_q.size()<1){ // if imp_q is empty
        cerr <<" parameters.txt: WARNING: imposed discharge not specified in the left boundary condition."<< endl;
      }
      if (left_imp_h[1] < 0){ // if negative imp_h
        cerr <<" parameters.txt: ERROR: negative imposed water height in the left boundary condition."<< endl;
        exit(EXIT_FAILURE);
      }
      break;

      // case 2: // wall condition
      // nothing to test on these values

      // case 3: // Neumann condition
      // nothing to test on these values

      case 4:  // periodic case
      if (Rbound[1] != 4){  // if the other bound is not periodic
        cerr <<" parameters.txt: ERROR: you must choose a periodic right condition."<< endl;
        exit(EXIT_FAILURE);
      }
      break ;

      case 5:  // imposed discharge
      if (parser_Limp_q.size()<1){ // if imp_q is empty
        cerr <<" parameters.txt: ERROR: imposed discharge not specified in the left boundary condition."<< endl;
        exit(EXIT_FAILURE);
      }
      if (parser_Limp_h.size()<1){ // if imp_h is empty
        cerr <<" parameters.txt: WARNING: imposed water height not specified in the left boundary condition."<< endl;
      }
      if (left_imp_h[1] < 0 ){ // if negative imp_h
        cerr <<" parameters.txt: ERROR: negative imposed water height in the left boundary condition."<< endl;
        exit(EXIT_FAILURE);
      }
      break;

    }


    switch (Rbound[1]){    // RIGHT BOUNDARY CONDITION

      case 1:  // imposed water height
      if (parser_Rimp_q.size()<1){ // if imp_q is empty
        cerr <<" parameters.txt: WARNING: imposed discharge not specified in the right boundary condition."<< endl;
      }
      if (right_imp_h[1] < 0){ // if negative imp_h
        cerr <<" parameters.txt: ERROR: negative imposed water height in the right boundary condition."<< endl;
        exit(EXIT_FAILURE);
      }
      break;

      // case 2: // wall condition
      // nothing to test on these values

      // case 3: // Neumann condition
      // nothing to test on these values

      case 4:  // periodic case
      if (Lbound[1] != 4){  // if the other bound is not periodic
        cerr <<" parameters.txt: ERROR: you must choose a periodic left condition."<< endl;
        exit(EXIT_FAILURE);
      }
      break ;

      case 5:  // imposed discharge
      if (parser_Rimp_q.size()<1){ // if imp_q is empty
        cerr <<" parameters.txt: ERROR: imposed discharge not specified in the right boundary condition."<< endl;
        exit(EXIT_FAILURE);
      }
      if (parser_Rimp_h.size()<1){ // if imp_h is empty
        cerr <<" parameters.txt: WARNING: imposed water height not specified in the right boundary condition."<< endl;
      }
      if (right_imp_h[1] < 0 ){ // if negative imp_h
        cerr <<" parameters.txt: ERROR: negative imposed water height in the right boundary condition."<< endl;
        exit(EXIT_FAILURE);
      }
      break;

    }
  }//end Left and Right boundary condition are constant coefficients



  /*-----------------------------Bbound:--------------------------------------------------*/
  // type of left boundary condition
  string Btype_bc=fileParser.getValue("<B_bc_init>").c_str();
  if (Btype_bc.size()<1){
    cerr <<" parameters.txt: ERROR: the type of bottom boundary condition must be specified." << endl;
    exit(EXIT_FAILURE);
  }
  Bbound_type=atoi(Btype_bc.c_str());


  switch (Bbound_type){ //In case of constant coefficient
    case 1 :
    break;
    case 2 :
    parser_output=fileParser.getValue("<Bbound>").c_str();
    if (parser_output.size()<1){
      cerr <<" parameters.txt: ERROR: the bottom boundary condition must be specified." << endl;
      exit(EXIT_FAILURE);
    }
    Bbound[1] = atoi(parser_output.c_str());
    if (Bbound[1] < 1 || Bbound[1] > 5){
      cerr <<" parameters.txt: ERROR: the bottom boundary condition " << Bbound[1]<< " does not exist." << endl;
      exit(EXIT_FAILURE);
    }

    /*-----------------------------bottom_imp_discharge:------------------------------------*/
    parser_Bimp_q=fileParser.getValue("<bottom_imp_discharge>").c_str();
    bottom_imp_discharge[1] = atof(parser_Bimp_q.c_str());
    /*-----------------------------bottom_imp_h:--------------------------------------------*/
    parser_Bimp_h=fileParser.getValue("<bottom_imp_h>").c_str();
    bottom_imp_h[1] = atof(parser_Bimp_h.c_str());
    break;
    default :
    cerr <<" parameters.txt: ERROR: the type of bottom boundary condition " << Bbound_type << " does not exist." << endl;
    exit(EXIT_FAILURE);
  }



  /*-----------------------------Tbound:--------------------------------------------------*/
  // type of top boundary condition
  string Ttype_bc=fileParser.getValue("<T_bc_init>").c_str();
  if (Ttype_bc.size()<1){
    cerr <<" parameters.txt: ERROR: the type of top boundary condition must be specified." << endl;
    exit(EXIT_FAILURE);
  }
  Tbound_type=atoi(Ttype_bc.c_str());

  switch (Tbound_type){
    case 1 :
    break;
    case 2://In case of constant coefficient
    parser_output=fileParser.getValue("<Tbound>").c_str();
    if (parser_output.size()<1){
      cerr <<" parameters.txt: ERROR: the top boundary condition must be specified." << endl;
      exit(EXIT_FAILURE);
    }
    Tbound[1] = atoi(parser_output.c_str());
    if (Tbound[1] < 1 || Tbound[1] > 5){
      cerr <<" parameters.txt: ERROR: the bottom boundary condition " << Tbound[1]<< " does not exist." << endl;
      exit(EXIT_FAILURE);
    }

    /*-----------------------------top_imp_discharge:---------------------------------------*/
    parser_Timp_q=fileParser.getValue("<top_imp_discharge>").c_str();
    top_imp_discharge[1] = atof(parser_Timp_q.c_str());

    /*-----------------------------top_imp_h:-----------------------------------------------*/
    parser_Timp_h=fileParser.getValue("<top_imp_h>").c_str();
    top_imp_h[1] = atof(parser_Timp_h.c_str());
    break;

    default :
    cerr <<" parameters.txt: ERROR: the type of top boundary condition " << Tbound_type << " does not exist." << endl;
    exit(EXIT_FAILURE);
  }



  // In case of Bottom and Top boundary condition are constant coefficients
  if((2==Tbound_type) && (2==Bbound_type)){
    switch (Bbound[1]){    // BOTTOM BOUNDARY CONDITION

      case 1:  // imposed water height
      if (parser_Bimp_q.size()<1){ // if imp_q is empty
        cerr <<" parameters.txt: WARNING: imposed discharge not specified in the bottom boundary condition."<< endl;
      }
      if (bottom_imp_h[1] < 0){ // if negative imp_h
        cerr <<" parameters.txt: ERROR: negative imposed water height in the bottom boundary condition."<< endl;
        exit(EXIT_FAILURE);
      }
      break;

      // case 2: // wall condition
      // nothing to test on these values

      // case 3: // Neumann condition
      // nothing to test on these values

      case 4:  // periodic case
      if (Tbound[1] != 4){  // if the other bound is not periodic
        cerr <<" parameters.txt: ERROR: you must choose a periodic top condition."<< endl;
        exit(EXIT_FAILURE);
      }
      break ;

      case 5:  // imposed discharge
      if (parser_Bimp_q.size()<1){ // if imp_q is empty
        cerr <<" parameters.txt: ERROR: imposed discharge not specified in the bottom boundary condition."<< endl;
        exit(EXIT_FAILURE);
      }
      if (parser_Bimp_h.size()<1){ // if imp_h is empty
        cerr <<" parameters.txt: WARNING: imposed water height not specified in the bottom boundary condition."<< endl;
      }
      if (bottom_imp_h[1] < 0 ){ // if negative imp_h
        cerr <<" parameters.txt: ERROR: negative imposed water height in the bottom boundary condition."<< endl;
        exit(EXIT_FAILURE);
      }
      break;

    }


    switch (Tbound[1]){    // TOP BOUNDARY CONDITION

      case 1:  // imposed water height
      if (parser_Timp_q.size()<1){ // if imp_q is empty
        cerr <<" parameters.txt: WARNING: imposed discharge not specified in the top boundary condition."<< endl;
      }
      if (top_imp_h[1] < 0){ // if negative imp_h
        cerr <<" parameters.txt: ERROR: negative imposed water height in the top boundary condition."<< endl;
        exit(EXIT_FAILURE);
      }
      break;

      // case 2: // wall condition
      // nothing to test on these values

      // case 3: // Neumann condition
      // nothing to test on these values

      case 4:  // periodic case
      if (Bbound[1] != 4){  // if the other bound is not periodic
        cerr <<" parameters.txt: ERROR: you must choose a periodic bottom condition."<< endl;
        exit(EXIT_FAILURE);
      }
      break;

      case 5:  // imposed discharge
      if (parser_Timp_q.size()<1){ // if imp_q is empty
        cerr <<" parameters.txt: ERROR: imposed discharge not specified in the top boundary condition."<< endl;
        exit(EXIT_FAILURE);
      }
      if (parser_Timp_h.size()<1){ // if imp_h is empty
        cerr <<" parameters.txt: WARNING: imposed water height not specified in the top boundary condition."<< endl;
      }
      if (top_imp_h[1] < 0 ){ // if negative imp_h
        cerr <<" parameters.txt: ERROR: negative imposed water height in the top boundary condition."<< endl;
        exit(EXIT_FAILURE);
      }
      break;

    }
  }//en Bottom and Top boundary condition are constant coefficients


  if(2==Lbound_type) { //In case of constant coefficient
    for (int j=2 ; j<=Nycell ; j++){
      left_imp_discharge[j] = left_imp_discharge[1];

      left_imp_h[j] = left_imp_h[1];

      Lbound[j] = Lbound[1] ;
    }
  }





  if(2==Rbound_type) { //In case of constant coefficient
    for (int j=2 ; j<=Nycell ; j++){
      right_imp_h[j] = right_imp_h[1];

      right_imp_discharge[j] = right_imp_discharge[1];

      Rbound[j] =  Rbound[1];
    }
  }

  if(2==Bbound_type) { //In case of constant coefficient
    for (int i=2 ; i<=Nxcell ; i++){
      bottom_imp_discharge[i] = bottom_imp_discharge[1];

      bottom_imp_h[i] = bottom_imp_h[1];

      Bbound[i] = Bbound[1];
    }
  }

  if(2==Tbound_type) { //In case of constant coefficient

    for (int i=2 ; i<=Nxcell ; i++){
      top_imp_h[i] = top_imp_h[1];

      top_imp_discharge[i] = top_imp_discharge[1];

      Tbound[i] = Tbound[1];
    }
  }


  /*-----------------------------Modification pour BC variable:--------------------------------------------------*/
  string  LBc_NF,BBc_NF,RBc_NF,TBc_NF;

  if(1==Lbound_type) {//initialization from a file
    LBc_NF = fileParser.getValue("<L_bc_NF>");
    left_times_files = verif_file_bc_inhomogeneous(LBc_NF, path_input_directory,'L',
                                                   Lbound, left_imp_discharge,  left_imp_h) ;

  }


  if(1==Rbound_type) {//initialization from a file

    RBc_NF = fileParser.getValue("<R_bc_NF>");
    right_times_files= verif_file_bc_inhomogeneous(RBc_NF, path_input_directory,'R',
                                                   Rbound, right_imp_discharge,  right_imp_h);

  }





  if(1==Bbound_type) {//initialization from a file

    BBc_NF = fileParser.getValue("<B_bc_NF>");
    bottom_times_files = verif_file_bc_inhomogeneous(BBc_NF, path_input_directory,'B',
                                                     Bbound,bottom_imp_discharge,  bottom_imp_h);

  }



  if(1==Tbound_type) {//initialization from a file

    TBc_NF = fileParser.getValue("<T_bc_NF>");
    top_times_files = verif_file_bc_inhomogeneous(TBc_NF, path_input_directory,'T',
                                                  Tbound,top_imp_discharge,  top_imp_h);

  }//


  // In case of Top or Bottom boundary condition is not initialized by constant coefficient
  if((2!=Tbound_type) || (2!=Rbound_type)){

    for (int i=1 ; i<=Nxcell ; i++){
      if (((Tbound[i] == 4) && (Bbound[i] != 4)) || ((Tbound[i] != 4) && (Bbound[i] == 4)) ){
        cerr <<" parameters.txt: ERROR: you must choose a periodic bottom condition."<< endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  // In case of Left or Right boundary condition is not initialized by constant coefficient
  if((2!=Rbound_type) || (2!=Lbound_type)){
    for (int j=1 ; j<=Nycell ; j++){
      if (((Rbound[j] == 4) && (Lbound[j] != 4)) || ((Rbound[j] != 4) && (Lbound[j] == 4)) ){
        cerr <<" parameters.txt: ERROR: you must choose a periodic bottom condition."<< endl;
        exit(EXIT_FAILURE);
      }
    }
  }



  /*-----------------------------fric:----------------------------------------------------*/

  parser_output=fileParser.getValue("<fric>").c_str();
  if (parser_output.size()<1){ // if fric is empty
    cerr <<" parameters.txt: ERROR: the friction law  must be specified." << endl;
    exit(EXIT_FAILURE);
  }

  fric = atoi(parser_output.c_str());

  if (fric < 0 || fric >3){
    cerr <<" parameters.txt: ERROR: the friction law " << fric << " does not exist."<< endl;
    exit(EXIT_FAILURE);
  }

  if (fric > 0){//if the case differs from no friction
    if (fileParser.getValue("<fric_init>").size()<1){  // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of friction term (<fric_init>)."<<endl;
      exit(EXIT_FAILURE);
    }

    fric_init = atoi(fileParser.getValue("<fric_init>").c_str());

    if (fric_init <1 || fric_init >2){
      cerr <<" parameters.txt: ERROR: initialization number " << fric_init << " for friction term (<fric_init>) does not exist."<< endl;
      exit(EXIT_FAILURE);
    }

    if (1==fric_init){//initialization from a file
      fric_NF = fileParser.getValue("<fric_NF>");
      if (fric_NF.size()<1){  // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of friction."<<endl;
        exit(EXIT_FAILURE);
      }
      fric_namefile = path_input_directory + fric_NF;
      if (-1 == access (fric_namefile.c_str(), R_OK) ) {
        cerr << " parameters.txt: ERROR: the "<< fric_NF<< " file does not exist in the directory Inputs." << endl ;
        exit(EXIT_FAILURE);

      }
    }
    else{

      /*-----------------------------fric_coef:------------------------------------------------*/
      friccoef = atof(fileParser.getValue("<friccoef>").c_str());

      if (friccoef < 0){//if negative
        cerr <<" parameters.txt: ERROR: negative friction coefficient."<< endl;
        exit(EXIT_FAILURE);
      }
      if (0. >= friccoef){
        cerr << " parameters.txt: WARNING: friction law with null friction coefficient " <<endl;
        cerr << "*********************************************************" << endl;
      }
    }
  }

  /*-----------------------------flux:----------------------------------------------------*/
  parser_output=fileParser.getValue("<flux>").c_str();
  if (parser_output.size()<1){ // if flux is empty
    cerr << " parameters.txt: WARNING: the numerical flux is empty, using HLL instead"<<endl;
    cerr << "*********************************************************" << endl;
    flux = 2;
  }
  else{
    flux = atoi(parser_output.c_str());
    if (flux < 1 || flux >5){
      cerr <<" parameters.txt: ERROR: the numerical flux " << flux << " does not exist."<< endl;
      exit(EXIT_FAILURE);
    }
  }

  /*-----------------------------order:---------------------------------------------------*/

  parser_output=fileParser.getValue("<order>").c_str();
  if (parser_output.size()<1){ // if order is empty
    cerr << " parameters.txt: WARNING: the order of the scheme is empty, using order 2 instead"<<endl;
    cerr << "*********************************************************" << endl;
    order = 2;
  }
  else{
    order = atoi(parser_output.c_str());

    if (order< 1 || order >2){
      cerr << "parameters.txt: ERROR: the order " << order  <<" does not exist for this scheme."<<endl;
      exit(EXIT_FAILURE);
    }
  }

  /*-----------------------------rec:-----------------------------------------------------*/
  if (2 == order){// in case of order 2 only
    parser_output=fileParser.getValue("<rec>").c_str();
    if (parser_output.size()<1){ // if rec is empty
      cerr << " parameters.txt: WARNING: the reconstruction is empty, using MUSCL instead"<<endl;
      cerr << "*********************************************************" << endl;
      rec = 1; //using MUSCL by default
    }
    else{
      rec = atoi(parser_output.c_str());
      if (rec< 1 || rec >3){
        cerr << " parameters.txt: ERROR: the reconstruction " << rec  <<" does not exist."<<endl;
        exit(EXIT_FAILURE);
      }

      if (rec >1){ //in case of ENO or ENOmod
        /*-----------------------------amortENO:------------------------------------------------*/
        parser_output2=fileParser.getValue("<amortENO>").c_str();
        if (parser_output2.size()<1){ // if amortENO is empty
          cerr << " parameters.txt: ERROR: the AmortENO is empty"<<endl;
          cerr << "*********************************************************" << endl;
          exit(EXIT_FAILURE);
        }
        else{
          amortENO = atof(parser_output2.c_str());
          if (amortENO< 0 || amortENO >1){
            cerr << "parameters.txt: ERROR: the AmortENO " << amortENO  <<" is not between 0 and 1."<<endl;
            exit(EXIT_FAILURE);
          }
        }
        if (3 == rec){//in case of  ENOmod
          /*-----------------------------modifENO:------------------------------------------------*/
          parser_output2=fileParser.getValue("<modifENO>").c_str();
          if (parser_output2.size()<1){ // if modifENO is empty
            cerr << " parameters.txt: ERROR: the modifENO is empty"<<endl;
            cerr << "*********************************************************" << endl;
            exit(EXIT_FAILURE);
          }
          else{
            modifENO = atof(parser_output2.c_str());
            if (modifENO< 0 || modifENO >1){
              cerr << " parameters.txt: ERROR: the modifENO " << modifENO <<" is not between 0 and 1."<<endl;
              exit(EXIT_FAILURE);
            }
          }
        }//endif rec==3

      }//endif in case of ENO or ENOmod
    }//end rec treatment

    /*-----------------------------lim:-----------------------------------------------------*/
    parser_output=fileParser.getValue("<lim>").c_str();
    if (parser_output.size()<1){ // if lim is empty
      cerr << " parameters.txt: WARNING: the limiter is empty, using Minmod instead"<<endl;
      cerr << "*********************************************************" << endl;
      lim = 1; //using Minmod by default
    }
    else{
      lim = atoi(parser_output.c_str());
      if (lim< 1 || lim >3){
        cerr << "parameters.txt: ERROR: the limiter " << lim  <<" does not exist."<<endl;
        exit(EXIT_FAILURE);
      }
    }
  }//endif order 2

  /*-----------------------------inf:-----------------------------------------------------*/
  parser_output=fileParser.getValue("<inf>").c_str();
  if (parser_output.size()<1){ // if inf is empty
    cerr <<" parameters.txt: ERROR: the infiltration model  must be specified." << endl;
    exit(EXIT_FAILURE);
  }

  inf = atoi(parser_output.c_str());

  if (inf < 0 || inf >1){
    cerr <<" parameters.txt: ERROR: the infiltration model " << inf << " does not exist."<< endl;
    exit(EXIT_FAILURE);
  }
  if (inf > 0){//if the case differs from no infiltration
    /*-----------------------------Kc:------------------------------------------------------*/
    if (fileParser.getValue("<Kc_init>").size()<1){  // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of hydraulic conductivity of the crust (<Kc_init>)."<<endl;
      exit(EXIT_FAILURE);
    }

    Kc_init = atoi(fileParser.getValue("<Kc_init>").c_str());



    switch (Kc_init){
      case 1://initialization from a file
      Kc_NF = fileParser.getValue("<Kc_NF>");
      if (Kc_NF.size()<1){  // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of Kc."<<endl;
        exit(EXIT_FAILURE);
      }
      Kc_namefile = path_input_directory + Kc_NF;
      if (-1 == access (Kc_namefile.c_str(), R_OK) ) {
        cerr << " parameters.txt: ERROR: the "<< Kc_NF<< " file does not exist in the directory Inputs." << endl ;
        exit(EXIT_FAILURE);
      }
      break;
      case 2:
      Kc_coef = atof(fileParser.getValue("<Kccoef>").c_str());

      if (Kc_coef < 0){//if hydraulic conductivity of the crust negative
        cerr <<" parameters.txt: ERROR: negative hydraulic conductivity of the crust."<< endl;
        exit(EXIT_FAILURE);
      }

      if (0. >= Kc_coef){
        cerr << " parameters.txt: WARNING: infiltration model with null hydraulic conductivity of the crust  " <<endl;
        cerr << "*********************************************************" << endl;
      }
      break;
      default :
      cerr <<" parameters.txt: ERROR: initialization number " << Kc_init << " for hydraulic conductivity of the crust (<Kc_init>) does not exist."<< endl;
      exit(EXIT_FAILURE);
    }


    /*-----------------------------Ks:------------------------------------------------------*/
    if (fileParser.getValue("<Ks_init>").size()<1){  // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of hydraulic conductivity of the soil (<Ks_init>)."<<endl;
      exit(EXIT_FAILURE);
    }

    Ks_init = atoi(fileParser.getValue("<Ks_init>").c_str());


    switch (Ks_init){
      case 1://initialization from a file
      Ks_NF = fileParser.getValue("<Ks_NF>");
      if (Ks_NF.size()<1){  // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of Ks."<<endl;
        exit(EXIT_FAILURE);
      }
      Ks_namefile = path_input_directory + Ks_NF;
      if (-1 ==  access (Ks_namefile.c_str(), R_OK)) {
        cerr << " parameters.txt: ERROR: the "<< Ks_NF<< " file does not exist in the directory Inputs." << endl ;
        exit(EXIT_FAILURE);
      }
      break;
      case 2:
      Ks_coef = atof(fileParser.getValue("<Kscoef>").c_str());
      if (Ks_coef < 0){//if hydraulic conductivity of the soil negative
        cerr <<" parameters.txt: ERROR: negative hydraulic conductivity of the soil."<< endl;
        exit(EXIT_FAILURE);
      }
      break;
      default:
      cerr <<" parameters.txt: ERROR: initialization number " << Ks_init << " for hydraulic conductivity of the soil (<Ks_init>) does not exist."<< endl;
      exit(EXIT_FAILURE);
    }

    /*-----------------------------dtheta:--------------------------------------------------*/
    if (fileParser.getValue("<dtheta_init>").size()<1){  // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of initial water deficit (<dtheta_init>)."<<endl;
      exit(EXIT_FAILURE);
    }

    dtheta_init = atoi(fileParser.getValue("<dtheta_init>").c_str());

    switch (dtheta_init){
      case 1: //initialization from a file
      dtheta_NF = fileParser.getValue("<dtheta_NF>");
      if (dtheta_NF.size()<1){  // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of dtheta."<<endl;
        exit(EXIT_FAILURE);
      }
      dtheta_namefile = path_input_directory + dtheta_NF;
      if (-1 ==  access (dtheta_namefile.c_str(), R_OK)) {
        cerr << " parameters.txt: ERROR: the "<< dtheta_NF<< " file does not exist in the directory Inputs." << endl ;
        exit(EXIT_FAILURE);
      }
      break;
      case 2:
      dtheta_coef = atof(fileParser.getValue("<dthetacoef>").c_str());
      if (dtheta_coef < 0){//if initial water deficit is negative
        cerr <<" parameters.txt: ERROR: negative initial water deficit."<< endl;
        exit(EXIT_FAILURE);
      }
      if (dtheta_coef > 1){//if initial water deficit is above 1
        cerr <<" parameters.txt: WARNING: the value of dtheta seems very large."<< endl;
      }
      break;
      default:
      cerr <<" parameters.txt: ERROR: initialization number " << dtheta_init << " for initial water deficit (<dtheta_init>) does not exist."<< endl;
      exit(EXIT_FAILURE);
    }

    /*-----------------------------Psi:-----------------------------------------------------*/
    parser_output = fileParser.getValue("<Psi_init>");
    if (parser_output.size()<1){  // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of load pressure (<Psi_init>)."<<endl;
      exit(EXIT_FAILURE);
    }

    Psi_init = atoi(parser_output.c_str());
    switch (Psi_init){
      case 1:    /* initialization from a file */
      Psi_NF = fileParser.getValue("<Psi_NF>");
      if (Psi_NF.size()<1){  // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of Psi."<<endl;
        exit(EXIT_FAILURE);
      }
      Psi_namefile = path_input_directory + Psi_NF;
      if ( -1 == access (Psi_namefile.c_str(), R_OK)) {
        cerr << " parameters.txt: ERROR: the "<< Psi_NF<< " file does not exist in the directory Inputs." << endl ;
        exit(EXIT_FAILURE);
      }
      break;

      case 2 : //constant coefficient
      Psi_coef = atof(fileParser.getValue("<Psicoef>").c_str());
      if (Psi_coef < 0){//if load pressure negative
        cerr <<" parameters.txt: ERROR: negative load pressure."<< endl;
        exit(EXIT_FAILURE);
      }
      break;

      default :
      cerr <<" parameters.txt: ERROR: initialization number " << Psi_init << " for load pressure (<Psi_init>) does not exist."<< endl;
      exit(EXIT_FAILURE);
    }

    /*-----------------------------zcrust:--------------------------------------------------*/
    if (fileParser.getValue("<zcrust_init>").size()<1){  // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of thickness of the crust (<zcrust_init>)."<<endl;
      exit(EXIT_FAILURE);
    }

    zcrust_init = atoi(fileParser.getValue("<zcrust_init>").c_str());

    switch (zcrust_init){
      case 1://initialization from a file
      zcrust_NF = fileParser.getValue("<zcrust_NF>");
      if (zcrust_NF.size()<1){  // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of zcrust."<<endl;
        exit(EXIT_FAILURE);
      }
      zcrust_namefile = path_input_directory + zcrust_NF;
      if (-1 ==  access (zcrust_namefile.c_str(), R_OK)) {
        cerr << " parameters.txt: ERROR: the "<< zcrust_NF<< " file does not exist in the directory Inputs." << endl ;
        exit(EXIT_FAILURE);
      }
      break;
      case 2:
      zcrust_coef = atof(fileParser.getValue("<zcrustcoef>").c_str());
      if (zcrust_coef < 0){//if thickness of the crust negative
        cerr <<" parameters.txt: ERROR: negative thickness of the crust."<< endl;
        exit(EXIT_FAILURE);
      }
      break;
      default:
      cerr <<" parameters.txt: ERROR: initialization number " << zcrust_init << " for thickness of the crust (<zcrust_init>) does not exist."<< endl;
      exit(EXIT_FAILURE);
    }

    /*-----------------------------imax:----------------------------------------------------*/
    if (fileParser.getValue("<imax_init>").size()<1){  // if the choice is empty
      cerr << " parameters.txt: ERROR: it is necessary to choose the initialization of Maximum infiltration rate (<imax_init>)."<<endl;
      exit(EXIT_FAILURE);
    }

    imax_init = atoi(fileParser.getValue("<imax_init>").c_str());

    switch (imax_init){
    case 1://initialization from a file
      imax_NF = fileParser.getValue("<imax_NF>");
      if (imax_NF.size()<1){  // if file's name is empty
        cerr << " parameters.txt: ERROR: it is necessary to have a file for this choice of initialization of imax."<<endl;
        exit(EXIT_FAILURE);
      }
      imax_namefile = path_input_directory + imax_NF;
      if (-1 == access (imax_namefile.c_str(), R_OK)) {
        cerr << " parameters.txt: ERROR: the "<< imax_NF<< " file does not exist in the directory Inputs." << endl ;
        exit(EXIT_FAILURE);
      }
      break;
      case 2:
      imax_coef = atof(fileParser.getValue("<imaxcoef>").c_str());
      if (imax_coef < 0){//if Maximum infiltration rate negative
        cerr <<" parameters.txt: ERROR: negative Maximum infiltration rate."<< endl;
        exit(EXIT_FAILURE);
      }
      break;
      default:
      cerr <<" parameters.txt: ERROR: initialization number " << imax_init << " for Maximum infiltration rate (<imax_init>) does not exist."<< endl;
      exit(EXIT_FAILURE);
    }
  }//endif the case differs from no infiltration

  /*-----------------------------topo:----------------------------------------------------*/
  topo = atoi(fileParser.getValue("<topo>").c_str());
  if (topo < 1 || topo > 3){
    cerr << " parameters.txt: ERROR: initialization number "<< topo <<" for the topography does not exist." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------topo_NF:-------------------------------------------------*/
  topo_NF = fileParser.getValue("<topo_NF>");
  topography_namefile = path_input_directory + topo_NF;
  if (1 == topo && -1 == access (topography_namefile.c_str(), R_OK) ) {
    cerr << " parameters.txt: ERROR: the topography file "<< topo_NF<< " does not exist in the directory Inputs." << endl ;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------huv_init:------------------------------------------------*/
  huv_init = atoi(fileParser.getValue("<huv_init>").c_str());
  if (huv_init < 1 || huv_init > 5){
    cerr << " parameters.txt: ERROR: initialization number "<< huv_init <<" for h, u does not exist." << endl;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------huv_NF:--------------------------------------------------*/
  huv_NF = fileParser.getValue("<huv_NF>");
  huv_namefile = path_input_directory + huv_NF;
  if (1 == huv_init && -1 == access (huv_namefile.c_str(), R_OK)) {
    cerr << " parameters.txt: ERROR: the huv file "<< huv_NF<< " does not exist in the directory Inputs." << endl ;
    exit(EXIT_FAILURE);
  }

  /*-----------------------------rain:----------------------------------------------------*/
  rain = atoi(fileParser.getValue("<rain>").c_str());

  switch (rain) {
    case 0 :
    break;
    case 1 :
    /*-----------------------------rain_NF:-------------------------------------------------*/
    rain_NF = fileParser.getValue("<rain_NF>");
    rain_namefile = path_input_directory + rain_NF;
    if (-1 == access (rain_namefile.c_str(), R_OK)) {
      cerr << " parameters.txt: ERROR: the rain file "<< rain_NF  << " does not exist in the directory Inputs" << endl ;
      exit(EXIT_FAILURE);
    }
    break;
    case 2 :
    break;
    default :
    cerr << " parameters.txt: ERROR: the choice "<< rain <<" for rain does not exist." << endl;
    exit(EXIT_FAILURE);
  }



  /*-----------------------------suffix_o:------------------------------------------------*/
  suffix_outputs = fileParser.getValue("<suffix_o>");

  /*-----------------------------output_f:--------------------------------------------------------*/
  if (0==nbtimes){
    output_format=0;
  }
  else{
    output_format = atoi(fileParser.getValue("<output_f>").c_str());
    if (output_format < 1 || output_format > 2){
      cerr << " parameters.txt: ERROR: the choice "<< output_format <<" for the format of the output file does not exist." << endl;
      exit(EXIT_FAILURE);
    }
  }



  /*-----------------------------Choice of the specific points--------------------------------------------------------*/
  /*defaults*/
  Choice_points = 0;
  Choice_dt_specific_points = 1;

  parser_specific_points = fileParser.getValue("<Choice_points>").c_str();
  if (parser_specific_points.size()<1){ // if the choice of specific_points is empty
    cerr <<" parameters.txt: ERROR: the  choice of specific points must be specified." << endl;
    exit(EXIT_FAILURE);
  }
  Choice_points = atoi(parser_specific_points.c_str());



  switch (Choice_points){
    case 0:
    break;
    /*A single point is given in parameters.txt file*/
    case 1 :
    if (fileParser.getValue("<x_coord>").size()<1){  //
      cerr << " parameters.txt: ERROR: it is necessary to put a value for x coordinate."<<endl;
      exit(EXIT_FAILURE);
    }
    x_coord = atof(fileParser.getValue("<x_coord>").c_str());
    if (fileParser.getValue("<y_coord>").size()<1){  //
      cerr << " parameters.txt: ERROR: it is necessary to put a value for y coordinate."<<endl;
      exit(EXIT_FAILURE);
    }
    y_coord = atof(fileParser.getValue("<y_coord>").c_str());
    if(!is_coord_in_file_valid(x_coord,y_coord,-1,"parameters.txt")){
      exit(EXIT_FAILURE);
    }
    break;

    /* A list of specific points is given in a file*/
    case 2:
    list_point_NF = fileParser.getValue("<list_point_NF>");
    if (list_point_NF.size()<1){  // if file's name is empty
      cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of list of points."<<endl;
      exit(EXIT_FAILURE);
    }
    list_point_namefile = path_input_directory + list_point_NF;
    if (-1 == access (list_point_namefile.c_str(), R_OK) ) {
      cerr << " parameters.txt: ERROR: the "<< list_point_NF<< " file does not exist in the directory Inputs." << endl ;
      exit(EXIT_FAILURE);
    }
    break;

    default :
    cerr <<" parameters.txt: ERROR: the choice " << Choice_points << " for the specific points which must be saved does not exist."<< endl;
    exit(EXIT_FAILURE);
  }

  if (Choice_points > 0) { // we don't choose a time saved in the case 0 (0 = no save)

    /*Choice of time step for specific points. */
    Choice_dt_specific_points=atoi(fileParser.getValue("<Choice_dt_specific_points>").c_str());

    switch (Choice_dt_specific_points){
      case 1: // User chooses to save each time
      dt_specific_points = 0.0;
      break;
      case 2: // User imposes a given time step
      if (fileParser.getValue("<dt_specific_points>").size()<1){  //
        cerr << " parameters.txt: ERROR: it is necessary to specify a value for time step to save the specific points."<<endl;
        exit(EXIT_FAILURE);
      }
      dt_specific_points = atof(fileParser.getValue("<dt_specific_points>").c_str());
      break;
      default:
      cerr <<" parameters.txt: ERROR: the choice " << Choice_dt_specific_points << " in the choice of times saved for the specific points does not exist."<< endl;
      exit(EXIT_FAILURE);
    }
  }


  output_directory = path_output_directory + suffix_outputs + "/";

  namefile = output_directory + "parameters.dat";

  ofstream param(namefile.c_str(),ios::out);
  if (!param){
    cerr << " parameters.txt: ERROR: Impossible to open the " << namefile.c_str() <<" file\n";
    cerr << " parameters.txt: ERROR: Verify if the directory " << output_directory <<" exists\n";
    exit(EXIT_FAILURE);
  }

  param << "#####################################################################" <<endl;
  param << "# Parameters used in FullSWOF_2D software."<<endl;
  param << "#####################################################################" <<endl;
  param << endl;
  param << "Number of grid cells with respect to x <Nxcell>:: "<< Nxcell << endl;
  param << "Number of grid cells with respect to y <Nycell>:: "<< Nycell << endl;
  param << "Length of the domain with respect to x (in meters) <L>:: "<< L << endl;
  param << "Length of the domain with respect to y (in meters) <l>:: "<< l << endl;
  param << endl;
  param << "Final time <T>:: " << T << endl;
  param << "Number of times saved <nbtimes>:: " << nbtimes << endl;
  param << "Choice of type of scheme (1=fixed cfl 2=fixed dt) <scheme_type>:: " << scheme_type << endl;
  if (2 == scheme_type){ // if fixed dt
    param << "Timestep (in seconds) <dtfix>:: "<< dt_fix << endl;
  }
  else{
    param << "Timestep (in seconds) <dtfix>:: " << endl;
  }
  param << "Value of the cfl <cflfix>:: "<< cfl_fix << endl;
  param << endl;


  param << "Left boundary condition choice (1=file 2=const_coef) <L_bc_init>:: "<< Lbound_type<< endl;

  if(1==Lbound_type) { //In case of constant coefficient
    param << "Name of the left boundary condition file <L_bc_NF>:: "<< LBc_NF << endl<< endl;
    param << "Left boundary condition (x = 0) (1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q) <Lbound>:: " << endl;
    param << "Imposed discharge in left bc <left_imp_discharge>:: " << endl; // empty as nothing is used
    param << "Imposed water height in left bc (if flow supercritical) <left_imp_h>:: " << endl; // empty as nothing is used
  }else{

    param << "Name of the left boundary condition file <L_bc_NF>:: " << endl<< endl;
    param << "Left boundary condition (x = 0) (1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q) <Lbound>:: " << Lbound[1] << endl;
    if (2 == Lbound[1] || 3 == Lbound[1]||4 == Lbound[1]){ // if wall, Neumann, or periodic
      param << "Imposed discharge in left bc <left_imp_discharge>:: " << endl; // empty as nothing is used
      param << "Imposed water height in left bc (if flow supercritical) <left_imp_h>:: " << endl; // empty as nothing is used
    }
    else{
      param << "Imposed discharge in left bc <left_imp_discharge>:: " << left_imp_discharge[1] << endl; // writes 0 if the value was empty
      param << "Imposed height in left bc (if flow supercritical) <left_imp_h>:: " << left_imp_h[1] << endl;  // writes 0 if the value was empty
    }
  }

  param << endl;


  param << "Right boundary condition choice (1=file 2=const_coef) <R_bc_init>:: "<< Rbound_type << endl;
  if(1==Rbound_type) { //In case of constant coefficient
    param << "Name of the right boundary condition file <R_bc_NF>:: "<< RBc_NF << endl << endl;
    param << "Right boundary condition (x = xmax)(1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q) <Rbound>:: " << endl;
    param << "Imposed discharge in right bc <right_imp_discharge>:: " << endl; // empty as nothing is used
    param << "Imposed water height in right bc (if flow supercritical) <right_imp_h>:: " << endl; // empty as nothing is used
  }else{
    param << "Name of the right boundary condition file <R_bc_NF>:: "<< endl << endl;
    param << "Right boundary condition (x = xmax)(1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q) <Rbound>:: " << Rbound[1] << endl;
    if (2 == Rbound[1] ||3 == Rbound[1] || 4 == Rbound[1]){ // if wall, Neumann, or periodic
      param << "Imposed discharge in right bc (if flow supercritical) <right_imp_discharge>:: " << endl; // empty as nothing is used
      param << "Imposed water height in right bc <right_imp_h>:: " << endl; // empty as nothing is used
    }
    else{
      param << "Imposed discharge in right bc <right_imp_discharge>:: " << right_imp_discharge[1] << endl; // writes 0 if the value was empty
      param << "Imposed height in right bc (if flow supercritical) <right_imp_h>:: " << right_imp_h[1] << endl; // writes 0 if the value was empty
    }
  }
  param << endl;



  param << "Bottom boundary condition choice (1=file 2=const_coef) <B_bc_init>:: "<< Bbound_type << endl;
  if(1==Bbound_type) { //In case of constant coefficient
    param << "Name of the bottom boundary condition file <B_bc_NF>:: "<< BBc_NF << endl << endl;
    param << "Bottom boundary condition (y = 0) (1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q) <Bbound>:: " << endl;
    param << "Imposed discharge in bottom bc <bottom_imp_discharge>:: " << endl; // empty as nothing is used
    param << "Imposed water height in bottom bc (if flow supercritical) <bottom_imp_h>:: " << endl; // empty as nothing is used
  }else{
    param << "Name of the bottom boundary condition file <B_bc_NF>:: "<< endl << endl;
    param << "Bottom boundary condition (y = 0) (1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q) <Bbound>:: " << Bbound[1] << endl;
    if (2 == Bbound[1] || 3 == Bbound[1] || 4 == Bbound[1]){ // if wall, Neumann, or periodic
      param << "Imposed discharge in bottom bc <bottom_imp_discharge>:: " << endl; // empty as nothing is used
      param << "Imposed water height in bottom bc (if flow supercritical) <bottom_imp_h>:: " << endl; // empty as nothing is used
    }
    else{
      param << "Imposed discharge in bottom bc <bottom_imp_discharge>:: " << bottom_imp_discharge[1] << endl;// empty as nothing is used
      param << "Imposed height in bottom bc (if flow supercritical) <bottom_imp_h>:: " << bottom_imp_h[1] << endl;// empty as nothing is used
    }

  }


  param << endl;

  param << "Top boundary condition choice (1=file 2=const_coef) <T_bc_init>:: " << Tbound_type << endl;
  if(1==Tbound_type) { //In case of constant coefficient
    param << "Name of the Top boundary condition file <T_bc_NF>:: " << TBc_NF << endl << endl;
    param << "Top boundary condition (y = ymax) (1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q) <Tbound>:: " << endl;
    param << "Imposed discharge in top bc <top_imp_discharge>:: " << endl; // empty as nothing is used
    param << "Imposed water height in top bc (if flow supercritical) <top_imp_h>:: " << endl; // empty as nothing is used
  }else{
    param << "Name of the Top boundary condition file <T_bc_NF>:: "<< endl << endl;
    param << "Top boundary condition (y = ymax) (1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q) <Tbound>:: " << Tbound[1] << endl;
    if (2 == Tbound[1] || 3 == Tbound[1] || 4 == Tbound[1]){ // if wall, Neumann, or periodic
      param << "Imposed discharge in top bc <top_imp_discharge>:: " << endl; // empty as nothing is used
      param << "Imposed water height in top bc (if flow supercritical) <top_imp_h>:: " << endl; // empty as nothing is used
    }
    else{
      param << "Imposed discharge in top bc <top_imp_discharge>:: " << top_imp_discharge[1] << endl;// empty as nothing is used
      param << "Imposed height in top bc (if flow supercritical) <top_imp_h>:: " << top_imp_h[1] << endl;// empty as nothing is used
    }

  }

  param << endl;
  param << "Friction law (0=NoFriction 1=Manning 2=Darcy-Weisbach 3=laminar) <fric>:: "<< fric << endl;
  if (fric > 0){
    param << "Initialization of Friction (1=file 2=const_coef) <fric_init>:: "<< fric_init << endl;
  }
  else{
    param << "Initialization of Friction (1=file 2=const_coef) <fric_init>:: " << endl;
  }
  if (fric_init==1){
    param <<"Name of the friction file <fric_NF>:: "<< fric_NF<< endl;
  }
  else{
    param <<"Name of the friction file <fric_NF>:: " << endl;
  }
  if (fric_init==2){
    param << "Friction coefficient <friccoef>:: "<< friccoef << endl;
  }
  else{
    param << "Friction coefficient <friccoef>:: "<< endl;
  }
  param << endl;

  param << "Numerical flux (1=Rusanov 2=HLL 3=HLL2 4=HLLC 5=HLLC2) <flux>:: " << flux << endl;
  param << endl;
  param << "Order of the scheme <order>:: " << order << endl;
  param << endl;
  if (2==order){// in case of order 2
    param << "Reconstruction (1=MUSCL 2=ENO 3=ENOmod) <rec>:: " << rec << endl;
    if (1 == rec){
      param << "AmortENO <amortENO>:: "<< endl;
      param << "ModifENO <modifENO>:: "<< endl;
    }
    else{
      param << "AmortENO <amortENO>:: " << amortENO << endl;
      if (3 == rec){
        param << "ModifENO <modifENO>:: " << modifENO << endl;
      }
      else{
        param << "ModifENO <modifENO>:: " << endl;
      }
    }
    param << "Limiter (1=Minmod 2=VanAlbada 3=VanLeer) <lim>:: " << lim << endl;
  }
  else{
    param << "Reconstruction (1=MUSCL 2=ENO 3=ENOmod) <rec>:: "<< endl;
    param << "AmortENO <amortENO>:: "<< endl;
    param << "ModifENO <modifENO>:: "<< endl;
    param << "Limiter (1=Minmod 2=VanAlbada 3=VanLeer) <lim>:: "<< endl;
  }

  param << endl;
  param << "Topography (1=file 2=flat 3=Thacker) <topo>:: " << topo << endl;
  param << "Name of the topography file <topo_NF>:: "<< topo_NF << endl;
  param << endl;
  param << "Initialization of h, u and v (1=file 2=h,u&v=0 3=Thacker 4=Radial_Dam_dry 5=Radial_Dam_wet) <huv_init>:: " << huv_init << endl;
  param << "Name of the huv initialization file <huv_NF>:: " << huv_NF << endl;
  param << endl;
  param << "Rain (0=no rain 1=file 2=function) <rain>:: " << rain << endl;
  if (1==rain){// in case of rain read in a file
    param << "Name of the rain file <rain_NF>:: "<< rain_NF <<endl;
  }
  else{
    param << "Name of the rain file <rain_NF>:: "<< endl;
  }
  param << endl;

  param << "Infiltration model (0=No infiltration 1=Green-Ampt) <inf>:: " << inf << endl;
  if (0 == inf){
    param << "zcrust, thickness of the crust (1=file 2=const_coef) <zcrust_init>:: "<< endl;
    param <<"zcrust coefficient <zcrustcoef>:: "<< endl;
    param <<"Name of the zcrust file <zcrust_NF>:: "<< endl;
    param << endl;

    param << "Ks, hydraulic conductivity (saturation) of the soil (1=file 2=const_coef) <Ks_init_init>:: "<< endl;
    param <<"Ks coefficient <Kscoef>:: "<< endl;
    param <<"Name of the Ks file <Ks_NF>:: "<< endl;
    param << endl;

    param << "Kc, hydraulic conductivity (saturation) of the crust (1=file 2=const_coef) <Kc_init>:: "<< endl;
    param <<"Kc coefficient <Kccoef>:: "<< endl;
    param <<"Name of the Kc file <Kc_NF>:: "<< endl;
    param << endl;
    
    param << "dtheta, initial water deficit (1=file 2=const_coef) <dtheta_init>:: "<< endl;
    param <<"dtheta coefficient <dthetacoef>:: "<< endl;
    param <<"Name of the dtheta file <dtheta_NF>:: "<< endl;
    param << endl;

    param << "Psi, load pressure (1=file 2=const_coef) (1=file 2=const_coef) <Psi_init>:: "<< endl;
    param <<"Psi coefficient <Psicoef>:: "<< endl;
    param <<"Name of the dtheta file <dtheta_NF>:: "<< endl;
    param << endl;

    param << "imax, maximum infiltration rate (1=file 2=const_coef) <imax_init>:: "<< endl;
    param <<"imax coefficient <imaxcoef>:: "<< endl;
    param <<"Name of the imax file <imax_NF>:: "<< endl;
  }
  else{

    param << "zcrust, thickness of the crust (1=file 2=const_coef) <zcrust_init>:: " << zcrust_init << endl;
    if (zcrust_init==1){
      param <<"zcrust coefficient <zcrustcoef>:: "<< endl;
      param <<"Name of the zcrust file <zcrust_NF>:: "<< zcrust_NF<< endl;
    }
    else{
      param <<"zcrust coefficient <zcrustcoef>:: "<< zcrust_coef <<endl;
      param <<"Name of the zcrust file <zcrust_NF>:: "<< endl;
    }
    param << endl;


    param << "Ks, hydraulic conductivity (saturation) of the soil (1=file 2=const_coef) <Ks_init>:: " << Ks_init << endl;
    if (Ks_init==1){
      param <<"Ks coefficient <Kscoef>:: "<< endl;
      param <<"Name of the Ks file <Ks_NF>:: "<< Ks_NF <<endl;
    }
    else{
      param <<"Ks coefficient <Kscoef>:: "<< Ks_coef <<endl;
      param <<"Name of the Ks file <Ks_NF>:: "<< endl;
    }
    param << endl;

    param << "Kc, hydraulic conductivity (saturation) of the crust (1=file 2=const_coef) <Kc_init>:: " << Kc_init << endl;
    if (Kc_init==1){
      param <<"Kc coefficient <Kccoef>:: "<< endl;
      param <<"Name of the Kc file <Kc_NF>:: "<< Kc_NF <<endl;
    }
    else{
      param <<"Kc coefficient <Kccoef>:: "<< Ks_coef <<endl;
      param <<"Name of the Kc file <Kc_NF>:: "<< endl;
    }
    param << endl;
    param << "dtheta, initial water deficit (1=file 2=const_coef) <dtheta_init>:: " << dtheta_init << endl;
    if (dtheta_init==1){
      param <<"dtheta coefficient <dthetacoef>:: "<< endl;
      param <<"Name of the dtheta file <dtheta_NF>:: "<< dtheta_NF <<endl;
    }
    else{
      param <<"dtheta coefficient <dthetacoef>:: "<< dtheta_coef <<endl;
      param <<"Name of the dtheta file <dtheta_NF>:: "<< endl;
    }
    param << endl;

    param << "Psi, load pressure (1=file 2=const_coef) <Psi_init>:: " << Psi_init << endl;
    switch (Psi_init){
      case 1 :
      param <<"Psi coefficient <Psicoef>:: "<< endl;
      param <<"Name of the Psi file <Psi_NF>:: "<< Psi_NF <<endl;
      break;
      case 2 :
      param <<"Psi coefficient <Psicoef>:: "<< Psi_coef <<endl;
      param <<"Name of the Psi file <Psi_NF>:: "<< endl;
      break;
    }
    param << endl;

    param << "imax, maximum infiltration rate (1=file 2=const_coef) <imax_init>:: " << imax_init << endl;
    switch (imax_init){
      case 1 :
      param <<"imax coefficient <imaxcoef>:: "<< endl;
      param <<"Name of the imax file <imax_NF>:: "<< imax_NF <<endl;
      break;
      case 2:
      param <<"imax coefficient <imaxcoef>:: "<< imax_coef <<endl;
      param <<"Name of the imax file <imax_NF>:: "<< endl;
      break;
    }
  }
  param << endl;
  param << "Suffix for the 'Outputs' directory <suffix_o>:: " << suffix_outputs << endl;
  param << endl;
  if (0==nbtimes){
    param << "Format of the Output file (1=gnuplot 2=vtk) <output_f>:: "<< endl;
  }
  else{
    param << "Format of the Output file (1=gnuplot 2=vtk) <output_f>:: " << output_format << endl;
  }
  param << endl;
  param << "Saving specific points (0=no save 1=one point 2=several points) <Choice_points>:: "<< Choice_points << endl;
  switch (Choice_points){
    case 0 :
    param << "x coordinate of the point <x_coord>:: "<< endl;
    param << "y coordinate of the point <y_coord>:: "<< endl;
    param << "Name of file containing the list of the points <list_point_NF>:: "<< endl;
    param << "Choice of times saved (1=all 2=dt value) <Choice_dt_specific_points>:: " << endl;
    param << "Value of time step <dt_specific_points>:: "<< endl;
    break;
    case 1 :
    param << "x coordinate of the point <x_coord>:: "<< x_coord << endl;
    param << "y coordinate of the point <y_coord>:: "<< y_coord << endl;
    param << "Name of file containing the list of the points <list_point_NF>:: "<< endl;
    break;
    case 2 :
    param << "x coordinate of the point <x_coord>:: "<< endl;
    param << "y coordinate of the point <y_coord>:: "<< endl;
    param << "Name of file containing the list of the points <list_point_NF>:: "<< list_point_NF << endl;
    break;
  }
  if (Choice_points>0){
    param << "Choice of times saved (1=all 2=dt value) <Choice_dt_specific_points>:: " << Choice_dt_specific_points << endl;
    switch (Choice_dt_specific_points){
      case 1 :
      param << "Value of time step <dt_specific_points>:: "<< endl;
      break;
      case 2 :
      param << "Value of time step <dt_specific_points>:: "<< dt_specific_points << endl;
      break;
    }
  }

  param.close();

  // print the values in the terminal

  cout << "*********************************************************" << endl;
  cout << "number of grid cells Nxcell = " << Nxcell << endl;
  cout << "number of grid cells Nycell = " << Nycell << endl;
  cout << endl;
  cout << "length of the domain L = "<< L << endl;
  cout << "length of the domain l = "<< l << endl;
  cout << endl;
  cout << "length cells dx = " << dx << endl;
  cout << "length cells dy = " << dy << endl;
  cout << "---------------------------------------------------------" << endl;
  cout << "final time T= " << T << endl;
  cout << "number of times saved= " << nbtimes << endl;

  switch (scheme_type) {
    case 1:
    cout << "choice of type of scheme: fixed cfl" << endl;
    cout << "value of the cfl <cflfix>:: "<< cfl_fix << endl;
    break;
    case 2:
    cout << "choice of type of scheme: fixed dt" << endl;
    cout << "timestep (in seconds) <dtfix>:: "<< dt_fix << endl;
    cout << "value of the cfl <cflfix>:: "<< cfl_fix << endl;
  }

  cout << "---------------------------------------------------------" << endl;

  if(2==Lbound_type) { //In case of constant coefficient
    switch (Lbound[1]) {
      case 1:
      cout << "left condition = imposed h (and q if supercritical) "<< endl;
      cout << " imposed discharge = "<< left_imp_discharge[1] << endl;
      cout << " imposed height = "<< left_imp_h[1] << endl;
      break;
      case 2:
      cout << "left condition = wall " << endl;
      break;
      case 3:
      cout << "left condition = neumann " << endl;
      break;
      case 4:
      cout << "left condition = periodic "<< endl;
      break;
      case 5:
      cout << "left condition = imposed discharge "<< endl;
      cout << " imposed discharge = "<< left_imp_discharge[1] << endl;
      cout << " imposed height = "<< left_imp_h[1] << endl;
    }
  } else {
    cout << "left condition = file reading " << endl;
  }

  if(2==Rbound_type) { //In case of constant coefficient
    switch (Rbound[1]) {
      case 1:
      cout << "right condition = imposed h (and q if supercritical) "<< endl;
      cout << " imposed discharge = "<< right_imp_discharge[1] << endl;
      cout << " imposed height = "<< right_imp_h[1] << endl;
      break;
      case 2:
      cout << "right condition = wall " << endl;
      break;
      case 3:
      cout << "right condition = neumann " << endl;
      break;
      case 4:
      cout << "right condition = periodic "<< endl;
      break;
      case 5:
      cout << "right condition = imposed discharge "<< endl;
      cout << " imposed discharge = "<< right_imp_discharge[1] << endl;
      cout << " imposed height = "<< right_imp_h[1] << endl;
    }
  } else {
    cout << "right condition = file reading " << endl;
  }

  if(2==Bbound_type) { //In case of constant coefficient
    switch (Bbound[1]) {
      case 1:
      cout << "bottom condition = imposed h (and q if supercritical) "<< endl;
      cout << " imposed discharge (inflow supercritical) = "<< bottom_imp_discharge[1] << endl;
      cout << " imposed height = "<< bottom_imp_h[1] << endl;
      break;
      case 2:
      cout << "bottom condition = wall " << endl;
      break;
      case 3:
      cout << "bottom condition = neumann " << endl;
      break;
      case 4:
      cout << "bottom condition = periodic "<< endl;
      break;
      case 5:
      cout << "bottom condition = imposed discharge "<< endl;
      cout << " imposed discharge = "<< bottom_imp_discharge[1] << endl;
      cout << " imposed height = "<< bottom_imp_h[1] << endl;
    }
  } else {
    cout << "bottom condition = file reading " << endl;
  }

  if(2==Tbound_type) { //In case of constant coefficient
    switch (Tbound[1]) {
      case 1:
      cout << "top condition = imposed h (and q if supercritical) "<< endl;
      cout << " imposed discharge (inflow supercritical) = "<< top_imp_discharge[1] << endl;
      cout << " imposed height = "<< top_imp_h[1] << endl;
      break;
      case 2:
      cout << "top condition = wall " << endl;
      break;
      case 3:
      cout << "top condition = neumann " << endl;
      break;
      case 4:
      cout << "top condition = periodic "<< endl;
      break;
      case 5:
      cout << "top condition = imposed discharge "<< endl;
      cout << " imposed discharge = "<< top_imp_discharge[1] << endl;
      cout << " imposed height = "<< top_imp_h[1] << endl;
    }
  } else {
    cout << "top condition = file reading " << endl;
  }
  cout << "---------------------------------------------------------" << endl;

  switch (flux){
    case 1:
    cout << "flux condition = Rusanov"<< endl;
    break;
    case 2:
    cout << "flux condition = HLL "<< endl;
    break;
    case 3:
    cout << "flux condition = HLL2" << endl;
    break;
    case 4:
    cout << "flux condition = HLLC" << endl;
    break;
    case 5:
    cout << "flux condition = HLLC2" << endl;
  }

  cout << "---------------------------------------------------------" << endl;
  cout << "order of the scheme = " <<order << endl;
  cout << "---------------------------------------------------------" << endl;

  switch (fric){
    case 0:
    cout << "no friction"<< endl;
    break;
    case 1:
    cout << "friction condition = Manning"<< endl;
    if (fric_init==1){
      cout << "parameter friction = file reading"<< endl;
    }
    else{
      cout << "parameter friction = "<< friccoef << endl;
    }
    break;
    case 2:
    cout << "friction condition = Darcy-Weisbach"<< endl;
    if (fric_init==1){
      cout << "parameter friction = file reading" << endl;
    }
    else{
      cout << "parameter friction = "<< friccoef << endl;
    }
    break;
    case 3:
    cout << "friction law = laminar"<< endl;
    if (fric_init==1){
      cout << "parameter friction = file reading"<< endl;
    }
    else{
      cout << "parameter friction = "<< friccoef << endl;
    }
  }
  cout << "---------------------------------------------------------" << endl;

  if (2 == order){// in case of order 2
    switch (rec){
      case 1:
      cout << "reconstruction condition = MUSCL" << endl;
      break;
      case 2:
      cout << "reconstruction condition = ENO " << endl;
      cout << "parameter amortENO = "<< amortENO << endl;
      break;
      case 3:
      cout << "reconstruction condition = ENO_mod "<< endl;
      cout << "parameter amortENO = "<< amortENO << endl;
      cout << "parameter modifENO = "<< modifENO<< endl;
    }

    switch (lim){
      case 1:
      cout << "limitation condition = Minmod " << endl;
      break;
      case 2:
      cout << "limitation condition = VanAlbada " << endl;
      break;
      case 3:
      cout << "limitation condition = VanLeer " << endl;
    }

    cout << "---------------------------------------------------------" << endl;
  }

  switch (inf){
    case 0:
    cout << "no infiltration"<< endl;
    break;
    case 1:
    cout << "infiltration model = Green-Ampt"<< endl;
    if (Kc_init==1){
      cout << "parameter Kc     = file reading "<< endl;
    }
    else{
      cout << "parameter Kc     = " <<   Kc_coef << endl;
    }
    if (Ks_init==1){
      cout << "parameter Ks     = file reading "<< endl;
    }
    else{
      cout << "parameter Ks     = " <<   Ks_coef << endl;
    }
    if (dtheta_init==1){
      cout << "parameter dtheta = file reading " << endl;
    }
    else{
      cout << "parameter dtheta = " <<   dtheta_coef << endl;
    }
    if (Psi_init==1){
      cout << "parameter Psi    = file reading "<< endl;
    }
    else{
      cout << "parameter Psi    = " <<   Psi_coef << endl;
    }
    if (zcrust_init==1){
      cout << "parameter zcrust = file reading "  << endl;
    }
    else{
      cout << "parameter zcrust = " <<   zcrust_coef << endl;
    }
    if (imax_init==1){
      cout << "parameter imax   = file reading " << endl;
    }
    else{
      cout << "parameter imax   = " <<  imax_coef << endl;
    }

    break;
  }

  cout << "---------------------------------------------------------" << endl;
  switch (huv_init){
    case 1:
    cout << "hu initial condition = file reading" << endl;
    break;
    case 2:
    cout << "hu initial condition = h, u & v = 0 " << endl;
    break;
    case 3:
    cout << "hu initial condition = function Thacker" << endl;
    break;
    case 4:
    cout << "hu initial condition = function Radial_Dam_dry" << endl;
    break;
    case 5:
    cout << "hu initial condition = function Radial_Dam_wet" << endl;

  }
  cout << "---------------------------------------------------------" << endl;

  switch (topo){
    case 1:
    cout << "topography = file reading" << endl;
    break;
    case 2:
    cout << "topography = function generated" << endl;
    break;
    case 3:
    cout << "topography = function Thacker" << endl;
    break;
  }
  cout << "---------------------------------------------------------" << endl;
  switch (rain){
    case 0:
    cout << "no rain"<< endl;
    break;
    case 1:
    cout << "rain condition = file reading" << endl;
    break;
    case 2:
    cout << "rain condition = rain intensity = 0.00001 m/s " << endl;
  }

  switch (output_format){
    case 1:
    cout << "---------------------------------------------------------" << endl;
    cout << "format of the output files = gnuplot " << endl;
    break;
    case 2:
    cout << "---------------------------------------------------------" << endl;
    cout << "format of the output files = vtk " << endl;
    break;
  }

  switch (Choice_points){
    case 0:
    cout << "---------------------------------------------------------" << endl;
    cout << "no specific point in saving"<< endl;
    break;
    case 1:
    cout << "---------------------------------------------------------" << endl;
    cout << "x coordinate of the specific point in saving "<< x_coord << endl;
    cout << "y coordinate of the specific point in saving "<< y_coord << endl;
    break;
    case 2:
    cout << "---------------------------------------------------------" << endl;
    cout << "Saving specific points = file reading " << endl;
    break;
  }

  if (Choice_points>0){
  	switch (Choice_dt_specific_points){
    	case 1:
    	cout << "---------------------------------------------------------" << endl;
    	cout << "Choice of times saved = all "<< endl;
    	break;
    	case 2:
    	cout << "---------------------------------------------------------" << endl;
    	cout << "Choice of times saved =  "<< Choice_dt_specific_points << endl;
    	cout << "Value of time step =  "<< dt_specific_points << endl;
    	break;
  }
  }

  cout << "---------------------------------------------------------" << endl;
  cout << "entries ok" << endl;
  cout << "*********************************************************" << endl;

  /*
   The value of the imposed discharge per cell in the boundary condition is
   in m2/s
   */
  if(2==Lbound_type) { //In case of constant coefficient
    for (int j=1 ; j<=Nycell ; j++)
    left_imp_discharge[j] /= l;
  }else {
    for (int j=1 ; j<=Nycell ; j++)
    left_imp_discharge[j] /= dy;
  }

  if(2==Rbound_type) { //In case of constant coefficient
    for (int j=1 ; j<=Nycell ; j++)
    right_imp_discharge[j] /= l;
  }else {
    for (int j=1 ; j<=Nycell ; j++)
    right_imp_discharge[j] /= dy;
  }

  if(2==Bbound_type) { //In case of constant coefficient
    for (int i=1 ; i<=Nxcell ; i++)
    bottom_imp_discharge[i] /= L;
  }else {
    for (int i=1 ; i<=Nxcell ; i++)
    bottom_imp_discharge[i] /= dx;
  }

  if(2==Tbound_type) { //In case of constant coefficient
    for (int i=1 ; i<=Nxcell ; i++)
    top_imp_discharge[i] /= L;
  }else {
    for (int i=1 ; i<=Nxcell ; i++)
    top_imp_discharge[i] /= dx;
  }



}


int Parameters::get_Nxcell() const {

  /**
   * @details
   * @return The number of cells in space in the first (x) direction Parameters#Nxcell.
   */

  return Nxcell;
}

int Parameters::get_Nycell() const {

  /**
   * @details
   * @return The number of cells in space in the second (y) direction Parameters#Nycell.
   */

  return Nycell;
}

SCALAR Parameters::get_T() const {

  /**
   * @details
   * @return The final time Parameters#T.
   */

  return T;
}

int Parameters::get_nbtimes() const {

  /**
   * @details
   * @return The number of times saved Parameters#nbtimes.
   */

  return nbtimes;
}

int Parameters::get_scheme_type() const {

  /**
   * @details
   * @return The type of scheme Parameters#scheme_type.
   */

  return scheme_type;
}

SCALAR Parameters::get_dtfix() const {

  /**
   * @details
   * @return The fixed space step Parameters#dx_fix.
   */

  return dt_fix;
}

SCALAR Parameters::get_cflfix() const {

  /**
   * @details
   * @return The fixed cfl Parameters#cfl_fix.
   */

  return cfl_fix;
}

SCALAR Parameters::get_dx() const {

  /**
   * @details
   * @return The space step in the first (x) direction Parameters#dx.
   */

  return dx;
}

SCALAR Parameters::get_dy() const {

  /**
   * @details
   * @return The space step in the second (y) direction Parameters#dy.
   */

  return dy;
}

SCALAR Parameters::get_L() const {

  /**
   * @details
   * @return the length of the domain in the first (x) direction.
   */

  return L;
}

SCALAR Parameters::get_l() const {

  /**
   * @details
   * @return the length of the domain in the second (y) direction.
   */

  return l;
}

int Parameters::get_type_Lbound()  const  {
  /**
   * @details
   * @return The value corresponding to the left boundary condition Parameters#Lbound.
   */

  return Lbound_type;
}



map<int,int> & Parameters::get_Lbound()  {

  /**
   * @details
   * @return The value corresponding to the left boundary condition Parameters#Lbound.
   */

  return Lbound;
}

map<SCALAR,string> & Parameters::get_times_files_Lbound()  {

  /**
   * @details
   * @return The value corresponding to the bottom boundary condition Parameters#Bbound.
   */

  return left_times_files;
}
//SCALAR Parameters::get_left_imp_discharge() const {
map<int,SCALAR> & Parameters::get_left_imp_discharge()  {
  /**
   * @details
   * @return The value of the imposed discharge per cell in the left boundary condition, that is Parameters#left_imp_discharge / Parameters#l.
   */


  return left_imp_discharge;
}

//SCALAR Parameters::get_left_imp_h() const {
map<int,SCALAR> & Parameters::get_left_imp_h()  {
  /**
   * @details
   * @return The value of the imposed water height in the left boundary condition Parameters#left_imp_h.
   */

  return left_imp_h;
}

int Parameters::get_type_Rbound() const {
  /**
   * @details
   * @return The value corresponding to the left boundary condition Parameters#Lbound.
   */

  return Rbound_type;
}

map<int,int> & Parameters::get_Rbound()  {

  /**
   * @details
   * @return The value corresponding to the right boundary condition Parameters#Rbound.
   */

  return Rbound;
}

map<SCALAR,string> & Parameters::get_times_files_Rbound()  {

  /**
   * @details
   * @return The value corresponding to the bottom boundary condition Parameters#Bbound.
   */

  return right_times_files;
}
//SCALAR Parameters::get_right_imp_discharge() const {
map<int,SCALAR> & Parameters::get_right_imp_discharge()  {
  /**
   * @details
   * @return The value of the imposed discharge per cell in the right boundary condition, that is Parameters#right_imp_discharge / Parameters#l.
   */

  return right_imp_discharge;
}

//SCALAR Parameters::get_right_imp_h() const {
map<int,SCALAR> & Parameters::get_right_imp_h()  {
  /**
   * @details
   * @return The value of the imposed water height in the right boundary condition Parameters#right_imp_h.
   */

  return right_imp_h;
}

int Parameters::get_type_Bbound() const {
  /**
   * @details
   * @return The value corresponding to the left boundary condition Parameters#Lbound.
   */

  return Bbound_type;
}



map<int,int> & Parameters::get_Bbound() {

  /**
   * @details
   * @return The value corresponding to the bottom boundary condition Parameters#Bbound.
   */

  return Bbound;
}

map<SCALAR,string> & Parameters::get_times_files_Bbound()  {

  /**
   * @details
   * @return The value corresponding to the bottom boundary condition Parameters#Bbound.
   */

  return bottom_times_files;
}



//SCALAR Parameters::get_bottom_imp_discharge() const {
map<int,SCALAR> & Parameters::get_bottom_imp_discharge()  {
  /**
   * @details
   * @return The value of the imposed discharge per cell in the bottom boundary condition, that is Parameters#bottom_imp_discharge / Parameters#L.
   */

  return bottom_imp_discharge;
}

//SCALAR Parameters::get_bottom_imp_h() const {
map<int,SCALAR> & Parameters::get_bottom_imp_h()  {
  /**
   * @details
   * @return The value of the imposed water height in the bottom boundary condition Parameters#bottom_imp_h.
   */

  return bottom_imp_h;
}

int Parameters::get_type_Tbound() const {
  /**
   * @details
   * @return The value corresponding to the left boundary condition Parameters#Lbound.
   */

  return Tbound_type;
}



map<int,int> & Parameters::get_Tbound()  {

  /**
   * @details
   * @return The value corresponding to the top boundary condition Parameters#Tbound.
   */

  return Tbound;
}

map<SCALAR,string> & Parameters::get_times_files_Tbound()  {

  /**
   * @details
   * @return The value corresponding to the bottom boundary condition Parameters#Bbound.
   */

  return top_times_files;
}

//SCALAR Parameters::get_top_imp_discharge() const {
map<int,SCALAR> & Parameters::get_top_imp_discharge()  {
  /**
   * @details
   * @return The value of the imposed discharge per cell in the top boundary condition, that is Parameters#top_imp_discharge / Parameters#L.
   */

  return top_imp_discharge;
}

//SCALAR Parameters::get_top_imp_h() const {
map<int,SCALAR> & Parameters::get_top_imp_h()   {
  /**
   * @details
   * @return The value of the imposed water height in the bottom boundary condition Parameters#top_imp_h.
   */

  return top_imp_h;
}

int Parameters::get_flux() const {

  /**
   * @details
   * @return The value corresponding to the flux Parameters#flux.
   */

  return flux;
}

int Parameters::get_order() const {

  /**
   * @details
   * @return The order of the scheme Parameters#order.
   */

  return order;
}

int Parameters::get_rec() const {

  /**
   * @details
   * @return The value corresponding to the reconstruction Parameters#rec.
   */

  return rec;
}

int Parameters::get_fric() const {

  /**
   * @details
   * @return The value corresponding to the friction law Parameters#fric.
   */

  return fric;
}


int Parameters::get_lim() const {

  /**
   * @details
   * @return The value corresponding to the limiter Parameters#lim.
   */

  return lim;
}

SCALAR Parameters::get_amortENO() const {

  /**
   * @details
   * @return The value of the amortENO parameter Parameters#amortENO.
   */

  return amortENO;
}

SCALAR Parameters::get_modifENO() const {

  /**
   * @details
   * @return The value of the modifENO parameter Parameters#modifENO.
   */

  return modifENO;
}

int Parameters::get_inf() const {

  /**
   * @details
   * @return The value corresponding to the infiltration Parameters#inf.
   */

  return inf;
}

int Parameters::get_fric_init() const {

  /**
   * @details
   * @return The value corresponding to the friction coefficient Parameters#fric_init.
   */

  return fric_init;
}

string Parameters::get_frictionNameFile(void) const{

  /**
   * @details
   * @return The friction coefficient path + Input directory Parameters#fric_namefile.
   */

  return fric_namefile;
}

string Parameters::get_frictionNameFileS(void) const{

  /**
   * @details
   * @return The friction coefficient namefile (inside the Input directory) Parameters#fric_NF.
   */

  return fric_NF;
}

SCALAR Parameters::get_friccoef() const {

  /**
   * @details
   * @return The value of the friction coefficient Parameters#friccoef.
   */

  return friccoef;
}

string Parameters::get_topographyNameFile(void) const {

  /**
   * @details
   * @return The topography path + Input directory Parameters#topography_namefile.
   */

  return topography_namefile;
}

string Parameters::get_topographyNameFileS(void) const {

  /**
   * @details
   * @return The topography namefile (inside the Input directory) Parameters#topo_NF.
   */

  return topo_NF;
}

int Parameters::get_topo() const {

  /**
   * @details
   * @return The value corresponding to the topography Parameters#topo.
   */

  return topo;
}

int Parameters::get_huv() const {

  /**
   * @details
   * @return The value corresponding to the initialization of h and u,v Parameters#huv_init.
   */

  return huv_init;
}

string Parameters::get_huvNameFile(void) const {

  /**
   * @details
   * @return The h and u,v path for the initialization + Input directory Parameters#huv_namefile.
   */

  return huv_namefile;
}

string Parameters::get_huvNameFileS(void) const {

  /**
   * @details
   * @return The h and u namefile for the initialization (inside the Input directory) Parameters#huv_NF.
   */

  return huv_NF;
}

int Parameters::get_rain() const {

  /**
   * @details
   * @return The value corresponding to the initialization of the rain Parameters#rain.
   */

  return rain;
}

string Parameters::get_rainNameFile(void) const {

  /**
   * @details
   * @return The rain path for the initialization + Input directory Parameters#rain_namefile.
   */

  return rain_namefile;
}

string Parameters::get_rainNameFileS(void) const {

  /**
   * @details
   * @return The rain namefile for the initialization (inside the Input directory) Parameters#rain_NF.
   */

  return rain_NF;
}

string Parameters::get_outputDirectory(void) const {

  /**
   * @details
   * @return The output directory with the suffix Parameters#output_directory.
   */

  return output_directory;
}

int Parameters::get_Kc_init() const {

  /**
   * @details
   * @return The value corresponding to the initialization of Kc Parameters#Kc_init.
   */

  return Kc_init;
}

SCALAR Parameters::get_Kc_coef() const {

  /**
   * @details
   * @return The value of Kc Parameters#Kc_coef.
   */

  return Kc_coef;
}

string Parameters::get_KcNameFile(void) const {

  /**
   * @details
   * @return The Kc path for the initialization + Input directory Parameters#Kc_namefile.
   */

  return Kc_namefile;
}

string Parameters::get_KcNameFileS() const {

  /**
   * @details
   * @return The Kc namefile for the initialization (inside the Input directory) Parameters#Kc_NF.
   */

  return Kc_NF;
}

int Parameters::get_Ks_init() const {

  /**
   * @details
   * @return The value corresponding to the initialization of Ks Parameters#Ks_init.
   */

  return Ks_init;
}

SCALAR Parameters::get_Ks_coef() const {

  /**
   * @details
   * @return The value of Ks Parameters#Ks_coef.
   */

  return Ks_coef;
}

string Parameters::get_KsNameFile(void) const {

  /**
   * @details
   * @return The Ks path for the initialization + Input directory Parameters#Ks_namefile.
   */

  return Ks_namefile;
}

string Parameters::get_KsNameFileS() const {

  /**
   * @details
   * @return The Ks namefile for the initialization (inside the Input directory) Parameters#Ks_NF.
   */

  return Ks_NF;
}

int Parameters::get_dtheta_init() const {

  /**
   * @details
   * @return The value corresponding to the initialization of dtheta Parameters#dtheta_init.
   */

  return dtheta_init;
}

SCALAR Parameters::get_dtheta_coef() const {

  /**
   * @details
   * @return The value of dtheta Parameters#dtheta_coef.
   */

  return dtheta_coef;
}

string Parameters::get_dthetaNameFile(void) const {

  /**
   * @details
   * @return The dtheta path for the initialization + Input directory Parameters#dtheta_namefile.
   */

  return dtheta_namefile;
}

string Parameters::get_dthetaNameFileS() const {

  /**
   * @details
   * @return The dtheta namefile for the initialization (inside the Input directory) Parameters#dtheta_NF.
   */

  return dtheta_NF;
}

int Parameters::get_Psi_init() const {

  /**
   * @details
   * @return The value corresponding to the initialization of Psi Parameters#Psi_init.
   */

  return Psi_init;
}

SCALAR Parameters::get_Psi_coef() const {

  /**
   * @details
   * @return The value of Psi Parameters#Psi_coef.
   */

  return Psi_coef;
}

string Parameters::get_PsiNameFile(void) const {

  /**
   * @details
   * @return The Psi path for the initialization + Input directory Parameters#Psi_namefile.
   */

  return Psi_namefile;
}

string Parameters::get_PsiNameFileS() const {

  /**
   * @details
   * @return The Psi namefile for the initialization (inside the Input directory) Parameters#Psi_NF.
   */

  return Psi_NF;
}

int Parameters::get_zcrust_init() const {

  /**
   * @details
   * @return The value corresponding to the initialization of zcrust Parameters#zcrust_init.
   */

  return zcrust_init;
}

SCALAR Parameters::get_zcrust_coef() const {

  /**
   * @details
   * @return The value of zcrust Parameters#zcrust_coef.
   */

  return zcrust_coef;
}

string Parameters::get_zcrustNameFile(void) const {

  /**
   * @details
   * @return The zcrust path for the initialization + Input directory Parameters#zcrust_namefile.
   */

  return zcrust_namefile;
}

string Parameters::get_zcrustNameFileS() const {

  /**
   * @details
   * @return The zcrust namefile for the initialization (inside the Input directory) Parameters#zcrust_NF.
   */

  return zcrust_NF;
}

int Parameters::get_imax_init() const {

  /**
   * @details
   * @return The value corresponding to the initialization of imax Parameters#imax_init.
   */

  return imax_init;
}

SCALAR Parameters::get_imax_coef() const {

  /**
   * @details
   * @return The value of imax Parameters#imax_coef.
   */

  return imax_coef;
}

string Parameters::get_imaxNameFile(void) const {

  /**
   * @details
   * @return The imax path for the initialization + Input directory Parameters#imax_namefile.
   */

  return imax_namefile;
}

string Parameters::get_imaxNameFileS() const {

  /**
   * @details
   * @return The imax namefile for the initialization (inside the Input directory) Parameters#imax_NF.
   */

  return imax_NF;
}

string Parameters::get_suffix(void) const{

  /**
   * @details
   * @return The suffix (for the output directory) Parameters#suffix_outputs.
   */

  return suffix_outputs;
}

int Parameters::get_output() const {

  /**
   * @details
   * @return The type of output Parameters#output_format.
   */

  return output_format;
}



void Parameters::fill_array(TAB & myarray,const SCALAR myvalue) const {

  /**
   * @details Fills an array with a constant value.
   * @param[in, out] myarray array to fill.
   * @param[in] myvalue value.
   */

  for (unsigned int i=0;i<myarray.size();i++){
    //fill is a function provided by the STL that initializes a container with a value.
    fill(myarray[i].begin(),myarray[i].end(), myvalue);
  } //end for i
}

void Parameters::fill_array(TAB & myarray,string namefile) const {

  /**
   * @details Fills an array with the values given in the file
   * @param[in, out] myarray array to fill.
   * @param[in] namefile name of the file containing the values to be inserted into the array.
   * @warning ***: ERROR: cannot open the file.
   * @warning ***: ERROR: the number of data in this file is too big/small.
   * @warning ***: ERROR: line ***.
   * @warning ***: ERROR: the value for the point x = *** y = *** is missing.
   * @warning ***: WARNING: line *** ; a commentary should begin with the # symbol.
   * @note If the array cannot be filled correctly, the code will exit with failure termination code.
   */

  SCALAR x, y,value;
  int row, column;
  int iter=0; //iterator
  int line_number=0; //iterator
  string line;// string to store a line of the input file
  char car;//First character of a commentary line

  //Opening the file and verification if it exists.
  ifstream getdata(namefile.c_str(),ios::in);
  if (!getdata){
    //The input file cannot be opened.
    cerr << namefile  << ": ERROR: cannot open the file\n";
    exit(EXIT_FAILURE);
  }

  //Initialization of myarray to the largest finite representable floating-point number.
  //So that we will be able to check that the user fills the variable properly.
  for (int i=1 ; i<Nxcell+1 ; i++){
    for (int j=1 ; j<Nycell+1 ; j++){
      myarray[i][j]=MAX_SCAL;
    } //end for j
  } //end for i

  // As long as the end of the file is not reached, the next line is read.
  while (!getdata.eof()){
    line_number++;
    getline (getdata, line, '\n');       // read a line
    istringstream entree (line);
    if (entree >> x >> y >> value) {

      //Check that the input file contains the expected number of data.
      if(++iter > Nxcell*Nycell){
        cerr << namefile << ": ERROR: the number of data in this file is too big!"<<endl;
        exit(EXIT_FAILURE);
      }
      //We compute the index of the array from the space variable in the input file.
      row=(int)ceil(x/dx); column=(int)ceil(y/dy); //The index corresponds to the smallest integer superior or equal to x/dx or y/dy.

      //Error if the x or y value of the input file is out of the domain.
      if(x < 0){
        cerr<< namefile <<": ERROR: line "<< line_number<<", x = "<< x <<" must be positive."<< endl;
        exit(EXIT_FAILURE);
      }
      if(y < 0){
        cerr<< namefile <<": ERROR: line "<< line_number<<", y = "<< y <<" must be positive."<< endl;
        exit(EXIT_FAILURE);
      }
      if(x > Nxcell*dx){
        cerr<< namefile <<": ERROR: line "<< line_number<<", x = "<< x <<" must be lower than "<< Nxcell*dx <<"."<< endl;
        exit(EXIT_FAILURE);
      }
      if(y > Nycell*dy){
        cerr<< namefile <<": ERROR: at line "<< line_number<<", y = "<< y <<" must be lower than "<< Nycell*dy <<"."<< endl;
        exit(EXIT_FAILURE);
      }

      //Error if the x or y value of the input file is not near the center of the cell.
      if (fabs(x - (row-0.5)*dx) > dx*RATIO_CLOSE_CELL){
        cout<< namefile <<": ERROR: line "<< line_number<<"; x = "<< x <<";\n This value is not close enough to the center of the cell.\n You may want to replace it with " << (row-0.5)*dx <<endl;
        exit(EXIT_FAILURE);
      }
      if (fabs(y - (column-0.5)*dy) > dy*RATIO_CLOSE_CELL){
        cout<<"column ="<<column<<endl;
        cout<< namefile <<": ERROR: line "<< line_number<<"; y = "<< y <<";\n This value is not close enough to the center of the cell.\n You may want to replace it with " << (column-0.5)*dy <<endl;
        exit(EXIT_FAILURE);
      }

      //Store the input values into the array.
      myarray[row][column]=value;
    }
    else{
      car='#'; //Initialization of the character used to identity the beginning of a comment
      istringstream entree_car (line);
      entree_car >> car;
      if (car != '#'){
        cout << namefile << ": WARNING: line "<< line_number <<"; a commentary should begin with the # symbol " <<endl;
      }
    }
  }

  //Closing the input file
  getdata.close();

  //Check that the input file contains the expected number of data.
  if(iter < Nxcell*Nycell){
    cerr << namefile<< ": ERROR: the number of data in this file is too small!"<<endl;
    exit(EXIT_FAILURE);
  }

  //Final check: Does all the grid cells were filled with a value?
  for (int i=1 ; i<Nxcell+1 ; i++){
    for (int j=1 ; j<Nycell+1 ; j++){
      if(myarray[i][j]>=MAX_SCAL){
        cerr<< namefile<< ": ERROR: the value for the point x ="<< (i-0.5)*dx <<" y = "<< (j-0.5)*dy <<" is missing!"<<endl;
        exit(EXIT_FAILURE);
      }
    }//end for j
  } //end for i

}

bool Parameters::is_coord_in_file_valid(const SCALAR & x, const SCALAR & y, int line_number,string namefile) const{
  /**
   * @details Fills an array with the values given in the file
   * @param[in] x coordinate of the specific point.
   * @param[in] y coordinate of the specific point.
   * @param[in] line_number the line number where the error appears.
   * @param[in] namefile name of the file containing the values to be verified.
   * @warning ***: ERROR: .
   * @warning ***: ERROR: x / y must be positive / lower than ***.
   * @warning ***: ERROR: line ***.
   * @note -1 is used to not display the line number in the error message .
   * @note if the coordinates are not valid, the code will return with false.
   */

  //Error if the x or y value of the input file is out of the domain.
  if(x < 0){
    if (-1 == line_number){
      cerr<< namefile <<": ERROR: x = "<< x <<" must be positive."<< endl;
    }
    else {
      cerr<< namefile <<": ERROR: line "<< line_number<<", x = "<< x <<" must be positive."<< endl;
    }
    return false;
  }
  if(y < 0){
    if (-1 == line_number){
      cerr<< namefile <<": ERROR: y = "<< y <<" must be positive."<< endl;
    }
    else {
      cerr<< namefile <<": ERROR: line "<< line_number<<", y = "<< y <<" must be positive."<< endl;
    }
    return false;
  }

  if(x > Nxcell*dx){
    if (-1 == line_number){
      cout << "Nxcell = " << Nxcell << "  , dx = "<< dx << endl;
      cerr<< namefile <<": ERROR:  x = "<< x <<" must be lower than "<< Nxcell*dx <<"."<< endl;
    }
    else {
      cerr<< namefile <<": ERROR: line "<< line_number<<", x = "<< x <<" must be lower than "<< Nxcell*dx <<"."<< endl;
    }
    return false;
  }

  if(y > Nycell*dy){
    if (-1 == line_number){
      cerr<< namefile <<": ERROR: y = "<< y <<" must be lower than "<< Nycell*dy <<"."<< endl;

    }
    else {
      cerr<< namefile <<": ERROR: at line "<< line_number<<", y = "<< y <<" must be lower than "<< Nycell*dy <<"."<< endl;
    }
    return false;
  }

  //We compute the index of the array from the space variable in the input file.
  int row=(int)ceil(x/dx); //The index corresponds to the smallest integer superior or equal to x/dx
  int column=(int)ceil(y/dy);//The index corresponds to the smallest integer superior or equal to y/dy.

  //Error if the x or y value of the input file is not near the center of the cell.
  if (fabs(x - (row-0.5)*dx) > dx*RATIO_CLOSE_CELL){
    if (-1 == line_number){
      cout<< namefile <<": ERROR: x = "<< x <<";\n This value is not close enough to the center of the cell.\n You may want to replace it with " << fabs((row-0.5)*dx) <<endl;
    }
    else {
      cout<< namefile <<": ERROR: line "<< line_number<<"; x = "<< x <<";\n This value is not close enough to the center of the cell.\n You may want to replace it with " << fabs((row-0.5)*dx) <<endl;
    }
    return false;
  }
  if (fabs(y - (column-0.5)*dy) > dy*RATIO_CLOSE_CELL){
    if (-1 == line_number){
      cout<<"column ="<<column<<endl;
      cout<< namefile <<": ERROR: y = "<< y <<";\n This value is not close enough to the center of the cell.\n You may want to replace it with " << fabs((column-0.5)*dy) <<endl;
    }
    else {
      cout<<"column ="<<column<<endl;
      cout<< namefile <<": ERROR: line "<< line_number<<"; y = "<< y <<";\n This value is not close enough to the center of the cell.\n You may want to replace it with " << fabs((column-0.5)*dy) <<endl;
    }
    return false;
  }

  return true;
}

SCALAR Parameters::get_x_coord(void) const{

  /**
   * @details
   * @return x coordinate of the specific point to be saved.
   */

  return x_coord;
}

SCALAR Parameters::get_y_coord(void) const{

  /**
   * @details
   * @return y coordinate of the specific point to be saved.
   */

  return y_coord;
}

string Parameters::get_list_pointNameFile(void) const {

  /**
   * @details
   * @return The name of file containing the list of the specific points path for the initialization + Input directory Parameters#list_point_namefile.
   */

  return list_point_namefile;
}
string Parameters::get_list_pointNameFileS(void) const {

  /**
   * @details
   * @return The name of file containing list of the specific points for the initialization (inside the Input directory) Parameters#list_point_NF.
   */

  return list_point_NF;
}

int Parameters::get_choice_specific_point() const {

  /**
   * @details
   * @return The choice for saving specific points.
   */

  return Choice_points;
}

int Parameters::get_choice_dt_specific_points() const {

  /**
   * @details
   * @return The choice of times saved for the specific points.
   */

  return Choice_dt_specific_points;
}

SCALAR Parameters::get_dt_specific_points() const {

  /**
   * @details
   * @return The time step for saving specific points.
   */

  return dt_specific_points;
}

map<int,int>  Parameters::fill_array_bc_inhomogeneous(string Bc_NF, string path_input_directory,char BC,map<int,SCALAR> & imp_q, map<int,SCALAR> & imp_h){


  /**
   * @details Extracts from the Bc_NF file the type of boundary condition,
   *          discharge and water heights imposed for each point.
   * @param[in] Bc_NF name of the file containing the values to be verified.
   * @param[in] path_input_directory path of directory containing Bc_NF file.
   * @param[in] BC represents the boundary condition which we handle
   * @param[out] imp_q container containing discharges [m3/s]
   * @param[out] imp_h container containing water heights [m]
   * @warning ***: ERROR: cannot open the file.
   * @warning ***: ERROR: it is necessary to specify a file for this choice of initialization of boundary condition.
   * @warning ***: ERROR: the file does not exist in the directory Inputs.
   * @warning ***: ERROR: the number of data in this file is not good.
   * @warning ***: ERROR: the format of the file is not good.
   * @warning ***: ERROR: line ***.
   * @warning ***: ERROR: the boundary does not exist!
   * @return the vector containing the choice of boundary conditions.
   * @note If the containers cannot be filled correctly, the code will exit with failure termination code.
   */

  int number_of_points = 0;
  int line_number = 0;
  /* vtab contains the elements of line of Bc_NF file:
   if the line is valid, then:
   vtab[0] contains the x coordinate of the boundary condition
   vtab[1] contains the type of boundary condition (ex: 1=imp.h 2=wall 3=neumann 4=periodic 5=imp.q)
   vtab[2] contains discharge if vtab[1] = 1 or vtab[1]= 5
   vtab[3] contains water height if vtab[1] = 1 or vtab[1]= 5

   */
  vector<double> vtab;

  /* vchoice contains the type of boundary condition.
   vchoice[0] contains the index of x coordinate of the boundary condition: vchoice[0] = (int) vtab[0]/dx
   vchoice[1] contains the type of boundary condition: vchoice[1] =  vtab[1].
   */
  map<int,int> vchoice;



  int choice_bc; // represents the index of x coordinate: choice_bc is between 1 and (Nxcell or Nycell).

  /* Initialisation */
  imp_q.clear();
  imp_h.clear();

  if (Bc_NF.size()<1){  // if file's name is empty
    cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of boundary condition ."<<endl;
    exit(EXIT_FAILURE);
  }

  /* if file is exist */
  if (-1 == access ((path_input_directory+Bc_NF).c_str(), R_OK) ) {
    cerr << " parameters.txt: ERROR: the "<< Bc_NF<< " file does not exist in the directory Inputs." << endl ;
    exit(EXIT_FAILURE);
  }

  ifstream inputfile((path_input_directory+Bc_NF).c_str());
  if (!inputfile){
    cout << " parameters.txt : ERROR: cannot open "<< path_input_directory+Bc_NF << " file\n";
    exit(EXIT_FAILURE);
  }


  int index_x;// represents the index of x coordinate
  double number;
  char car;

  string line;

  while(getline(inputfile,line)){
    ++line_number;

    istringstream ss(line);
    ss>>car;

    if(isdigit(car)){ //Verify if the first caracter of the line is a digit

      ss.putback(car); //Put character back in the line to have a complete line.

      while (ss>>number) //fills vtab with all elements of the line
      vtab.push_back(number);

      if(vtab.size()<2){ // Each line must contain at least two elements (x coordinate and
                         // choice of boundary condition)
        cerr << " parameters.txt: ERROR: the format of "<< Bc_NF<< " file is not good." << endl ;
        exit(EXIT_FAILURE);
      }

      choice_bc = int(vtab[1]);

      switch(BC){
        case 'B':
        case 'T':
        if(!is_coord_in_file_valid(vtab[0],0.5*dy,line_number,(path_input_directory+Bc_NF).c_str())){
          exit(EXIT_FAILURE);
        }
        index_x=(int)ceil(vtab[0]/dx);
        break;
        case 'L':
        case 'R':
        if(!is_coord_in_file_valid(0.5*dx,vtab[0],line_number,(path_input_directory+Bc_NF).c_str())){
          exit(EXIT_FAILURE);
        }
        index_x=(int)ceil(vtab[0]/dy);
        break;
        default:
        cerr << "fill_array_bc_inhomogeneous: ERROR: the boundary " << BC << "does not exist!" << endl;
        exit(EXIT_FAILURE);
      }

      vchoice[index_x]  = choice_bc;
      imp_q[index_x]    = MAX_SCAL;
      imp_h[index_x]    = MAX_SCAL;

      if((1==choice_bc)||(5==choice_bc)){ // in case of imp.h or imp.q it's necessary to have q and h
        if(vtab.size()>3){
          imp_q[index_x]=vtab[2];

          if (vtab[3]< 0.){
            cerr << path_input_directory+Bc_NF <<": ERROR: line "<< line_number << ", negative imposed water height in boundary condition."<< endl;
            exit(EXIT_FAILURE);
          }

          imp_h[index_x]=vtab[3];

        }else{
          cerr << " parameters.txt: ERROR: the format of "<< Bc_NF<< " file is not good." << endl ;
          exit(EXIT_FAILURE);
        }
      }// end in case of imp.h or imp.q

      ++number_of_points;
      vtab.clear();
    }//end if digit(car)
  }//end while


  switch(BC){
    case 'B':
    case 'T':
    if(Nxcell != number_of_points){
      cerr << " parameters.txt: ERROR: the number of data in << "<< Bc_NF<< " >> file is not good." << endl ;
      cerr << " total number is " << number_of_points << " instead of " << Nxcell << endl;
      exit(EXIT_FAILURE);
    }
    break;
    case 'L':
    case 'R':
    if(Nycell != number_of_points){
      cerr << " parameters.txt: ERROR: the number of data in << "<< Bc_NF<< " >> file is not good." << endl ;
      cerr << " total number is " << number_of_points << " instead of " << Nycell << endl;
      exit(EXIT_FAILURE);
    }
    break;
  }

  return vchoice;
  /*----------------------------------------------------------:--------------------------------------------------*/
}

map<double,string>  Parameters::verif_file_bc_inhomogeneous(string Bc_NF, string path_input_directory,char BC,
                                                            map<int,int> &  vchoice,map<int,SCALAR> & imp_q,
                                                            map<int,SCALAR> & imp_h)
{


  /**
   * @details Extracts from the Bc_NF file the list of files and time values corresponding to the boundary condition.
   * @param[in] Bc_NF name of the file containing the values to be verified.
   * @param[in] path_input_directory path of directory containing Bc_NF file.
   * @param[in] BC represents the boundary condition which we handle
   * @param[out] vchoice container containing the choice of boundary conditions at time equal to 0s.
   * @param[out] imp_q container containing discharges [m3/s] at time equal to 0s.
   * @param[out] imp_h container containing water heights [m] at time equal to 0s.
   * @warning ***: ERROR: cannot open the file.
   * @warning ***: ERROR: line ***.
   * @warning ***: ERROR: the first time must be t = 0.
   * @warning ***: ERROR: it is necessary to specify a file for this choice of initialization of boundary condition.
   * @warning ***: ERROR: the file does not exist in the directory Inputs.
   * @return the list of files and time values corresponding to the boundary condition.
   * @note If the containers cannot be filled correctly, the code will exit with failure termination code.
   */


  if (Bc_NF.size()<1){  // if file's name is empty
    cerr << " parameters.txt: ERROR: it is necessary to specify a file for this choice of initialization of boundary condition ."<<endl;
    exit(EXIT_FAILURE);
  }

  /* if file is exist */
  if (-1 == access ((path_input_directory+Bc_NF).c_str(), R_OK) ) {
    cerr << " parameters.txt: ERROR: the "<< Bc_NF<< " file does not exist in the directory Inputs." << endl ;
    exit(EXIT_FAILURE);
  }

  ifstream inputfile((path_input_directory+Bc_NF).c_str());
  if (!inputfile){
    cout << " parameters.txt : ERROR: cannot open "<< path_input_directory+Bc_NF << " file\n";
    exit(EXIT_FAILURE);
  }

  /* vbound contains the type of boundary condition.
   vbound[0] contains the index of x coordinate of the boundary condition: vchoice[0] = (int) vtab[0]/dx
   vbound[1] contains the type of boundary condition: vchoice[1] =  vtab[1].
   */
  map<int,int> vbound;

  /* vtab contains the list of files and time values corresponding to the boundary condition.
   vtab[0] contains time values.
   vtab[1] contains the name of the file containing the values of the boundary condition.
   */
  map<double,string> vtab;

  map<int,SCALAR>  vq; // container containing discharges [m3/s]
  map<int,SCALAR>  vh; // container containing water heights [m]

  int line_number = 0;
  string line;
  string file_name;
  SCALAR time_value;
  SCALAR time_value_prev=0.0;
  char car;// character used to identify the beginning of a comment

  // As long as the end of the file is not reached, the next line is read.
  while (!inputfile.eof()){
    line_number++;

    getline (inputfile, line, '\n');       // read a line
    istringstream entree (line);
    if (entree >> time_value >> file_name){

      //Error if the times are decreasing.
      if(time_value < time_value_prev){
        cerr<< Bc_NF <<": ERROR: line "<< line_number <<", t = "<< time_value <<" must be greater than " << time_value_prev <<"."<< endl;
        inputfile.close();
        exit(EXIT_FAILURE);
      }

      /* Opens the file corresponding to the time changing and
       verifies that it contains all information to describe the boundary condition */
      vbound  =  fill_array_bc_inhomogeneous(file_name, path_input_directory,BC, vq , vh);
      vq.clear();
      vh.clear();
      vtab[time_value] = file_name;
    }
    else{
      car='#';
      istringstream entree_car (line);
      entree_car >> car;
      if (car != '#'){
        cout << Bc_NF  << ": WARNING: line "<< line_number <<"; a commentary should begin with the # symbol " <<endl;
      }

    }//end else
  }//end while

  if(0 == vtab.size()) {
    cerr << Bc_NF <<": there is not valid data in the file." << endl;
    inputfile.close();
    exit(EXIT_FAILURE);
  }

  //The last time is the Maximum finite representable floating-point number,
  //so it's not necessary to write the final time in the file.
  vtab[DBL_MAX]=string("last_time");

  if ((vtab.begin())->first > 0.0){
    cerr<< Bc_NF <<": ERROR: the first time must be t = 0." << endl;
    inputfile.close();
    exit(EXIT_FAILURE);
  }

  /* Retrieves the values for time equal to 0 */
  vchoice =  fill_array_bc_inhomogeneous((vtab.begin())->second, path_input_directory,BC, imp_q , imp_h);
  inputfile.close();

  return vtab;
}

string Parameters::get_path_input_directory(void) const{

  /**
   * @details
   * @return The Name of the input directory.
   */

  return path_input_directory;
}



Parameters::Parameters(){
}

Parameters::~Parameters(){
}
