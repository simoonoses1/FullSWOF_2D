/**
 * @file scheme.cpp
 * @author Olivier Delestre <olivierdelestre41@yahoo.fr> (2008)
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @author Carine Lucas <carine.lucas@univ-orleans.fr> (2020-2023)
 * @version 1.10.00
 * @date 2023-03-17
 *
 * @brief Numerical scheme
 * @details
 * Common part for all the numerical schemes.
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

#include "scheme.hpp"


Scheme::Scheme(Parameters & par):NXCELL(par.get_Nxcell()),NYCELL(par.get_Nycell()),ORDER(par.get_order()),T(par.get_T()),NBTIMES(par.get_nbtimes()),SCHEME_TYPE(par.get_scheme_type()),DX(par.get_dx()),DY(par.get_dy()),CFL_FIX(par.get_cflfix()),DT_FIX(par.get_dtfix()),FRICCOEF(par.get_friccoef()),L_IMP_Q(par.get_left_imp_discharge()),L_IMP_H(par.get_left_imp_h()),R_IMP_Q(par.get_right_imp_discharge()),R_IMP_H(par.get_right_imp_h()),B_IMP_Q(par.get_bottom_imp_discharge()),B_IMP_H(par.get_bottom_imp_h()),T_IMP_Q(par.get_top_imp_discharge()),T_IMP_H(par.get_top_imp_h()),fs2d_par(par){

  /**
   * @details
   * Initializations and allocations.
   * @param[in] par parameter, contains all the values from the parameters file.
   */

  allocation();

  cur_time=0;

  n=0;//initialization of the variable for the time loop

  /*----------------------------------------------------------------------- */
  Prain = new Choice_rain(par);
  //Initialization Rain
  Prain->rain_func(cur_time,Rain);
  Volrain_Tot=0.;

  // Initialization of the topography
  topo = new Choice_init_topo(par);
  topo->initialization(z);

  // Initialization of h, u, v
  huv_init = new Choice_init_huv(par);
  huv_init->initialization(h,u,v);
  Vol_of_tot=0.;

  flux_num = new Choice_flux(par.get_flux());

  for (int i=1 ; i<=NXCELL ; i++){
    for (int j=1 ; j<=NYCELL ; j++){
      q1[i][j] = u[i][j]*h[i][j];
      q2[i][j] = v[i][j]*h[i][j];
    } //end for j
  } //end for i

  Total_volume_outflow = 0.;

  fric = new Choice_friction(par);

  I = new Choice_infiltration(par);

  for (int i=1 ; i<=NXCELL ; i++){
    for (int j=1 ; j<=NYCELL ; j++){
      Vin_tot[i][j] = 0.;
    } //end for j
  } //end for i
  Vol_inf_tot_cumul=0.;


  // Left Boundary condition
  nchoice_Lbound = par.get_Lbound();
  Lbound = new Choice_condition(nchoice_Lbound,par,z,-1,0);
  left_times_files = par.get_times_files_Lbound();
  p_left_times_files = left_times_files.begin();
  Lbound_type = par.get_type_Lbound();
  is_Lbound_changed = false;
  // Right Boundary condition
  nchoice_Rbound = par.get_Rbound();
  Rbound = new Choice_condition(nchoice_Rbound,par,z,1,0);
  right_times_files = par.get_times_files_Rbound();
  p_right_times_files = right_times_files.begin();
  Rbound_type = par.get_type_Rbound();
  is_Rbound_changed = false;
  // Bottom Boundary condition
  nchoice_Bbound = par.get_Bbound();
  Bbound = new Choice_condition(nchoice_Bbound,par,z,0,-1);
  bottom_times_files = par.get_times_files_Bbound();
  p_bottom_times_files = bottom_times_files.begin();
  Bbound_type = par.get_type_Bbound();
  is_Bbound_changed = false;
  // Top Boundary condition
  nchoice_Tbound = par.get_Tbound();
  Tbound = new Choice_condition(nchoice_Tbound,par,z,0,1);
  top_times_files = par.get_times_files_Tbound();
  p_top_times_files = top_times_files.begin();
  Tbound_type = par.get_type_Tbound();
  is_Tbound_changed = false;


  dt_max=min(DX*CFL_FIX,DY*CFL_FIX);
  dt1=0.;

  // initialization of time value for the output
  if (0 == NBTIMES){     //in this case we don't call any function to store the values of variables (h,u ..)
    dt_output=2*T;       //computation of the step of time to save the variables in huz_evolution.dat file
    T_output=dt_output;  //Initialization of the variable used to save the variables in huz_evolution.dat file
                         // T_output=2*T but it can take any value greater than T, because we don't want to call
                         // the method to save picture in the final time (T).
  }
  else{
    dt_output=T/(NBTIMES-1); //computation of the step of time to save the variables in huz_evolution.dat file
    T_output=dt_output;      //Initialization of the variable used to save the variables in huz_evolution.dat file
  }

  out = new Choice_output(par);

  //storage of the topography
  out->initial(z, h, u,v);
  //storage the initialization of the main variables
  out->write(h,u,v,z,cur_time);

  string suffix_outputs = par.get_suffix();

  verif=1;

  out_specific_points = new Choice_save_specific_points(par);

}

void Scheme:: maincalcflux(SCALAR cflfix, SCALAR T , SCALAR curtime , SCALAR dt_max , SCALAR dt , SCALAR & dt_cal){

  /**
   * @details
   * First part.
   * Construction of variables for hydrostatic reconstruction.
   * Fluxes in the two directions.
   * Computation of the time step for the fixed cfl.
   * This calculation is called once at the order 1, and twice at the second order.
   * @param[in] cflfix fixed cfl.
   * @param[in] T final time (unused).
   * @param[in] curtime current time.
   * @param[in] dt_max maximum value of the time step.
   * @param[in] dt time step.
   * @param[out] dt_cal effective time step.
   * @warning the CFL condition is not satisfied: CFL > ***
   */

  (void) T; //unused variable
  (void) curtime; //unused variable

  SCALAR dt_tmp,dtx,dty;
  SCALAR velocity_max_x,velocity_max_y;  //temporary velocity to verify if clf > cflfix
  dtx=dty=dt_max;
  velocity_max_x=velocity_max_y=-VE_CA;
  for (int i=1 ; i<=NXCELL+1 ; i++){
    for (int j=1 ; j<NYCELL+1 ; j++){
      rec_hydro.calcul(h1r[i-1][j],h1l[i][j],delz1[i-1][j]);
      h1right[i-1][j] = rec_hydro.get_hhydro_l();
      h1left[i][j] = rec_hydro.get_hhydro_r();
      flux_num->calcul(h1right[i-1][j],u1r[i-1][j],v1r[i-1][j],h1left[i][j],u1l[i][j],v1l[i][j]);
      f1[i][j] = flux_num->get_f1();
      f2[i][j] = flux_num->get_f2();
      f3[i][j] = flux_num->get_f3();

      if (fabs(flux_num->get_cfl()*dt/DX) < EPSILON){
        dt_tmp=dt_max;
      }
      else{
        //  dt_tmp=min(T-curtime,cflfix*DX/flux_num->get_cfl());
        dt_tmp=cflfix*DX/flux_num->get_cfl();
      }
      dtx=min(min(dt,dt_tmp),dtx);
      velocity_max_x=max(velocity_max_x,flux_num->get_cfl());
    } //end for j
  } //end for i

  for (int i=1 ; i<NXCELL+1 ; i++){
    for (int j=1 ; j<=NYCELL+1 ; j++){
      rec_hydro.calcul(h2r[i][j-1],h2l[i][j],delz2[i][j-1]);
      h2right[i][j-1] = rec_hydro.get_hhydro_l();
      h2left[i][j] = rec_hydro.get_hhydro_r();
      flux_num->calcul(h2right[i][j-1],v2r[i][j-1],u2r[i][j-1],h2left[i][j],v2l[i][j],u2l[i][j]);
      g1[i][j] = flux_num->get_f1();
      g2[i][j] = flux_num->get_f3();
      g3[i][j] = flux_num->get_f2();

      if (fabs(flux_num->get_cfl()*dt/DY) < EPSILON){
        dt_tmp=dt_max;
      }
      else{
        dt_tmp=cflfix*DY/flux_num->get_cfl();
      }
      dty=min(min(dt,dt_tmp),dty);
      velocity_max_y=max(velocity_max_y,flux_num->get_cfl());
    } //end for j
  } //end for i

  if (1 == SCHEME_TYPE){
    dt_cal=min(dtx,dty);
  }
  else{
    if ((velocity_max_x*DT_FIX/DX>cflfix)||(velocity_max_y*DT_FIX/DY>cflfix)){
      cout << " the CFL condition is not satisfied: CFL >"<<cflfix << endl;
      exit(1);
    } //end if
    dt_cal=DT_FIX;
  }
}

void Scheme:: maincalcscheme(TAB & he, TAB & ve1, TAB & ve2, TAB & qe1, TAB & qe2, TAB & hes, TAB & ves1, TAB & ves2, TAB & qes1, TAB & qes2,TAB & Vin, SCALAR curtime, SCALAR dt, int n){

  /**
   * @details
   * Second part.
   * Computation of h, u and v.
   * This calculation is called once at the order 1, and twice at the second order.
   * @param[in] he water height.
   * @param[in] ve1 first component of the velocity.
   * @param[in] ve2 second component of the velocity.
   * @param[in] qe1 first component of the discharge (unused).
   * @param[in] qe2 second component of the discharge (unused).
   * @param[out] hes water height.
   * @param[out] ves1 first component of the velocity.
   * @param[out] ves2 second component of the velocity.
   * @param[out] qes1 first component of the discharge.
   * @param[out] qes2 second component of the discharge.
   * @param[out] Vin infiltrated volume
   * @param[in] curtime current time.
   * @param[in] dt time step.
   * @param[in] n number of iterations (unused).
   * @note In DEBUG mode, the programme will save three other files with boundaries fluxes.
   */

  (void) qe1; //unused variable
  (void) qe2; //unused variable
  (void) n; //unused variable

  /*-------------- Periodic boundary conditions ----------------------------------------------------------*/
  //In case of periodic boundary conditions,  the flux at the boundaries
  // (Left and Right, Bottom and Top) must be the same.
  //Moreover we need to consider the direction of the discharge in order to exchange the flux.

  for (int i=1 ; i<NXCELL+1 ; i++){
    if ((4==nchoice_Tbound[i]) && (4==nchoice_Bbound[i])){
      if ((ve2[i][1] > 0.) && (ve2[i][NYCELL]> 0.)){ //the direction of flow is Bottom to the Top
        g1[i][1]= g1[i][NYCELL+1];
      }
      if ((ve2[i][1] < 0.) && (ve2[i][NYCELL]< 0.)){ //the direction of flow is Top to the Bottom
        g1[i][NYCELL+1]  = g1[i][1];
      }
    }
  }

  for (int j=1 ; j<NYCELL+1 ; j++){
    if ((4==nchoice_Rbound[j]) && (4==nchoice_Lbound[j])){
      if ((ve1[1][j] > 0.) && (ve1[NXCELL][j]> 0.)){ //the direction of flow is Left to the Right
        f1[1][j]= f1[NXCELL+1][j];
      }
      //    for (int j=1 ; j<NYCELL+1 ; j++){
      if ((ve1[1][j] < 0.) && (ve1[NXCELL][j]< 0.)){ //the direction of flow is Right to the Left
	f1[NXCELL+1][j] = f1[1][j];
      }
    }
  }

  /*-------------- Rainfall and infiltration --------------------------------------------------------------*/

  Prain->rain_func(curtime,Rain);

  tx=dt/DX;
  ty=dt/DY;

  /*-------------- Main computation ------------------------------------------------------------------------*/

  if (1 == verif){
    dt_first=dt;
  }

  for (int i=1 ; i<NXCELL+1 ; i++){
    for (int j=1 ; j<NYCELL+1 ; j++){
      // Solution of the equation of mass conservation (First equation of Shallow-Water)
      hes[i][j] = he[i][j]-tx*(f1[i+1][j]-f1[i][j])-ty*(g1[i][j+1]-g1[i][j])+Rain[i][j]*dt;

    } //end for j
  } //end for i


  //Infiltration
  I->calcul(hes,Vin,dt);
  hes = I->get_hmod();
  Vin = I->get_Vin();


  for (int i=1 ; i<NXCELL+1 ; i++){
    for (int j=1 ; j<NYCELL+1 ; j++){

      //Solution of the equation of momentum (Second and third equation of Shallow-Water)
      // This expression for the flux (instead of the differences of the squares) avoids numerical errors
      // see http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html section "Cancellation".

      qes1[i][j] = (he[i][j]*ve1[i][j]-tx*(f2[i+1][j]-f2[i][j]+GRAV_DEM*((h1left[i][j]-h1l[i][j])*(h1left[i][j]+h1l[i][j])+(h1r[i][j]-h1right[i][j])*(h1r[i][j]+h1right[i][j])+(h1l[i][j]+h1r[i][j])*delzc1[i][j]))-ty*(g2[i][j+1]-g2[i][j]));
      
      qes2[i][j] = (he[i][j]*ve2[i][j]-tx*(f3[i+1][j]-f3[i][j])-ty*(g3[i][j+1]-g3[i][j]+GRAV_DEM*((h2left[i][j]-h2l[i][j])*(h2left[i][j]+h2l[i][j])+(h2r[i][j]-h2right[i][j])*(h2r[i][j]+h2right[i][j])+(h2l[i][j]+h2r[i][j])*delzc2[i][j])));

    } //end for j
  } //end for i

  //Friction
  fric->calcul(ve1,ve2,hes,qes1,qes2,dt);
  qes1 = fric->get_q1mod();
  qes2 = fric->get_q2mod();

  for (int i=1 ; i<NXCELL+1 ; i++){
    for (int j=1 ; j<NYCELL+1 ; j++){
      if (hes[i][j] > HE_CA){
        ves1[i][j] = qes1[i][j]/hes[i][j];
        ves2[i][j] = qes2[i][j]/hes[i][j];
      }
      else{ // Case of height of water is zero.
        ves1[i][j] = 0.;
        ves2[i][j] = 0.;
      }
    } //end for j
  } //end for i

  // The total cumulated rain's computed
  for (int i=1 ; i<NXCELL+1 ; i++){
    for (int j=1 ; j<NYCELL+1 ; j++){
      Volrain_Tot += Rain[i][j]*(dt-dt_first*(1-verif))*(1./ORDER)*DX*DY;
    }
  }

  Total_volume_outflow =  out->boundaries_flux(curtime,f1,g1, dt, dt_first,ORDER,verif);

#ifdef DEBUG
  out->boundaries_flux_LR(curtime,f1);
  out->boundaries_flux_BT(curtime,g1);
#endif

}


void Scheme::boundary(TAB & h_tmp,TAB & u_tmp ,TAB & v_tmp,SCALAR time_tmp,const int NODEX, const int NODEY){

  /**
   * @details
   * @param[in, out] h_tmp water height.
   * @param[in, out] u_tmp first component of the velocity.
   * @param[in, out] v_tmp second component of the velocity.
   * @param[in] time_tmp current time.
   * @param[in] NODEX number of space cells in the first (x) direction.
   * @param[in] NODEY number of space cells in the second (y) direction.
   */


  /* If we are in the case where the users have chosen to use a file for boundary condition  (1==Lbound_type):*/
  if (1==Lbound_type){
    if (time_tmp >= p_left_times_files->first){
    /* We extract from file corresponding to the current time the choice of boundary conditions, the discharges [m3/s]
       and water heights [m]. */
      L_choice_bound = fs2d_par.fill_array_bc_inhomogeneous(p_left_times_files->second, fs2d_par.get_path_input_directory(),'L', L_IMP_Q , L_IMP_H);

      for (int j=1 ; j<=NODEY ; j++){
	/*
	  The value of the imposed discharge per cell in the boundary condition is
	  in m2/s
	*/
	L_IMP_Q[j]/= DY;
      }

      /* We update the choices of boundary conditions */
      Lbound->setChoice(L_choice_bound);
      p_left_times_files++;
      is_Lbound_changed = true;
    }
  }

  /* If we are in the case where the users have chosen to use a file for boundary condition  (1==Rbound_type):*/
  if (1==Rbound_type){
    if (time_tmp >= p_right_times_files->first){
     /* We extract from file corresponding to the current time the choice of boundary conditions, the discharges [m3/s]
       and water heights [m]. */
     R_choice_bound  =  fs2d_par.fill_array_bc_inhomogeneous(p_right_times_files->second, fs2d_par.get_path_input_directory(),'R', R_IMP_Q , R_IMP_H);

      for (int j=1 ; j<=NODEY ; j++){
	/*
	  The value of the imposed discharge per cell in the boundary condition is
	  in m2/s
	*/
	R_IMP_Q[j]/= DY;

      }

      /* We update the choices of boundary conditions */
      Rbound->setChoice(R_choice_bound);
      p_right_times_files++;
      is_Rbound_changed = true;
    }
  }

  /* if the boundary conditions has been updated then we verify the periodic condition is compatible on each side */
  if (is_Lbound_changed || is_Rbound_changed  ){
    for (int j=1 ; j<=NODEY ; j++){
      if (((R_choice_bound[j] == 4) && (L_choice_bound[j] != 4)) || ((R_choice_bound[j] != 4) && (L_choice_bound[j] == 4)) ){
	cerr << p_left_times_files->second << " : ERROR: you must choose a periodic bottom condition."<< endl;
	exit(EXIT_FAILURE);
      }
    }
    is_Lbound_changed = false;
    is_Rbound_changed = false;
  }



  /* If we are in the case where the users have chosen to use a file for boundary condition  (1==Bbound_type):*/
  if (1==Bbound_type){
    if (time_tmp >= p_bottom_times_files->first){
     /* We extract from file corresponding to the current time the choice of boundary conditions, the discharges [m3/s]
       and water heights [m]. */
	B_choice_bound  = fs2d_par.fill_array_bc_inhomogeneous(p_bottom_times_files->second, fs2d_par.get_path_input_directory(),'B', B_IMP_Q , B_IMP_H);

	for (int i=1 ; i<=NODEX ; i++){
	  /*
	    The value of the imposed discharge per cell in the boundary condition is
	    in m2/s
	  */
	  B_IMP_Q[i]/= DX;

	}

      /* We update the choices of boundary conditions */
	Bbound->setChoice(B_choice_bound);
	p_bottom_times_files++;
	is_Bbound_changed = true;
      }
    }

  /* If we are in the case where the users have chosen to use a file for boundary condition  (1==Tbound_type):*/
  if (1==Tbound_type){
    if (time_tmp >= p_top_times_files->first){
      /* We extract from file corresponding to the current time the choice of boundary conditions, the discharges [m3/s]
	 and water heights [m]. */
      T_choice_bound = fs2d_par.fill_array_bc_inhomogeneous(p_top_times_files->second, fs2d_par.get_path_input_directory(),'T', T_IMP_Q , T_IMP_H);

      for (int i=1 ; i<=NODEX ; i++){
	/*
	  The value of the imposed discharge per cell in the boundary condition is
	  in m2/s
	*/
	T_IMP_Q[i]/= DX;

      }

      /* We update the choices of boundary conditions */
      Tbound->setChoice(T_choice_bound);
      p_top_times_files++;
      is_Tbound_changed = true;
    }
  }

  /* if the boundary conditions has been updated then we verify the periodic condition is compatible on each side */
  if (is_Bbound_changed || is_Tbound_changed){
    for (int i=1 ; i<=NODEX ; i++){
      if (((T_choice_bound[i] == 4) && (B_choice_bound[i] != 4)) || ((T_choice_bound[i] != 4) && (B_choice_bound[i] == 4)) ){
	cerr <<" parameters.txt: ERROR: you must choose a periodic bottom condition."<< endl;
	exit(EXIT_FAILURE);
      }
    }
  }




  for (int j=1 ; j<NODEY+1 ; j++){
    Lbound->setXY(0,j);
    Lbound->calcul(h_tmp[1][j],u_tmp[1][j],v_tmp[1][j],L_IMP_H[j],L_IMP_Q[j],h_tmp[NODEX][j],u_tmp[NODEX][j],v_tmp[NODEX][j], time_tmp,-1,0);
    h_tmp[0][j] = Lbound->get_hbound();
    u_tmp[0][j] = Lbound->get_unormbound();
    v_tmp[0][j] = Lbound->get_utanbound();

    Rbound->setXY(NODEX+1,j);
    Rbound->calcul(h_tmp[NODEX][j],u_tmp[NODEX][j],v_tmp[NODEX][j],R_IMP_H[j],R_IMP_Q[j],h_tmp[1][j],u_tmp[1][j],v_tmp[1][j], time_tmp,1,0);
    h_tmp[NODEX+1][j] = Rbound->get_hbound();
    u_tmp[NODEX+1][j] = Rbound->get_unormbound();
    v_tmp[NODEX+1][j] = Rbound->get_utanbound();
  } //end for j

  for (int i=1 ; i<NODEX+1 ; i++){
    Bbound->setXY(i,0);
    Bbound->calcul(h_tmp[i][1],v_tmp[i][1],u_tmp[i][1],B_IMP_H[i],B_IMP_Q[i],h_tmp[i][NODEY],v_tmp[i][NODEY],u_tmp[i][NODEY], time_tmp,0,-1);
    h_tmp[i][0] = Bbound->get_hbound();
    u_tmp[i][0] = Bbound->get_utanbound();
    v_tmp[i][0] = Bbound->get_unormbound();

    Tbound->setXY(i,NODEY+1);
    Tbound->calcul(h_tmp[i][NODEY],v_tmp[i][NODEY],u_tmp[i][NODEY],T_IMP_H[i],T_IMP_Q[i],h_tmp[i][1],v_tmp[i][1],u_tmp[i][1], time_tmp,0,1);
    h_tmp[i][NODEY+1] = Tbound->get_hbound();
    u_tmp[i][NODEY+1] = Tbound->get_utanbound();
    v_tmp[i][NODEY+1] = Tbound->get_unormbound();
  } //end for i

}


SCALAR Scheme::froude_number(TAB h_s,TAB u_s,TAB v_s){

  /**
   * @details
   * Mean value in space of the Froude number at the final time.
   * @param[in] h_s water height.
   * @param[in] u_s first component of the velocity.
   * @param[in] v_s second component of the velocity.
   * @return The mean Froude number \f$ \displaystyle \frac{\sqrt{u_s^2+v_s^2}}{\sqrt{gh_s}} \f$.
   */

  SCALAR Fr;
  int i,j;
  SCALAR stock_u,stock_v,stock_h;

  stock_u=0.;
  stock_v=0.;
  stock_h=0.;

  for(j=1;j<=NYCELL;j++){
    for(i=1;i<=NXCELL;i++){
      stock_u+=u_s[i][j];
      stock_v+=v_s[i][j];
      stock_h+=h_s[i][j];
    }
  }

  stock_u=stock_u/(NXCELL*NYCELL);
  stock_v=stock_v/(NXCELL*NYCELL);
  stock_h=stock_h/(NXCELL*NYCELL);

  Fr=sqrt((stock_u*stock_u+stock_v*stock_v)/(GRAV*stock_h));

  return Fr;
}


void Scheme::allocation(){

  /**
   * @details
   * Allocation of Scheme#z, Scheme#h, Scheme#u, Scheme#v, Scheme#q1, Scheme#q2, Scheme#Vin_tot,
   * Scheme#hs, Scheme#us, Scheme#vs, Scheme#qs1, Scheme#qs2, Scheme#f1, Scheme#f2, Scheme#f3, Scheme#g1, Scheme#g2, Scheme#g3
   * Scheme#h1left, Scheme#h1l, Scheme#u1l, Scheme#v1l, Scheme#h1right, Scheme#h1r, Scheme#u1r, Scheme#v1r,
   * Scheme#h2left, Scheme#h2l, Scheme#u2l, Scheme#v2l, Scheme#h2right, Scheme#h2r, Scheme#u2r, Scheme#v2r,
   * Scheme#delz1, Scheme#delz2, Scheme#delzc1, Scheme#delzc2, Scheme#Rain.
   */

  Rain.resize(NXCELL+2); // i : 0->NXCELL+1

  z.resize(NXCELL+2); // i : 0->NXCELL+1
  h.resize(NXCELL+2); // i : 0->NXCELL+1
  u.resize(NXCELL+2); // i : 0->NXCELL+1
  v.resize(NXCELL+2); // i : 0->NXCELL+1
  q1.resize(NXCELL+1); // i : 1->NXCELL
  q2.resize(NXCELL+1); // i : 1->NXCELL

  Vin_tot.resize(NXCELL+1); // i : 1->NXCELL
  hs.resize(NXCELL+2); // i : 0->NXCELL+1
  us.resize(NXCELL+2); // i : 0->NXCELL+1
  vs.resize(NXCELL+2); // i : 0->NXCELL+1
  qs1.resize(NXCELL+1); // i : 1->NXCELL
  qs2.resize(NXCELL+1); // i : 1->NXCELL
  f1.resize(NXCELL+2); // i : 1->NXCELL+1
  f2.resize(NXCELL+2); // i : 1->NXCELL+1
  f3.resize(NXCELL+2); // i : 1->NXCELL+1
  g1.resize(NXCELL+1); // i : 1->NXCELL
  g2.resize(NXCELL+1); // i : 1->NXCELL
  g3.resize(NXCELL+1); // i : 1->NXCELL
  h1left.resize(NXCELL+2); // i : 1->NXCELL+1
  h1l.resize(NXCELL+2); // i : 1->NXCELL+1
  u1l.resize(NXCELL+2); // i : 1->NXCELL+1
  v1l.resize(NXCELL+2); // i : 1->NXCELL+1
  h1right.resize(NXCELL+1); // i : 0->NXCELL
  h1r.resize(NXCELL+1); // i : 0->NXCELL
  u1r.resize(NXCELL+1); // i : 0->NXCELL
  v1r.resize(NXCELL+1); // i : 0->NXCELL
  h2left.resize(NXCELL+1); // i : 1->NXCELL
  h2l.resize(NXCELL+1); // i : 1->NXCELL
  u2l.resize(NXCELL+1); // i : 1->NXCELL
  v2l.resize(NXCELL+1); // i : 1->NXCELL
  h2right.resize(NXCELL+1); // i : 1->NXCELL
  h2r.resize(NXCELL+1); // i : 1->NXCELL
  u2r.resize(NXCELL+1); // i : 1->NXCELL
  v2r.resize(NXCELL+1); // i : 1->NXCELL
  delz1.resize(NXCELL+1); // i : 0->NXCELL
  delz2.resize(NXCELL+1); // i : 1->NXCELL
  delzc1.resize(NXCELL+1); // i : 1->NXCELL
  delzc2.resize(NXCELL+1); // i : 1->NXCELL

  Rain[0].resize(NYCELL+2); // j : 0->NYCELL+1
  z[0].resize(NYCELL+2); // j : 0->NYCELL+1
  h[0].resize(NYCELL+2); // j : 0->NYCELL+1
  u[0].resize(NYCELL+2); // j : 0->NYCELL+1
  v[0].resize(NYCELL+2); // j : 0->NYCELL+1
  hs[0].resize(NYCELL+2); // j : 0->NYCELL+1
  us[0].resize(NYCELL+2); // j : 0->NYCELL+1
  vs[0].resize(NYCELL+2); // j : 0->NYCELL+1
  h1right[0].resize(NYCELL+1); // j : 1->NYCELL
  h1r[0].resize(NYCELL+1); // j : 1->NYCELL
  u1r[0].resize(NYCELL+1); // j : 1->NYCELL
  v1r[0].resize(NYCELL+1); // j : 1->NYCELL
  delz1[0].resize(NYCELL+1); // j : 1->NYCELL

  for (int i=1 ; i<=NXCELL ; i++){
    Rain[i].resize(NYCELL+2); // j : 0->NYCELL+1
    z[i].resize(NYCELL+2); // j : 0->NYCELL+1
    h[i].resize(NYCELL+2); // j : 0->NYCELL+1
    u[i].resize(NYCELL+2); // j : 0->NYCELL+1
    v[i].resize(NYCELL+2); // j : 0->NYCELL+1

    q1[i].resize(NYCELL+1); // j : 1->NYCELL
    q2[i].resize(NYCELL+1); // j : 1->NYCELL

    Vin_tot[i].resize(NYCELL+1); // j : 1->NYCELL

    hs[i].resize(NYCELL+2); // j : 0->NYCELL+1
    us[i].resize(NYCELL+2); // j : 0->NYCELL+1
    vs[i].resize(NYCELL+2); // j : 0->NYCELL+1

    qs1[i].resize(NYCELL+1); // j : 1->NYCELL
    qs2[i].resize(NYCELL+1); // j : 1->NYCELL
    f1[i].resize(NYCELL+1); // j : 1->NYCELL
    f2[i].resize(NYCELL+1); // j : 1->NYCELL
    f3[i].resize(NYCELL+1); // j : 1->NYCELL

    g1[i].resize(NYCELL+2); // j : 1->NYCELL+1
    g2[i].resize(NYCELL+2); // j : 1->NYCELL+1
    g3[i].resize(NYCELL+2); // j : 1->NYCELL+1

    h1left[i].resize(NYCELL+1); // j : 1->NYCELL
    h1l[i].resize(NYCELL+1); // j : 1->NYCELL
    u1l[i].resize(NYCELL+1); // j : 1->NYCELL
    v1l[i].resize(NYCELL+1); // j : 1->NYCELL
    h1right[i].resize(NYCELL+1); // j : 1->NYCELL
    h1r[i].resize(NYCELL+1); // j : 1->NYCELL
    u1r[i].resize(NYCELL+1); // j : 1->NYCELL
    v1r[i].resize(NYCELL+1); // j : 1->NYCELL

    h2left[i].resize(NYCELL+2); // j : 1->NYCELL+1
    h2l[i].resize(NYCELL+2); // j : 1->NYCELL+1
    u2l[i].resize(NYCELL+2); // j : 1->NYCELL+1
    v2l[i].resize(NYCELL+2); // j : 1->NYCELL+1

    h2right[i].resize(NYCELL+1); // j : 0->NYCELL
    h2r[i].resize(NYCELL+1); // j : 0->NYCELL
    u2r[i].resize(NYCELL+1); // j : 0->NYCELL
    v2r[i].resize(NYCELL+1); // j : 0->NYCELL
    delz1[i].resize(NYCELL+1); // j : 1->NYCELL
    delz2[i].resize(NYCELL+1); // j : 0->NYCELL
    delzc1[i].resize(NYCELL+1); // j : 1->NYCELL
    delzc2[i].resize(NYCELL+1); // j : 1->NYCELL

  }

  Rain[NXCELL+1].resize(NYCELL+2); // j : 0->NYCELL+1
  z[NXCELL+1].resize(NYCELL+2); // j : 0->NYCELL+1
  h[NXCELL+1].resize(NYCELL+2); // j : 0->NYCELL+1
  u[NXCELL+1].resize(NYCELL+2); // j : 0->NYCELL+1
  v[NXCELL+1].resize(NYCELL+2); // j : 0->NYCELL+1
  hs[NXCELL+1].resize(NYCELL+2); // j : 0->NYCELL+1
  us[NXCELL+1].resize(NYCELL+2); // j : 0->NYCELL+1
  vs[NXCELL+1].resize(NYCELL+2); // j : 0->NYCELL+1

  f1[NXCELL+1].resize(NYCELL+1); // j : 1->NYCELL
  f2[NXCELL+1].resize(NYCELL+1); // j : 1->NYCELL
  f3[NXCELL+1].resize(NYCELL+1); // j : 1->NYCELL
  h1left[NXCELL+1].resize(NYCELL+1); // j : 1->NYCELL
  h1l[NXCELL+1].resize(NYCELL+1); // j : 1->NYCELL
  u1l[NXCELL+1].resize(NYCELL+1); // j : 1->NYCELL
  v1l[NXCELL+1].resize(NYCELL+1); // j : 1->NYCELL

}


void Scheme::deallocation(){

  /**
   * @details
   * Deallocation of Scheme#z, Scheme#h, Scheme#u, Scheme#v, Scheme#q1, Scheme#q2, Scheme#Vin_tot,
   * Scheme#hs, Scheme#us, Scheme#vs, Scheme#qs1, Scheme#qs2, Scheme#f1, Scheme#f2, Scheme#f3, Scheme#g1, Scheme#g2, Scheme#g3
   * Scheme#h1left, Scheme#h1l, Scheme#u1l, Scheme#v1l, Scheme#h1right, Scheme#h1r, Scheme#u1r, Scheme#v1r,
   * Scheme#h2left, Scheme#h2l, Scheme#u2l, Scheme#v2l, Scheme#h2right, Scheme#h2r, Scheme#u2r, Scheme#v2r,
   * Scheme#delz1, Scheme#delz2, Scheme#delzc1, Scheme#delzc2, Scheme#Rain.
   */

  delete Prain ;
  delete topo;
  delete huv_init;
  delete flux_num;
  delete fric;
  delete I;
  delete out;
  delete Lbound;
  delete Rbound;
  delete Bbound;
  delete Tbound;

  Rain[0].clear();
  z[0].clear();
  h[0].clear();
  u[0].clear();
  v[0].clear();
  hs[0].clear();
  us[0].clear();
  vs[0].clear();
  h1right[0].clear();
  h1r[0].clear();
  u1r[0].clear();
  v1r[0].clear();
  delz1[0].clear();

  for (int i=1 ; i<=NXCELL ; i++){

    Rain[i].clear();
    z[i].clear();
    Vin_tot[i].clear();
    h[i].clear();
    u[i].clear();
    v[i].clear();
    q1[i].clear();
    q2[i].clear();
    hs[i].clear();
    us[i].clear();
    vs[i].clear();
    qs1[i].clear();
    qs2[i].clear();
    f1[i].clear();
    f2[i].clear();
    f3[i].clear();
    g1[i].clear();
    g2[i].clear();
    g3[i].clear();
    h1left[i].clear();
    h1l[i].clear();
    u1l[i].clear();
    v1l[i].clear();
    h1right[i].clear();
    h1r[i].clear();
    u1r[i].clear();
    v1r[i].clear();
    h2left[i].clear();
    h2l[i].clear();
    u2l[i].clear();
    v2l[i].clear();
    h2right[i].clear();
    h2r[i].clear();
    u2r[i].clear();
    v2r[i].clear();
    delz1[i].clear();
    delz2[i].clear();
    delzc1[i].clear();
    delzc2[i].clear();
  }

  Rain[NXCELL+1].clear();
  z[NXCELL+1].clear();
  h[NXCELL+1].clear();
  u[NXCELL+1].clear();
  v[NXCELL+1].clear();
  hs[NXCELL+1].clear();
  us[NXCELL+1].clear();
  vs[NXCELL+1].clear();
  f1[NXCELL+1].clear();
  f2[NXCELL+1].clear();
  f3[NXCELL+1].clear();
  h1left[NXCELL+1].clear();
  h1l[NXCELL+1].clear();
  u1l[NXCELL+1].clear();
  v1l[NXCELL+1].clear();

  z.clear();
  h.clear();
  u.clear();
  v.clear();
  q1.clear();
  q2.clear();
  hs.clear();
  us.clear();
  vs.clear();
  h1right.clear();
  h1r.clear();
  h1left.clear();
  h1l.clear();
  u1l.clear();
  v1l.clear();
  u1r.clear();
  v1r.clear();
  h2left.clear();
  h2l.clear();
  u2l.clear();
  v2l.clear();
  h2right.clear();
  h2r.clear();
  u2r.clear();
  v2r.clear();
  f1.clear();
  f2.clear();
  f3.clear();
  g1.clear();
  g2.clear();
  g3.clear();
  delz1.clear();
  delz2.clear();
  delzc1.clear();
  delzc2.clear();
  qs1.clear();
  qs2.clear();
  Vin_tot.clear();

}

Scheme::~Scheme(){
  deallocation();
  delete out_specific_points;
#ifdef DEBUG
  cout<<"Deallocation of objects is finished"<< endl;
#endif
}
