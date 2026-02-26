/**
 * @file choice_save_specific_points.cpp
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2017)
 * @version 1.07.01
 * @date 2017-02-03
 *
 * @brief Choice of the output of the specific points
 * @details 
 * From the value of the corresponding parameter,
 * calls the savings at the specific points.
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

#include "choice_save_specific_points.hpp"

Choice_save_specific_points::Choice_save_specific_points(Parameters & par){
  
  /**
   * @details
   * Defines the number (0,1,>1) of points to save from the value given in the parameters file.
   * @param[in] par parameter, contains all the values from the parameters file.
   */
  
  switch (par.get_choice_specific_point()){
    case 0:
      out = new No_save(par);
      break;
    
    case 1:
      out = new One_point(par);
      break;
    
    case 2:
      out = new Several_points(par);
      break;
  }
}

void Choice_save_specific_points::save(const TAB & h, const TAB & u, const TAB & v, const SCALAR time){
  
  /**
   * @details
   * Calls the saving of the current time.
   * @param[in] h water height.
   * @param[in] u first component of the velocity.
   * @param[in] v second component of the velocity.
   * @param[in] time value of the current time.
   */
  
  out->save(h, u, v, time);
}

Choice_save_specific_points::~Choice_save_specific_points(){
  if (out != NULL){
    delete out;
    out = NULL;
  }
}
