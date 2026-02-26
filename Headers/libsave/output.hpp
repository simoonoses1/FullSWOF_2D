/**
 * @file output.hpp
 * @author Christian Laguerre <christian.laguerre@math.cnrs.fr> (2012-2015)
 * @author Carine Lucas <carine.lucas@univ-orleans.fr> (2023)
 * @version 1.10.00
 * @date 2023-03-17
 *
 * @brief %Output format
 * @details 
 * Common part for all the output formats.
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

#ifndef PARAMETERS_HPP
#include "parameters.hpp"
#endif

#ifndef OUTPUT_HPP
#define OUTPUT_HPP

/** @class Output
 * @brief %Output format
 * @details 
 * Class that contains all the common declarations for the output formats.
 */


class Output{
  
  public :
    
    /** @brief Constructor */
    Output(Parameters &);
    
    /** @brief Function to be specified in each output format */
    virtual void write(const TAB &, const TAB &, const TAB &, const TAB &, const SCALAR &) = 0;
    
    /** @brief Saves the infiltrated and rain volumes */
    void check_vol(const SCALAR &, const SCALAR &, const SCALAR &, const SCALAR & , const SCALAR & , const SCALAR &) const;
    
    /** @brief Saves global values */
    void result(const SCALAR &, const clock_t &, const SCALAR & , const SCALAR & , const SCALAR & , const SCALAR & , const int &, const SCALAR &) const;
    
    /** @brief Saves the initial time */
    void initial(const TAB &, const TAB &, const TAB &, const TAB &) const;
    
    /** @brief Saves the final time */
    void final(const TAB &, const TAB &, const TAB &, const TAB &) const;
    
    /** @brief Saves the cumulated fluxes on the boundaries */
    SCALAR boundaries_flux(const SCALAR &, const TAB &, const TAB &, const SCALAR &, const SCALAR &, const int &, const int &);
  
    /** @brief Saves the fluxes on the left and right boundaries */
    void boundaries_flux_LR(const SCALAR &, const TAB &) const;
  
    /** @brief Saves the fluxes on the bottom and top boundaries */
    void boundaries_flux_BT(const SCALAR &, const TAB &) const;
    
    /** @brief Rounded value up to 10^-5*/
    SCALAR round5(const SCALAR &) const;
    
    /** @brief Destructor */
    virtual ~Output();
    
  protected :
    /** Number of cells in space in the first (x) direction. */
    const int NXCELL;
    /** Number of cells in space in the second (y) direction. */
    const int NYCELL;
    /** Space step in the first (x) direction. */
    const SCALAR DX;
    /** Space step in the second (y) direction. */
    const SCALAR DY;
    /** Name of the output directory.*/
    string outputDirectory;
    /** Name of the file where the verification of volumes is saved.*/
    string namefile_check_volume;
    /** Name of the file where the global results are saved.*/
    string namefile_res;
    /** Name of the file where the initialization is saved.*/
    string namefile_init;
    /** Name of the file where the final time is saved.*/
    string namefile_final;
    /** Name of the file where the cumulated boundary fluxes are saved.*/
    string namefile_Bound_flux;
    /** Name of the file where the bottom and top boundary fluxes are saved.*/
    string namefile_Bound_flux_BT;
    /** Name of the file where the left and right boundary fluxes are saved.*/
    string namefile_Bound_flux_LR;
    
  private:
    SCALAR fluxNycell_cum_T,fluxNxcell_cum_T, fluxy0_cum_T, fluxx0_cum_T; // values of cumulated flux at the boundaries
};
#endif

