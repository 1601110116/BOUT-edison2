
/*******************************************************************************
 * Integral operators
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 *******************************************************************************/

#ifndef __INTEGROPS_H__
#define __INTEGROPS_H__

#include "field3d.hxx"
#include "field2d.hxx"

#include "bout_types.hxx"

// Parallel derivative (central differencing)
// const Field2D IGrad_par(const Field2D &var, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
// const Field2D IGrad_par(const Field2D &var, DIFF_METHOD method, CELL_LOC outloc=CELL_DEFAULT);

// const Field3D IGrad_par(const Field3D &var, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
// const Field3D IGrad_par(const Field3D &var, DIFF_METHOD method, CELL_LOC outloc=CELL_DEFAULT);


const Field3D iSign_kpar(const Field3D &var, BoutReal k0=1.0, const int N=8);
const Field3D iSign_kpar_wcoll(const Field3D &vin, const Field3D &k0, const int N=0);
const Field2D iSign_kpar_wcoll(const Field2D &vin, const Field2D &k0, const int N=0);

#endif /* __INTEGROPS_H__ */
