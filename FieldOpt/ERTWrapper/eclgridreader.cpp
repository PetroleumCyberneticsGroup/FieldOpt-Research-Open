/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

 Modified by Alin G. Chitu (2016-2017)
<alin.chitu@tno.nl, chitu_alin@yahoo.com>

Modified 2017-2020 Mathias Bellout
<chakibbb-pcg@gmail.com>

This file is part of the FieldOpt project.

FieldOpt is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version
3 of the License, or (at your option) any later version.

FieldOpt is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
the GNU General Public License for more details.

You should have received a copy of the
GNU General Public License along with FieldOpt.
If not, see <http://www.gnu.org/licenses/>.
***********************************************************/

#include "eclgridreader.h"
#include <iostream>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/lexical_cast.hpp>

#include "ertwrapper_exceptions.h"

namespace ERTWrapper {
namespace ECLGrid {
Eigen::Vector3d ECLGridReader::GetCellCenter(int global_index) {
  double cx, cy, cz;
  ecl_grid_get_xyz1(ecl_grid_, global_index, &cx, &cy, &cz);
  return Eigen::Vector3d(cx, cy, cz);
}

std::vector<Eigen::Vector3d> ECLGridReader::GetCellCorners(int global_index) {
  std::vector<Eigen::Vector3d> corners;
  for (int i = 0; i < 8; ++i) {
    double x, y, z;
    ecl_grid_get_cell_corner_xyz1(ecl_grid_, global_index, i, &x, &y, &z);
    corners.push_back(Eigen::Vector3d(x, y, z));
  }
  return corners;
}

double ECLGridReader::GetCellVolume(int global_index) {
  return ecl_grid_get_cell_volume1(ecl_grid_, global_index);
}

ECLGridReader::ECLGridReader() {
  ecl_grid_ = 0;
  ecl_file_init_ = 0;
  poro_kw_ = 0;
  permx_kw_ = 0;
  permy_kw_ = 0;
  permz_kw_ = 0;
}

ECLGridReader::~ECLGridReader() {
  if (ecl_grid_ != 0)
    ecl_grid_free(ecl_grid_);
  if (ecl_file_init_ != 0)
    ecl_file_close(ecl_file_init_);
}

void ECLGridReader::ReadEclGrid(std::string file_name) {
  file_name_ = file_name;
  init_file_name_ = file_name;
  if (boost::algorithm::ends_with(file_name, ".EGRID"))
    boost::replace_all(init_file_name_, ".EGRID", ".INIT");

  else if (boost::algorithm::ends_with(file_name, ".GRID"))
    boost::replace_all(init_file_name_, ".GRID", ".INIT");

  if (ecl_grid_ == 0) {
    ecl_grid_ = ecl_grid_alloc(file_name_.c_str());
  } else {
    ecl_grid_free(ecl_grid_);
    ecl_grid_ = ecl_grid_alloc(file_name_.c_str());
  }

  if (ecl_file_init_ == 0) {
    ecl_file_init_ = ecl_file_open(init_file_name_.c_str(), 0);
    poro_kw_ = ecl_file_iget_named_kw(ecl_file_init_, "PORO", 0);
    permx_kw_ = ecl_file_iget_named_kw(ecl_file_init_, "PERMX", 0);
    permy_kw_ = ecl_file_iget_named_kw(ecl_file_init_, "PERMY", 0);
    permz_kw_ = ecl_file_iget_named_kw(ecl_file_init_, "PERMZ", 0);
  } else {
    ecl_file_close(ecl_file_init_);
    ecl_file_init_ = ecl_file_open(init_file_name_.c_str(), 0);
    poro_kw_ = ecl_file_iget_named_kw(ecl_file_init_, "PORO", 0);
    permx_kw_ = ecl_file_iget_named_kw(ecl_file_init_, "PERMX", 0);
    permy_kw_ = ecl_file_iget_named_kw(ecl_file_init_, "PERMY", 0);
    permz_kw_ = ecl_file_iget_named_kw(ecl_file_init_, "PERMZ", 0);
  }
}

int ECLGridReader::ConvertIJKToGlobalIndex(ECLGridReader::IJKIndex ijk) {
  if (ecl_grid_ == 0) throw GridNotReadException("Grid must be read before converting indices.");
  return ecl_grid_get_global_index3(ecl_grid_, ijk.i, ijk.j, ijk.k);
}

int ECLGridReader::ConvertIJKToGlobalIndex(int i, int j, int k) {
  if (ecl_grid_ == 0) throw GridNotReadException("Grid must be read before converting indices.");
  return ecl_grid_get_global_index3(ecl_grid_, i, j, k);
}

ECLGridReader::IJKIndex ECLGridReader::ConvertGlobalIndexToIJK(int global_index) {
  if (ecl_grid_ == 0) throw GridNotReadException("Grid must be read before converting indices.");
  int i, j, k;
  ecl_grid_get_ijk1(ecl_grid_, global_index, &i, &j, &k);
  ECLGridReader::IJKIndex ijk;
  ijk.i = i; ijk.j = j; ijk.k = k;
  return ijk;
}

int ECLGridReader::ConvertMatrixActiveIndexToGlobalIndex(int index) {
  if (ecl_grid_ == 0) throw GridNotReadException("Grid must be read before converting indices.");
  else return ecl_grid_get_global_index1A(ecl_grid_, index);
}

ECLGridReader::Dims ECLGridReader::Dimensions() {
  ECLGridReader::Dims dims;
  if (ecl_grid_ == 0) throw GridNotReadException("Grid must be read before getting dimensions.");
  int x, y, z;
  ecl_grid_get_dims(ecl_grid_, &x, &y, &z, NULL);
  dims.nx = x; dims.ny = y; dims.nz = z;
  return dims;
}

int ECLGridReader::NumActiveMatrixCells() {
  if (ecl_grid_ == 0) throw GridNotReadException("Grid must be read before getting the number of active cells.");
  else return ecl_grid_get_nactive(ecl_grid_);
}

int ECLGridReader::NumActiveFractureCells() {
  if (ecl_grid_ == 0) throw GridNotReadException("Grid must be read before getting the number of active cells.");
  else return ecl_grid_get_nactive_fracture(ecl_grid_);
}

bool ECLGridReader::IsCellActive(int global_index) {
  if (ecl_grid_ == 0) throw GridNotReadException("Grid must be read before getting the active status of cells.");
  else return (IsCellMatrixActive(global_index) || IsCellFractureActive(global_index));
}

bool ECLGridReader::IsCellMatrixActive(int global_index) {
  if (ecl_grid_ == 0) throw GridNotReadException("Grid must be read before getting the active status of cells.");
  else
  {
    int active_matrix_index = ecl_grid_get_active_index1(ecl_grid_, global_index);
    if (active_matrix_index < 0) return false;
    else return true;
  }
}

bool ECLGridReader::IsCellFractureActive(int global_index) {
  if (ecl_grid_ == 0) throw GridNotReadException("Grid must be read before getting the active status of cells.");
  else
  {
    int active_fracture_index = ecl_grid_get_active_fracture_index1(ecl_grid_, global_index);
    if (active_fracture_index < 0) return false;
    else return true;
  }
}

std::vector<double> ECLGridReader::GetCellDxDyDz(int global_index) {
  std::vector<double> dxdydz(3);
  dxdydz[0] = ecl_grid_get_cell_dx1(ecl_grid_, global_index);
  dxdydz[1] = ecl_grid_get_cell_dy1(ecl_grid_, global_index);
  dxdydz[2] = ecl_grid_get_cell_dz1(ecl_grid_, global_index);
  return dxdydz;
}

ECLGridReader::Cell ECLGridReader::GetGridCell(int global_index) {
  if (!GlobalIndexIsInsideGrid(global_index))
    throw InvalidIndexException("The global index "
                                  + boost::lexical_cast<std::string>(global_index)
                                  + " is outside the grid.");
  if (ecl_grid_ == 0) throw GridNotReadException("Grid must be read before getting grid cells.");
  ECLGridReader::Cell cell;
  cell.global_index = global_index;
  cell.volume = GetCellVolume(global_index);
  cell.corners = GetCellCorners(global_index);
  cell.center = GetCellCenter(global_index);
  cell.matrix_active = IsCellMatrixActive(global_index);
  cell.fracture_active = IsCellFractureActive(global_index);
  auto dxdydz = GetCellDxDyDz(global_index);
  cell.dx = dxdydz[0];
  cell.dy = dxdydz[1];
  cell.dz = dxdydz[2];

  // Get properties from the INIT file - only possible if the cell is active
  // Matrix grid
  int active_index = ecl_grid_get_active_index1(ecl_grid_, global_index);
  if (active_index >= 0)
  {
    cell.porosity.push_back(ecl_kw_iget_as_double(poro_kw_, active_index));
    cell.permx.push_back(ecl_kw_iget_as_double(permx_kw_, active_index));
    cell.permy.push_back(ecl_kw_iget_as_double(permy_kw_, active_index));
    cell.permz.push_back(ecl_kw_iget_as_double(permz_kw_, active_index));
  }

  // Fracture grid
  active_index = NumActiveMatrixCells() + ecl_grid_get_active_fracture_index1(ecl_grid_, global_index);
  if (active_index >= NumActiveMatrixCells())
  {
    cell.porosity.push_back(ecl_kw_iget_as_double(poro_kw_, active_index));
    cell.permx.push_back(ecl_kw_iget_as_double(permx_kw_, active_index));
    cell.permy.push_back(ecl_kw_iget_as_double(permy_kw_, active_index));
    cell.permz.push_back(ecl_kw_iget_as_double(permz_kw_, active_index));
  }

  return cell;
}

int ECLGridReader::GlobalIndexOfCellEnvelopingPoint(double x, double y, double z, int initial_guess) {
  if (ecl_grid_ == 0) throw GridNotReadException("Grid must be read before searching for cells.");
  return ecl_grid_get_global_index_from_xyz(ecl_grid_, x, y, z, initial_guess);
}

bool ECLGridReader::GlobalIndexIsInsideGrid(int global_index) {
  Dims dims = Dimensions();
  return global_index < dims.nx * dims.ny * dims.nz;
}

ECLGridReader::Cell ECLGridReader::FindSmallestCell() {
  auto dims = Dimensions();
  int max_index = dims.nx * dims.ny * dims.nz;
  int index_with_smallest_volume = 0;
  double smallest_volume = 1e7;
  double volume;
  for (int global_index = 0; global_index < max_index; ++global_index) {
    if (IsCellActive(global_index)) {
      volume = GetCellVolume(global_index);
      if (volume < smallest_volume) {
        index_with_smallest_volume = global_index;
        smallest_volume = volume;
      }
    }
  }
  return GetGridCell(index_with_smallest_volume);
}

}
}
