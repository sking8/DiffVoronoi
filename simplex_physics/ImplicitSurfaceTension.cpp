//////////////////////////////////////////////////////////////////////////
// Implicitly solve surface tension
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifdef USE_CPX
#include "ImplicitSurfaceTension.h"
#include "ArrayFunc.h"

template<class T, int d>
void ImplicitSurfaceTension<T, d>::Solve(const T dt, FaceField<T, d>& velocity, const T sigma, const T narrow_band_width, const T dirac_band_width)
{
	static Timer timer;
	static FaceField<int, d> face_idx;
	face_idx.Resize(mac_grid.grid.cell_counts);

	for (int axis = 0;axis < d;axis++) {
		timer.Begin_Loop();

		//Grid<d> face_grid = mac_grid.face_grids[axis];
		Grid<d> face_grid = Grid<d>(mac_grid.grid.cell_counts + VectorDi::Unit(axis));
		
		//we call the index in serialized matrix "idx"
		//the coordinate of face node is "face_cell"
		//the serial id in face_grid of a face node is "c"
		const int cell_num = face_grid.Number_Of_Cells();
		//temporary structure for boosting, c_idx[c]=1 or corresponding idx
		static Array<int> c_idx;c_idx.resize(cell_num);
		face_grid.Exec_Each(
			[&](const VectorDi& face_cell) {
				int c = face_grid.Cell_Index(face_cell);
				if (Is_Levelset_Interface_Face(axis, face_cell, narrow_band_width)) c_idx[c] = 1;
				else c_idx[c] = -1;
			}
		);
		//timer.Record("check all availability");
		int N = 0;
		//[SERIAL] indexing every cell
		for (int i = 0;i < cell_num;i++) {
			if (c_idx[i] != -1) {
				c_idx[i] = N++;
			}
		}
		static Array<VectorDi> idx_coord; idx_coord.resize(N);
		face_grid.Exec_Each(
			[&](const VectorDi& face_cell) {
				int c = face_grid.Cell_Index(face_cell);
				int idx = c_idx[c];
				if (idx != -1) {
					face_idx(axis, face_cell) = idx;
					idx_coord[idx] = face_cell;
				}
				else face_idx(axis, face_cell) = -1;
			}
		);
		static Array<int> idx_row_nnz;idx_row_nnz.resize(N);
		//calculate non-zero numbers
		ArrayFunc::Calc_Each(
			[&](const int idx) {
				int row_nnz = 1;//diagonal element itself
				const VectorDi face_cell = idx_coord[idx];
				for (int i = 0;i < Grid<d>::Number_Of_Nb_C();i++) {
					VectorDi nb_face = Grid<d>::Nb_C(face_cell, i);
					if (face_grid.Valid_Cell(nb_face) && face_idx(axis, nb_face) != -1) {
						row_nnz++;
					}
				}
				return row_nnz;
			},
			idx_row_nnz
		);
		int non_zero_numbers = 0;
#pragma omp parallel for reduction(+:non_zero_numbers)
		for (int idx = 0;idx < N;idx++) {
			non_zero_numbers += idx_row_nnz[idx];
		}
		Eigen::SparseMatrix<T, Eigen::RowMajor, int> B;
		B.resize(N, N);
		B.resizeNonZeros(non_zero_numbers);

		//[SERIAL] iterate all elements, calculate a prefix sum
		int* outer = B.outerIndexPtr();outer[0] = 0;
		for (int idx = 0;idx < N;idx++) {
			outer[idx + 1] = outer[idx] + idx_row_nnz[idx];
		}

		//fill col ptr, val ptr, u_old
		static Array<T> u_old;u_old.resize(N);
		static Array<T> u_new;u_new.resize(N);
		ArrayFunc::Exec_Each(
			[&](const int idx) {
				VectorDi face_cell = idx_coord[idx];
				VectorD pos = mac_grid.Face_Center(axis, face_cell);
				auto kappa = levelset->Curvature(pos);
				VectorD normal = levelset->Normal(pos);
				u_old[idx] = velocity(axis, face_cell) - sigma * Levelset_Dirac(levelset->Phi(pos), dirac_band_width) * dt * kappa * normal[axis];
				T dia_coef = 1;
				Array<std::pair<int, T>> col_val_pairs;col_val_pairs.clear();
				for (int i = 0;i < Grid<d>::Number_Of_Nb_C();i++) {
					VectorDi nb_face = Grid<d>::Nb_C(face_cell, i);
					if (!mac_grid.Valid_Face(axis, nb_face) || bc->Is_Psi_N(axis, nb_face)) continue;
					VectorD nb_pos = mac_grid.Face_Center(axis, nb_face);
					T a = sigma * Levelset_Dirac(levelset->Phi((pos + nb_pos) * .5), dirac_band_width) * dt * dt / (mac_grid.grid.dx * mac_grid.grid.dx);
					dia_coef += a;
					if (face_idx(axis, nb_face) == -1) {
						//put to rhs
						u_old[idx] += a * velocity(axis, nb_face);
					}
					else {
						int col = face_idx(axis, nb_face);
						col_val_pairs.push_back(std::make_pair(col, -a));
					}
				}
				col_val_pairs.push_back(std::make_pair(idx, dia_coef));
				int* col = B.innerIndexPtr() + outer[idx];
				T* val = B.valuePtr() + outer[idx];
				std::sort(col_val_pairs.begin(), col_val_pairs.end());
				for (int k = 0;k < col_val_pairs.size();k++) {
					col[k] = col_val_pairs[k].first;
					val[k] = col_val_pairs[k].second;
				}
			},
			idx_coord
				);

		//timer.Record("fill matrix main part");

		SparseSolverCPX<Scalar> cpx_solver;
		cpx_solver.Init(B, 3000, 1e-5);
		cpx_solver.Solve(u_new.data(), u_old.data());
		ArrayFunc::Exec_Each(
			[&](const int idx) {
				VectorDi face_cell = idx_coord[idx];
				velocity(axis, face_cell) = u_new[idx];
			},
			idx_coord
				);

		//timer.Record("solve system");
	}
	//timer.Output_Profile(std::cout);
	timer.Elapse_And_Output_And_Reset("Implicit Surface Tension: ");
}



template class ImplicitSurfaceTension<real, 2>;
template class ImplicitSurfaceTension<real, 3>;

#endif