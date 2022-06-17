#pragma once
#include "Common.h"
#include "Constants.h"
#include "SparseFunc.h"

using namespace Meso;

namespace LinearFEMFunc
{
	////////////////////////////////////////////////////////////////////////
	// E: strain -> stress
	template<int d> void Strain_Stress_Matrix_Linear(const real youngs, const real poisson, MatrixX& E);

	////////////////////////////////////////////////////////////////////////
	//Hex element
	template<int d> void dNde(const Vector<real, d>& natural_coord, MatrixX& dNde);							////dNde, return a dxn matrix, n is the number of basis functions
	template<int d> Vector<real, d> dNde(const Vector<real, d>& natural_coord, const int idx);				////dNde for the idx'th basis function, return a d-dim vector

	template<int d> void Cell_dNdX(const Vector<real, d>& natural_coord, const real dx, MatrixX& dNdX);		////dNdX

	template<int d> real Cell_Strain_Displacement_Matrix(const Vector<real, d>& natural_coord, const real dx, MatrixX& B);			//B: displacement->strain
	template<int d> void Cell_Stiffness_Matrix(const real youngs, const real poisson, const real dx, MatrixX& K_e);

	//Hex helper functions
	template<int d> void Add_Cell_Stiffness_Matrix(SparseMatrix<real>& K, const MatrixX& K_e, const Array<int>& nodes);
	template<int d> void Set_Cell_B_Elements_Helper(const int r, const int c, const VectorX& dN, MatrixX& B);

	//Gaussian integration
	template<int d> void Initialize_Gaussian_Integration_Points(ArrayF2P<Vector<real, d>, d>& points, ArrayF2P<real, d>& weights); // Gauss product rule p = 2, 3rd order accuracy

	////////////////////////////////////////////////////////////////////////
	//Helper functions
	void Set_Dirichlet_Boundary_Helper(SparseMatrix<real>& K, VectorX& b, const int i, const real psi_D_value);
};
