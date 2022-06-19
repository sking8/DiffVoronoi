#include "Grid.h"
#include "AuxFunc.h"
#include "LinearFEMFunc.h"

using namespace Meso;
namespace LinearFEMFunc {
	////////////////////////////////////////////////////////////////////////
	//Material model
	template<> void Strain_Stress_Matrix_Linear<2>(const real youngs, const real poisson, MatrixX& E)
	{
		////zero stress in z
		E.resize(3, 3); E.fill((real)0); real e = youngs / ((real)1 - pow(poisson, 2)); real G = youngs / ((real)2 * ((real)1 + poisson));
		E(0, 0) = E(1, 1) = e; E(0, 1) = E(1, 0) = e * poisson; E(2, 2) = G;
	}

	template<> void Strain_Stress_Matrix_Linear<3>(const real youngs, const real poisson, MatrixX& E)
	{
		E.resize(6, 6); E.fill((real)0); real e = youngs / (((real)1 - (real)2 * poisson) * ((real)1 + poisson)); real G = youngs / ((real)2 * ((real)1 + poisson));
		E(0, 0) = E(1, 1) = E(2, 2) = e * ((real)1 - poisson); E(0, 1) = E(1, 0) = E(0, 2) = E(2, 0) = E(1, 2) = E(2, 1) = e * poisson; E(3, 3) = E(4, 4) = E(5, 5) = G;
	}

	////////////////////////////////////////////////////////////////////////
	//Hex element
	template<> void dNde<2>(const Vector2& natural_coord, MatrixX& dNde)
	{
		const real s = natural_coord[0]; const real t = natural_coord[1];
		dNde << t - 1, -t - 1, -t + 1, t + 1,
			s - 1, -s + 1, -s - 1, s + 1;
		dNde *= (real).25;
	}

	template<> void dNde<3>(const Vector3& natural_coord, MatrixX& dNde)
	{
		const real s = natural_coord[0]; const real t = natural_coord[1]; const real q = natural_coord[2];
		dNde << -(1 - t) * (1 - q), -(1 - t) * (1 + q), -(1 + t) * (1 - q), -(1 + t) * (1 + q), (1 - t)* (1 - q), (1 - t)* (1 + q), (1 + t)* (1 - q), (1 + t)* (1 + q),
			-(1 - s) * (1 - q), -(1 - s) * (1 + q), (1 - s)* (1 - q), (1 - s)* (1 + q), -(1 + s) * (1 - q), -(1 + s) * (1 + q), (1 + s)* (1 - q), (1 + s)* (1 + q),
			-(1 - s) * (1 - t), (1 - s)* (1 - t), -(1 - s) * (1 + t), (1 - s)* (1 + t), -(1 + s) * (1 - t), (1 + s)* (1 - t), -(1 + s) * (1 + t), (1 + s)* (1 + t);
		dNde *= (real).125;
	}

	template<> Vector2 dNde<2>(const Vector2& natural_coord, const int idx)
	{
		const real s = natural_coord[0]; const real t = natural_coord[1];
		switch (idx) {
		case 0:return (real).25 * Vector2(t - 1, s - 1);
		case 1:return (real).25 * Vector2(-t - 1, -s + 1);
		case 2:return (real).25 * Vector2(-t + 1, -s - 1);
		case 3:return (real).25 * Vector2(t + 1, s + 1);
		default:return Vector2::Zero();
		}
	}

	template<> Vector3 dNde<3>(const Vector3& natural_coord, const int idx)
	{
		const real s = natural_coord[0]; const real t = natural_coord[1]; const real q = natural_coord[2];
		switch (idx) {
		case 0:return (real).125 * Vector3(-(1 - t) * (1 - q), -(1 - s) * (1 - q), -(1 - s) * (1 - t));
		case 1:return (real).125 * Vector3(-(1 - t) * (1 + q), -(1 - s) * (1 + q), (1 - s) * (1 - t));
		case 2:return (real).125 * Vector3(-(1 + t) * (1 - q), (1 - s) * (1 - q), -(1 - s) * (1 + t));
		case 3:return (real).125 * Vector3(-(1 + t) * (1 + q), (1 - s) * (1 + q), (1 - s) * (1 + t));
		case 4:return (real).125 * Vector3((1 - t) * (1 - q), -(1 + s) * (1 - q), -(1 + s) * (1 - t));
		case 5:return (real).125 * Vector3((1 - t) * (1 + q), -(1 + s) * (1 + q), (1 + s) * (1 - t));
		case 6:return (real).125 * Vector3((1 + t) * (1 - q), (1 + s) * (1 - q), -(1 + s) * (1 + t));
		case 7:return (real).125 * Vector3((1 + t) * (1 + q), (1 + s) * (1 + q), (1 + s) * (1 + t));
		default:return Vector3::Zero();
		}
	}

	template<int d> void Cell_dNdX(const Vector<real, d>& natural_coord, const real dx, MatrixX& dNdX)
	{
		dNde<d>(natural_coord, dNdX); dNdX *= ((real).5 / dx);
	}

	template<> real Cell_Strain_Displacement_Matrix<2>(const Vector2& natural_coord, const real dx, MatrixX& B)
	{
		const int d = 2; int tensor_n = d * (d + 1) / 2; int r = tensor_n; int vtx_n = (int)pow(2, d); int c = d * vtx_n;
		B.resize(r, c); B.fill(0);
		MatrixX dNdX(d, vtx_n); Cell_dNdX<d>(natural_coord, dx, dNdX);
		real J_det = pow((real).5 * dx, d);

		const int x = 0; const int y = 1;
		Set_Cell_B_Elements_Helper<d>(0, x, dNdX.row(0), B);
		Set_Cell_B_Elements_Helper<d>(1, y, dNdX.row(1), B);
		Set_Cell_B_Elements_Helper<d>(2, x, dNdX.row(1), B);
		Set_Cell_B_Elements_Helper<d>(2, y, dNdX.row(0), B);

		return J_det;
	}

	template<> real Cell_Strain_Displacement_Matrix<3>(const Vector3& natural_coord, const real dx, MatrixX& B)
	{
		const int d = 3; int tensor_n = d * (d + 1) / 2; int r = tensor_n; int vtx_n = (int)pow(2, d); int c = d * vtx_n;
		B.resize(r, c); B.fill(0);
		MatrixX dNdX(d, vtx_n); Cell_dNdX<d>(natural_coord, dx, dNdX);
		real J_det = pow((real).5 * dx, d);

		const int x = 0; const int y = 1; const int z = 2;
		Set_Cell_B_Elements_Helper<d>(0, x, dNdX.row(0), B);
		Set_Cell_B_Elements_Helper<d>(1, y, dNdX.row(1), B);
		Set_Cell_B_Elements_Helper<d>(2, z, dNdX.row(2), B);
		Set_Cell_B_Elements_Helper<d>(3, x, dNdX.row(1), B);
		Set_Cell_B_Elements_Helper<d>(3, y, dNdX.row(0), B);
		Set_Cell_B_Elements_Helper<d>(4, y, dNdX.row(2), B);
		Set_Cell_B_Elements_Helper<d>(4, z, dNdX.row(1), B);
		Set_Cell_B_Elements_Helper<d>(5, x, dNdX.row(2), B);
		Set_Cell_B_Elements_Helper<d>(5, z, dNdX.row(0), B);

		return J_det;
	}

	template<int d> void Cell_Stiffness_Matrix(const real youngs, const real poisson, const real dx, MatrixX& K_e)
	{
		int n = d * (int)pow(2, d); K_e.resize(n, n); K_e.fill((real)0);
		ArrayF2P<Vector<real, d>, d> points; ArrayF2P<real, d> weights; Initialize_Gaussian_Integration_Points<d>(points, weights);
		for (auto i = 0; i < points.size(); i++) {
			MatrixX B; real J = Cell_Strain_Displacement_Matrix<d>(points[i], dx, B);
			MatrixX K0; MatrixX E; Strain_Stress_Matrix_Linear<d>(youngs, poisson, E); K0 = B.transpose() * E * B * J * weights[i];
			K_e += K0;
		}
	}
	template void Cell_Stiffness_Matrix<2>(const real youngs, const real poisson, const real dx, MatrixX& K_e);
	template void Cell_Stiffness_Matrix<3>(const real youngs, const real poisson, const real dx, MatrixX& K_e);

	template<int d> void Set_Cell_B_Elements_Helper(const int r, const int c, const VectorX& dN, MatrixX& B)
	{
		for (int i = 0; i < (int)dN.size(); i++)B(r, c + i * d) = dN[i];
	}

	////////////////////////////////////////////////////////////////////////
	//Gaussian integration
	template<> void Initialize_Gaussian_Integration_Points<2>(ArrayF2P<Vector2, 2>& points, ArrayF2P<real, 2>& weights) //2^2=4 sample points
	{
		real c = (real)1 / sqrt((real)3); points[0] = Vector2(-c, -c); points[1] = Vector2(-c, c); points[2] = Vector2(c, -c); points[3] = Vector2(c, c); ArrayFunc::Fill(weights, (real)1);
	}

	template<> void Initialize_Gaussian_Integration_Points<3>(ArrayF2P<Vector3, 3>& points, ArrayF2P<real, 3>& weights) //2^3=8 sample points
	{
		real c = (real)1/sqrt((real)3); points[0] = Vector3(-c, -c, -c); points[1] = Vector3(-c, -c, c); points[2] = Vector3(-c, c, -c); points[3] = Vector3(-c, c, c);
		points[4] = Vector3(c, -c, -c); points[5] = Vector3(c, -c, c); points[6] = Vector3(c, c, -c); points[7] = Vector3(c, c, c); ArrayFunc::Fill(weights, (real)1);
	}

	////////////////////////////////////////////////////////////////////////
	//Operations on the global stiffness matrix
	template<int d> void Add_Cell_Stiffness_Matrix(/*rst*/SparseMatrix<real>& K, const MatrixX& K_e, const Array<int>& nodes)
	{
		const int node_n = (int)nodes.size(); for (int Ke_i = 0; Ke_i < node_n; Ke_i++) {
			int K_i = nodes[Ke_i]; for (int Ke_j = 0; Ke_j < node_n; Ke_j++) { int K_j = nodes[Ke_j]; SparseFunc::Add_Block<d>(K, K_i, K_j, K_e, Ke_i, Ke_j); }
		}
	}
	template void Add_Cell_Stiffness_Matrix<2>(/*rst*/SparseMatrix<real>& K, const MatrixX& K_e, const Array<int>& nodes);
	template void Add_Cell_Stiffness_Matrix<3>(/*rst*/SparseMatrix<real>& K, const MatrixX& K_e, const Array<int>& nodes);

	void Set_Dirichlet_Boundary_Helper(SparseMatrix<real>& K, VectorX& b, const int i, const real psi_D_value)
	{
		for (InnerIterator<real> iter(K, i); iter; ++iter) {
			int j = (int)iter.col();
			if (i == j) { b(j) = psi_D_value * K.coeff(i, j); }
			else {
				real K_ij = K.coeff(i, j); b(j) -= K_ij * psi_D_value;
				K.coeffRef(i, j) = (real)0; K.coeffRef(j, i) = (real)0;
			}
		}
	}
}
