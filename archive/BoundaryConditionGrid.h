#pragma once
#include "Hashtable.h"
using namespace Meso;

template<int d> class BoundaryConditionGrid
{
	Typedef_VectorD(d);
public:
	Hashtable<VectorDi, VectorD> psi_D_values;		////fixed points
	Hashtable<VectorDi, VectorD> forces;			////force

	BoundaryConditionGrid() {}
	void Set_Psi_D(const VectorDi& node, const VectorD& displacement)
	{
		psi_D_values[node] = displacement;
	}
	void Set_Force(const VectorDi& node, const VectorD& force)
	{
		forces[node] = force;
	}
	void Clear() { psi_D_values.clear(); forces.clear(); }
};