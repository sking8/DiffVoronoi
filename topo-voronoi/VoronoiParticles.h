#pragma once
#include "Points.h"

namespace Meso{
	template<int d>
	class VoronoiParticles : public Points {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:

		VoronoiParticles() {}

		Setup_Attribute(x, VectorD, VectorD::Zero());
		Setup_Attribute(A, MatrixD, MatrixD::Identity());
		Setup_Attribute(D, MatrixD, MatrixD::Identity());
	};
}