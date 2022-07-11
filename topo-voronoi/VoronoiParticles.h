#pragma once
#include "NAParticles.h"

namespace Meso{
	template<int d, class NeighborSearcherType=NeighborKDTree<d>>
	class VoronoiParticles : public NAParticles<d, NeighborSearcherType> {
		Typedef_VectorD(d); Typedef_MatrixD(d);
	public:

		VoronoiParticles() {}

		Setup_Attribute(A, MatrixD, MatrixD::Identity());
		Setup_Attribute(D, MatrixD, MatrixD::Identity());

		VoronoiParticles<d,NeighborSearcherType>& operator = (const VoronoiParticles<d, NeighborSearcherType>& voronoi_particles) {
			Shallow_Copy(voronoi_particles);
			_x = voronoi_particles._x;
			_A = voronoi_particles._A;
			_D = voronoi_particles._D;
			return *this; 
		}

		VoronoiParticles<d, NeighborSearcherType>& operator = (VoronoiParticles<d, NeighborSearcherType>& voronoi_particles) {
			Shallow_Copy(voronoi_particles);
			_x = voronoi_particles._x;
			_A = voronoi_particles._A;
			_D = voronoi_particles._D;
			return *this;
		}
	};
}