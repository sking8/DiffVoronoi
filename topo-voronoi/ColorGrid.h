//////////////////////////////////////////////////////////////////////////
// Color Grid
// Copyright (c) (2022-), Bo Zhu and Fan Feng
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"

using namespace Meso;
//color the grid using four colors in 2d or eight colors in 3d.
namespace ColorGrid
{
	template<int d> int Color(const Vector<int, d>& node)
	{
		int color=0;
		for (int i = 0; i < d; i++) { color += node[i] % 2; }
		return color;
	}

	////this function is used for generating regular domain colors
	template<int d> void Color(const Vector<int, d>& counts, Array<int>& colored_node_ptr, Array<int>& colored_node_indices)
	{
		Grid<d> grid(counts);
		int color_n = pow(2,d)-1;

		Array<Array<int> > buckets(color_n);
		for (int i = 0; i < color_n; i++)buckets[i].reserve(grid.Counts().prod() / color_n + 1);

		/*iterate_cell(iter, grid) {
			const VectorDi& cell = iter.Coord(); int c = Color(cell);
			buckets[c].push_back(grid.Cell_Index(cell));
		}*/

		grid.Iterate_Nodes(
			[&](const Vector<int,d> node) {
				int c = Color<d>(node);
				buckets[c].push_back(grid.Index(node));
			}
		);

		colored_node_indices.resize(grid.Counts().prod());

		int i = 0;
		for (int c = 0; c < color_n; c++) {
			colored_node_ptr[c] = i;
			std::copy(buckets[c].begin(), buckets[c].end(), &colored_node_indices[i]);
			i += (int)buckets[c].size();
		}
		colored_node_ptr[color_n] = i;
	}
};