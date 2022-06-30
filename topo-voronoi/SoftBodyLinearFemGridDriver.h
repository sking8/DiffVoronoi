//#####################################################################
// Soft body FEM driver
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __SoftBodyLinearFemGridDriver_h__
#define __SoftBodyLinearFemGridDriver_h__
#include <fstream>
#include "SPX_Common.h"
#include "SPX_Driver.h"
#include "SoftBodyLinearFemGrid.h"

template<int d> class SoftBodyLinearFemGridDriver : public Driver
{Typedef_VectorDii(d);Typedef_MatrixD(d);using Base=Driver;
public:
	Grid<d> grid;
	SoftBodyLinearFemGrid<d> soft_body;

	void Initialize()
	{
		Initialize_Grid();
		Initialize_Materials();
		Initialize_Boundary_Conditions();
	}

	void Initialize_Grid()
	{
		real dx=(real)1/(real)scale;VectorDi cell_counts=VectorDi::Ones()*scale;
		grid.Initialize(cell_counts,dx);
		soft_body.Initialize(grid);
	}

	void Initialize_Materials()
	{
		soft_body.Clear_Materials();
		soft_body.Add_Material((real)1,(real).35);
		soft_body.material_id.Fill(0);
	}

	void Initialize_Boundary_Conditions()
	{
		switch(test){
		case 1:{
			int axis=0;VectorD dis=VectorD::Unit(axis)*(real).1;
			iterate_node(iter,grid){const VectorDi& node=iter.Coord();
				if(node[axis]==0)soft_body.Set_Fixed(node);
				else if(node[axis]==grid.node_counts[axis]-1)soft_body.Set_Displacement(node,dis);}
		}break;
		case 2:{
			VectorD force=VectorD::Unit(1)*(real)-.1;
			iterate_node(iter,grid){const VectorDi& node=iter.Coord();
				if(node[0]==0)soft_body.Set_Fixed(node);
				else if(node[1]==0&&node[0]==grid.node_counts[0]-1)soft_body.Set_Force(node,force);}
		}break;
		}
	}

	void Run()
	{
		soft_body.Allocate_K();
		soft_body.Update_K_And_f();
		soft_body.Solve();

		Write_Output_Files(0);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		if(frame==0){
			std::string file_name=frame_dir+"/grid";
			soft_body.grid.Write_To_File_3d(file_name);
			std::cout<<"Write to file "<<file_name<<std::endl;}
		
		{std::string file_name=frame_dir+"/displacement";
		Write_To_Field_3d<VectorX,real,d>(soft_body.u,soft_body.grid.node_counts,file_name);}

		{std::string file_name=frame_dir+"/force";
		Write_To_Field_3d<VectorX,real,d>(soft_body.f,soft_body.grid.node_counts,file_name);}

		{std::string file_name=frame_dir+"/strain";
		Field<MatrixD,d> strain;soft_body.Compute_Strain(strain);
		strain.Write_To_File_3d(file_name);
		file_name=frame_dir+"/stress";
		Field<MatrixD,d> stress;soft_body.Compute_Stress(stress);
		stress.Write_To_File_3d(file_name);
		Field<real,d> von_mises;soft_body.Compute_Von_Mises_Stress(von_mises);
		file_name=frame_dir+"/von_mises_stress";
		von_mises.Write_To_File_3d(file_name);}
	}
};

#endif