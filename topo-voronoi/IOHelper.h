#pragma once
#include "Common.h"
#include "Field.h"
#include "AuxFunc.h"
#include "BoundaryConditionGrid.h"
#include <vtkNew.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

using namespace Meso;

namespace VTKFunc {
	template<class T, int d, DataHolder side>
	void Write_VTS(const Field<T, d, side>& rho, std::string file_name) { //field defined in cell
		Assert(!rho.Empty(), "VTKFunc::Output_VTS error: empty Field");
		Typedef_VectorD(d);
		Info("rho grid domain min {}:", rho.grid.Domain_Min());
		Info("rho grid domain max {}:", rho.grid.Domain_Max());
		const auto grid = Grid<d>(rho.grid.Counts()+VectorDi::Ones(), rho.grid.dx, rho.grid.Domain_Min(), CORNER);
		Info("vts grid domain min {}:", grid.Domain_Min());
		Info("vts grid domain max {}:", grid.Domain_Max());
		int nx, ny, nz;
		if constexpr (d == 2) {
			nx = grid.Counts()[0];
			ny = grid.Counts()[1];
			nz = 1;
		}
		else if constexpr (d == 3) {
			nx = grid.Counts()[0];
			ny = grid.Counts()[1];
			nz = grid.Counts()[2];
		}
		else {
			Error("Dimension is not supported!");
		}

		// setup VTK
		vtkNew<vtkXMLStructuredGridWriter> writer;
		vtkNew<vtkStructuredGrid> structured_grid;
		structured_grid->SetDimensions(nx, ny, nz);
		vtkNew<vtkPoints> nodes;
		nodes->Allocate(nx * ny * nz);
		vtkNew<vtkDoubleArray> rhoArray;
		rhoArray->SetName("rho");

		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					VectorDi cell = MathFunc::Vi<d>(i, j, k);
					Vector3 pos3 = MathFunc::V<3>(grid.Position(cell));
					nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);
					if constexpr (d == 2) {
						if ( j != (ny - 1) && i != (nx - 1)) {
							rhoArray->InsertNextTuple1(rho(cell));
						}
					}else if constexpr (d == 3) {
						if (k != (nz - 1) && j != (ny - 1) && i != (nx - 1)) {
							rhoArray->InsertNextTuple1(rho(cell));
						}
					}
				}
			}
		}

		structured_grid->SetPoints(nodes);
		structured_grid->GetCellData()->AddArray(rhoArray);
		structured_grid->GetCellData()->SetActiveScalars("rho");

		writer->SetInputData(structured_grid);
		writer->SetFileName(file_name.c_str());
		writer->SetDataModeToBinary();
		writer->Write();
	}

	template<class T, int d, DataHolder side>
	void Write_Vector_Field(const Field<Vector<T,d>, d, side>& f, std::string file_name) {
		Assert(!f.Empty(), "VTKFunc::Output_VTS error: empty Field");
		Typedef_VectorD(d);
		const auto grid = f.grid;
		int nx, ny, nz;
		if constexpr (d == 2) {
			nx = grid.Counts()[0];
			ny = grid.Counts()[1];
			nz = 1;
		}
		else if constexpr (d == 3) {
			nx = grid.Counts()[0];
			ny = grid.Counts()[1];
			nz = grid.Counts()[2];
		}
		else {
			Error("Dimension is not supported!");
		}

		// setup VTK
		vtkNew<vtkXMLStructuredGridWriter> writer;
		vtkNew<vtkStructuredGrid> structured_grid;
		structured_grid->SetDimensions(nx, ny, nz);
		vtkNew<vtkPoints> nodes;
		nodes->Allocate(nx * ny * nz);
		vtkNew<vtkDoubleArray> fArray;
		fArray->SetName("f");
		fArray->SetNumberOfComponents(3);

		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					VectorDi cell = MathFunc::Vi<d>(i, j, k);
					VectorD pos = grid.Position(cell);
					Vector3 pos3 = MathFunc::V<3>(pos);
					nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);
					Vector3 f3 = MathFunc::V<3>(f(cell));
					fArray->InsertNextTuple3(f3[0],f3[1],f3[2]);
				}
			}
		}

		structured_grid->SetPoints(nodes);
		structured_grid->GetPointData()->AddArray(fArray);
		structured_grid->GetPointData()->SetActiveScalars("f");

		writer->SetInputData(structured_grid);
		writer->SetFileName(file_name.c_str());
		writer->SetDataModeToBinary();
		writer->Write();
	}

	template<int d>
	void Write_Points(const Array<Vector<real,d>>& pos, std::string file_name) {
		Typedef_VectorD(d);

		// write psi_d
		vtkNew<vtkXMLUnstructuredGridWriter> writer;
		vtkNew<vtkUnstructuredGrid> unstructured_grid;

		vtkNew<vtkPoints> points;
		points->Allocate(pos.size());

		for (const auto& p : pos) {
			Vector3 pos3 = MathFunc::V<3>(p);
			points->InsertNextPoint(pos3[0], pos3[1], pos3[2]);
		}

		unstructured_grid->SetPoints(points);

		writer->SetFileName(file_name.c_str());
		writer->SetInputData(unstructured_grid);
		writer->Write();
	}

	template<int d>
	void Write_Boundary_Condition(const BoundaryConditionGrid<d>& bc, const Grid<d>& grid, bf::path base_path) {
		Typedef_VectorD(d);

		// write psi_d
		vtkNew<vtkXMLUnstructuredGridWriter> writer;
		vtkNew<vtkUnstructuredGrid> unstructured_grid;

		vtkNew<vtkPoints> psi_d_nodes;
		psi_d_nodes->Allocate(bc.psi_D_values.size());
		vtkNew<vtkDoubleArray> psi_d_array;
		psi_d_array->SetName("psi_d");
		psi_d_array->SetNumberOfComponents(3);

		for (const auto& b : bc.psi_D_values) {
			Vector3 pos3 = MathFunc::V<3>(grid.Position(b.first));
			psi_d_nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);
			Vector3 value3 = MathFunc::V<3>(b.second);
			psi_d_array->InsertNextTuple3(value3[0], value3[1], value3[2]);
		}

		unstructured_grid->SetPoints(psi_d_nodes);
		unstructured_grid->GetPointData()->AddArray(psi_d_array);
		unstructured_grid->GetPointData()->SetActiveVectors("psi_d");

		bf::path psi_d_path = base_path / bf::path("psi_d.vtu");
		writer->SetFileName(psi_d_path.string().c_str());
		writer->SetInputData(unstructured_grid);
		writer->SetDataModeToBinary();
		writer->Write();

		//write_psi_n
		vtkNew<vtkPoints> psi_n_nodes;
		psi_n_nodes->Allocate(bc.forces.size());
		vtkNew<vtkDoubleArray> psi_n_array;
		psi_n_array->SetName("psi_n");
		psi_n_array->SetNumberOfComponents(3);

		for (const auto& b : bc.forces) {
			Vector3 pos3 = MathFunc::V<3>(grid.Position(b.first));
			psi_n_nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);
			Vector3 value3 = MathFunc::V<3>(b.second);
			psi_n_array->InsertNextTuple3(value3[0], value3[1], value3[2]);
		}

		unstructured_grid->SetPoints(psi_n_nodes);
		unstructured_grid->GetPointData()->AddArray(psi_n_array);
		unstructured_grid->GetPointData()->SetActiveVectors("psi_n");

		bf::path psi_n_path = base_path / bf::path("psi_n.vtu");
		writer->SetFileName(psi_n_path.string().c_str());
		writer->SetInputData(unstructured_grid);
		writer->SetDataModeToBinary();
		writer->Write();
	}
}