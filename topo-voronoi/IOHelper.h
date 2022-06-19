#pragma once
#include "Common.h"
#include "Field.h"
#include "AuxFunc.h"
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
	void Write_VTS(const Field<T, d, side>& rho, std::string file_name) {
		Assert(!rho.Empty(), "VTKFunc::Output_VTS error: empty Field");
		Typedef_VectorD(d);
		const auto grid = rho.grid;
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
					VectorD pos = grid.Position(cell);
					Vector3 pos3 = MathFunc::V<3>(pos);
					nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);
					rhoArray->InsertNextTuple1(rho(cell));
				}
			}
		}

		structured_grid->SetPoints(nodes);
		structured_grid->GetPointData()->AddArray(rhoArray);
		structured_grid->GetPointData()->SetActiveScalars("rho");

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
}