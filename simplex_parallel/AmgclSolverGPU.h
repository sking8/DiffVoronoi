//////////////////////////////////////////////////////////////////////////
// Geometric Multigrid PCG GPU
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
//// This is the GPU implementation of the geometric multigrid class
//////////////////////////////////////////////////////////////////////////

#ifdef USE_CUDA
#ifndef __AmgclSolverGPU_h__
#define __AmgclSolverGPU_h__

#include "Timer.h"
#include "AmgclParam.h"
#include "ContextCuda.h"
#include "AuxFuncCuda.h"
#include "AuxFunc.h"

#include <amgcl/adapter/eigen.hpp>
#include <amgcl/backend/cuda.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>

//////////////////////////////////////////////////////////////////////////
//// GPU amgcl single function API
template<int d> bool AMGCL_GPU(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const AmgclSolver::Params params/*=Params()*/,bool verbose=false);

#endif
#endif