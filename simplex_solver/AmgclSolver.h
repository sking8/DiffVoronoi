//////////////////////////////////////////////////////////////////////////
// AMGCL Solver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __AmgclSolver_h__
#define __AmgclSolver_h__

#include <numeric>
#include <functional>
#include "AmgclParam.h"
#include "Timer.h"

//// include AMGCL headers
#include "assert.h"
#include <amgcl/amg.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/adapter/eigen.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_builder.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/profiler.hpp>

#ifdef USE_CUDA
#include "AmgclSolverGPU.h"
#endif

//// see https://github.com/ddemidov/amgcl/issues/103
//// to enable Eigen matrix with different backend
AMGCL_USE_EIGEN_VECTORS_WITH_BUILTIN_BACKEND()

namespace AmgclSolver{

template<int d> bool AMGPCG(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Params params/*=Params()*/,bool verbose/*=false*/)
{
#ifdef USE_CUDA
    if(params.use_gpu) return AMGCL_GPU<d>(A,x,b,params,true);
#endif
    Timer timer;timer.Reset();
    // Setup the solver:
    typedef amgcl::make_solver<
        amgcl::amg<
            amgcl::backend::builtin<real>,
            amgcl::coarsening::smoothed_aggregation,
            amgcl::relaxation::spai0
            >,
        amgcl::solver::bicgstab<amgcl::backend::builtin<real> >
        > Solver;

    Solver::params prm;
    prm.solver.maxiter=params.max_iter_num;
    prm.solver.tol=params.tolerance;
    // prm.solver.verbose=params.verbose;
    prm.solver.ns_search=params.ns_search;

    Solver solve(A,prm);
    if(verbose){timer.Elapse_And_Output_And_Reset("CPU AMGCL: allocate A");}
    if(verbose){std::cout<<"CPU AMGCL: solver info:"<<std::endl<<solve<<std::endl;}

    // Solve the system for the given RHS:
    int iters;
    real error;
    timer.Reset();
    std::tie(iters, error) = solve(b, x);
    if(verbose)timer.Elapse_And_Output_And_Reset("CPU AMGCL: cpu solving");
    if(verbose)std::cout <<"#Iters="<< iters << ", error=" << error << std::endl;
    return 0;
}

};
#endif