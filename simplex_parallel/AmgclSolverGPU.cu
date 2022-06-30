//////////////////////////////////////////////////////////////////////////
// Auxiliary Function CUDA
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include "AmgclSolverGPU.h"

//////////////////////////////////////////////////////////////////////////
//// single function API
////gpu amgcl
template<int d> bool AMGCL_GPU(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const AmgclSolver::Params params/*=Params()*/,bool verbose/*=false*/)
{
#ifndef USE_CUDA
	std::cerr<<"Error: [AMGCL_GPU] USE_CUDA disabled"<<std::endl;
	return false;
#else
    // ref to https://buildmedia.readthedocs.org/media/pdf/amgcl/latest/amgcl.pdf, Listing 2.6
	Timer timer;timer.Reset();
    typedef amgcl::make_solver<
        amgcl::amg<
            amgcl::backend::cuda<real>,
            amgcl::coarsening::smoothed_aggregation,
            amgcl::relaxation::spai0
            >,
        amgcl::solver::bicgstab<
            amgcl::backend::cuda<real>
            >
        > CudaSolver;
    
    amgcl::backend::cuda<real>::params bprm;
    cusparseCreate(&bprm.cusparse_handle);

	if(!Is_Cuda_Context_Initialized())Initialize_Cuda_Context();///one bug

	int n=A.rows();
    thrust::device_vector<real> b_dev(b.data(), b.data() + n);
    thrust::device_vector<real> x_dev(x.data(), x.data() + n);
	if(verbose){timer.Elapse_And_Output_And_Reset("GPU AMGCL: cpu to gpu transfer");}

    CudaSolver::params prm;
    prm.solver.maxiter=params.max_iter_num;
    prm.solver.tol=params.tolerance;
    prm.solver.ns_search=params.ns_search;

    CudaSolver solve(A, prm, bprm); // Eigen::SparseMatrix is used directly
	if(verbose){timer.Elapse_And_Output_And_Reset("GPU AMGCL: solver setup");}

    int iters;
    double error;
    std::tie(iters, error) = solve(b_dev, x_dev);
	if(verbose)timer.Elapse_And_Output_And_Reset("GPU AMGCL: gpu solving");
    if(verbose)std::cout <<"#Iters="<< iters << ", error=" << error << std::endl;

    thrust::copy(x_dev.begin(),x_dev.end(),x.data());
	if(verbose){timer.Elapse_And_Output_And_Reset("GPU AMGCL: gpu to cpu transfer");AuxFunc::Seperation_Line();}

	return true;
#endif
}


template bool AMGCL_GPU<2>(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const AmgclSolver::Params params/*=Params()*/,bool verbose/*=false*/);
template bool AMGCL_GPU<3>(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const AmgclSolver::Params params/*=Params()*/,bool verbose/*=false*/);