#pragma once
#include "Common.h"

//////////////////////////////////////////////////////////////////////////
////AmgclSolver parameters
namespace AmgclSolver{

//// see https://amgcl.readthedocs.io/_/downloads/en/latest/pdf/
enum struct SolverType : uint8_t {
    Default = 0,
    CG = 1, //amgcl::solver::cg
    BiCGStab, // amgcl::solver::bicgstab
    BiCGStabL, // amgcl::solver::bicgstabl
    GMRES, // amgcl::solver::gmres
    LGMRES, // amgcl::solver::lgmres
    FGMRES, // amgcl::solver::fgmres
    IDRs, // amgcl::solver::idrs
    Richardson, // amgcl::solver::richardson
};

class Params
{public:
    real tolerance=(real)1e-8;
    int max_iter_num=1000;
    bool verbose=false;
    bool use_gpu=false;
    bool ns_search=false;
};

}