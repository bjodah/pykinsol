#pragma once

// Thin C++11 wrapper around KINSOL v2.8.2 (SUNDIALS v2.6.2)
// far from all functionality is available yet.
// sundials-2.6.2.tar.gz (MD5: 3deeb0ede9f514184c6bd83ecab77d95)

// #include <assert.h>
// #include <cmath>
// #include <cstring>
// #include <memory>
// #include <utility>
// #include <vector>

#include "sundials_cxx.hpp" // sundials_cxx::nvector_serial::Vector
#include <kinsol/kinsol.h>
#include <kinsol/kinsol_direct.h>
#include <kinsol/kinsol_lapack.h>

#ifndef NDEBUG
#include <iostream> // DEBUG
#endif


namespace kinsol_cxx {

    using SVector = sundials_cxx::nvector_serial::Vector; // serial vector

    // enum class Strategy : int {NONE=KIN_NONE, LINESEARCH=KIN_LINESEARCH, FP=KIN_FP, PICARD=KIN_PICARD};
    // enum class EtaChoice : int {C1=KIN_ETACHOICE1, C2=KIN_ETACHOICE2, CONSTANT=KIN_ETACONSTANT}

    inline void check_flag(int flag) {
        switch (flag){
        case KIN_SUCCESS:
            break;
        case KIN_MEM_NULL:
            throw std::runtime_error("kin_mem is NULL");
        default:
            throw std::runtime_error("Unhandled flag");
        }
    }

    class Solver{ // Thin wrapper class of KINSolver in KINODES
    public:
        void *mem {nullptr};
        Solver() {
            this->mem = KINCreate();
        }
        ~Solver(){
            if (this->mem)
                KINFree(&(this->mem));
        }
        // init
        void init(KINSysFn cb, N_Vector tmpl) {
            int status = KINInit(this->mem, cb, tmpl);
            if (status == KIN_ILL_INPUT)
                throw std::runtime_error("KINInit failed (KIN_ILL_INPUT).");
            else if (status == KIN_MEM_FAIL)
                throw std::runtime_error("KINInit failed (allocation failed).");
            else
                check_flag(status);
        }
        void init(KINSysFn cb, SVector tmpl) {
            this->init(cb, tmpl.n_vec);
        }
        int solve(SVector u, int strategy, SVector u_scale, SVector f_scale){
            return KINSol(this->mem, u.n_vec, strategy, u_scale.n_vec, f_scale.n_vec);
        }
        void check_solve_flag(int flag, bool steptol_fail=true) {
            switch (flag){
            case KIN_SUCCESS:
                break;
            case KIN_INITIAL_GUESS_OK:
                break;
            case KIN_STEP_LT_STPTOL:
                if (steptol_fail)
                    throw std::runtime_error("KIN_STEP_LT_STPTOL");
                break;
            case KIN_MEM_NULL:
                throw std::runtime_error("KIN_MEM_NULL");
            case KIN_ILL_INPUT:
                throw std::runtime_error("KIN_ILL_INPUT");
            case KIN_NO_MALLOC:
                throw std::runtime_error("KIN_NO_MALLOC");
            case KIN_LINESEARCH_NONCONV:
                throw std::runtime_error("KIN_LINESEARCH_NONCONV");
            case KIN_MAXITER_REACHED:
                throw std::runtime_error("KIN_MAXITER_REACHED");
            case KIN_MXNEWT_5X_EXCEEDED:
                throw std::runtime_error("KIN_MXNEWT_5X_EXCEEDED");
            case KIN_LINESEARCH_BCFAIL:
                throw std::runtime_error("KIN_LINESEARCH_BCFAIL");
            case KIN_LINSOLV_NO_RECOVERY:
                throw std::runtime_error("KIN_LINSOLV_NO_RECOVERY");
            case KIN_LINIT_FAIL:
                throw std::runtime_error("KIN_LINIT_FAIL");
            case KIN_LSETUP_FAIL:
                throw std::runtime_error("KIN_LSETUP_FAIL");
            case KIN_LSOLVE_FAIL:
                throw std::runtime_error("KIN_LSOLVE_FAIL");
            case KIN_SYSFUNC_FAIL:
                throw std::runtime_error("KIN_SYSFUNC_FAIL");
            case KIN_FIRST_SYSFUNC_ERR:
                throw std::runtime_error("KIN_FIRST_SYSFUNC_ERR");
            case KIN_REPTD_SYSFUNC_ERR:
                throw std::runtime_error("KIN_REPTD_SYSFUNC_ERR");
            default:
                throw std::runtime_error("Unhandled flag");
            };
        }


        void set_num_max_iters(long int mxiter){
            int status = KINSetNumMaxIters(this->mem, mxiter);
            if (status == KIN_ILL_INPUT)
                throw std::runtime_error("The maximum number of iterations was non-positive.");
            check_flag(status);
        }
        void set_max_setup_calls(long int msbset){
            int status = KINSetMaxSetupCalls(this->mem, msbset);
            if (status == KIN_ILL_INPUT)
                throw std::runtime_error("The argument msbset was negative.");
            else
                check_flag(status);
        }
        void set_max_sub_setup_calls(long int msbsetsub){
            int status = KINSetMaxSubSetupCalls(this->mem, msbsetsub);
            if (status == KIN_ILL_INPUT)
                throw std::runtime_error("The argument msbsetsub was negative.");
            else
                check_flag(status);
        }
        void set_eta_form(int etachoice) {
            int status = KINSetEtaForm(this->mem, etachoice);
            if (status == KIN_ILL_INPUT)
                throw std::runtime_error("The argument etachoice had an illegal value.");
            else
                check_flag(status);
        }
        void set_eta_const_value(realtype eta) {
            int flag = KINSetEtaConstValue(this->mem, eta);
            if (flag == KIN_ILL_INPUT)
                throw std::runtime_error("The argument eta had an illegal value. (0 < eta <= 1)");
            else
                check_flag(flag);
        }
        void set_eta_params(realtype egamma, realtype ealpha) {
            int flag = KINSetEtaParams(this->mem, egamma, ealpha);
            if (flag == KIN_ILL_INPUT)
                throw std::runtime_error("One of the arguments egamma or ealpha had an illegal"
                                         " value. (0 < egamma <= 1, 1 < ealpha <= 2)");
            else
                check_flag(flag);
        }
        void set_res_mon_const_value(realtype omegaconst){
            int flag = KINSetResMonConstValue(this->mem, omegaconst);
            if (flag == KIN_ILL_INPUT)
                throw std::runtime_error("The argument omegaconst had an illegal value. (0 < omegaconst < 1)");
            else
                check_flag(flag);
        }
        void set_res_mon_params(realtype omegamin, realtype omegamax){
            int flag = KINSetResMonParams(this->mem, omegamin, omegamax);
            if (flag == KIN_ILL_INPUT)
                throw std::runtime_error("One of the arguments omegamin or omegamax had an illegal"
                                         " value. (0 < omegamin < omegamax < 1)");
            else
                check_flag(flag);
        }
        void set_func_norm_tol(realtype fnormtol){
            int flag = KINSetFuncNormTol(this->mem, fnormtol);
            if (flag == KIN_ILL_INPUT)
                throw std::runtime_error("The tolerance was negative.");
            else
                check_flag(flag);
        }
        void set_scaled_steptol(realtype scsteptol){
            int flag = KINSetScaledStepTol(this->mem, scsteptol);
            if (flag == KIN_ILL_INPUT)
                throw std::runtime_error("The tolerance was non-positive.");
            else
                check_flag(flag);
        }
        void set_constraints(SVector constraints){
            int flag = KINSetConstraints(this->mem, constraints.n_vec);
            if (flag == KIN_ILL_INPUT)
                throw std::runtime_error("The constraint vector contains illegal values.");
            else
                check_flag(flag);
        }


        // user data
        void set_user_data(void *user_data){
            int flag = KINSetUserData(this->mem, user_data);
            if (flag < 0)
                throw std::runtime_error("KINSetUserData failed.");
        }

        // dense jacobian
        void set_linear_solver_to_dense(int ny){
            int flag = KINLapackDense(this->mem, ny);
            if (flag != KINDLS_SUCCESS)
                throw std::runtime_error("KINLapackDense failed");
        }
        void set_dense_jac_fn(KINDlsDenseJacFn djac){
            int flag = KINDlsSetDenseJacFn(this->mem, djac);
            if (flag < 0)
                throw std::runtime_error("KINDlsSetDenseJacFn failed.");
        }

        // banded jacobian
        void set_linear_solver_to_banded(int N, int mupper, int mlower){
            int flag = KINLapackBand(this->mem, N, mupper, mlower);
            if (flag != KINDLS_SUCCESS)
                throw std::runtime_error("KINLapackBand failed");
        }
        void set_band_jac_fn(KINDlsBandJacFn djac){
            int flag = KINDlsSetBandJacFn(this->mem, djac);
            if (flag < 0)
                throw std::runtime_error("KINDlsSetBandJacFn failed.");
        }

        // Optional output
        long int get_num_func_evals(){
            long int res=0;
            int flag = KINGetNumFuncEvals(this->mem, &res);
            check_flag(flag);
            return res;
        }
        long int get_num_nonlin_solv_iters(){
            long int res=0;
            int flag = KINGetNumNonlinSolvIters(this->mem, &res);
            check_flag(flag);
            return res;
        }
        // KINDLS linear solvers
        long int get_num_jac_evals(){
            long int res=0;
            int flag = KINDlsGetNumJacEvals(this->mem, &res);
            check_flag(flag);
            return res;
        }
    };

    template<class NeqSys> // Matches KINSysFn, wraps NeqSys which with method: func(double *, double *)
    int f_cb(N_Vector u, N_Vector fval, void *user_data){
        NeqSys * neqsys = static_cast<NeqSys*>(user_data);
        neqsys->func(NV_DATA_S(u), NV_DATA_S(fval));
        return 0;
    }

    template <class NeqSys>
    int jac_dense_cb(long int N, N_Vector u, N_Vector fu, DlsMat Jac, void *user_data,
                     N_Vector tmp1, N_Vector tmp2){ // KINDlsDenseJacFn
        NeqSys * neqsys = (NeqSys*)user_data;
        neqsys->dense_jac_cmaj(NV_DATA_S(u), NV_DATA_S(fu), DENSE_COL(Jac, 0),
                               Jac->ldim);
        return 0;
    }

    template <typename NeqSys>
    int jac_band_cb(long int N, long int mupper, long int mlower,
                    N_Vector u, N_Vector fu, DlsMat Jac, void *user_data,
                    N_Vector tmp1, N_Vector tmp2){
        // callback of req. signature wrapping Neqsys method.
        NeqSys * neqsys = (NeqSys*)user_data;
        if (neqsys->mu != mupper)
            throw std::runtime_error("mupper mismatch");
        if (neqsys->ml != mlower)
            throw std::runtime_error("mlower mismatch");
        neqsys->banded_padded_jac_cmaj(NV_DATA_S(u), NV_DATA_S(fu), Jac->data, Jac->ldim);
        return 0;
    }

} // namespace cvodes_cxx
