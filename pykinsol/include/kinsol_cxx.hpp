#pragma once

#include <sundials/sundials_config.h>
#include "sundials_cxx.hpp" // sundials_cxx::nvector_serial::Vector
#include <kinsol/kinsol_direct.h>
#include <kinsol/kinsol.h>
#if SUNDIALS_VERSION_MAJOR >= 3
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunlinsol/sunlinsol_lapackdense.h>
#include <sunlinsol/sunlinsol_lapackband.h>
#else
#if defined(SUNDIALS_PACKAGE_VERSION)   /* == 2.7.0 */
#include <kinsol/kinsol_lapack.h>
#else
#error "Unkown sundials version"
#endif
#endif


namespace kinsol_cxx {

    using SVector = sundials_cxx::nvector_serial::Vector; // serial vector
    using SVectorView = sundials_cxx::nvector_serial::VectorView; // serial vector

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
#if SUNDIALS_VERSION_MAJOR >= 3
        SUNMatrix A_ = nullptr;
        SUNLinearSolver LS_ = nullptr;
        N_Vector y_ = nullptr;
#endif
    public:
        void *mem {nullptr};
        Solver() {
            this->mem = KINCreate();
        }
        ~Solver(){
            if (this->mem)
                KINFree(&(this->mem));
#if SUNDIALS_VERSION_MAJOR >= 3
            if (this->y_)
                N_VDestroy(this->y_);
            if (this->LS_)
                SUNLinSolFree(this->LS_);
            if (this->A_)
                SUNMatDestroy(this->A_);
#endif
        }
        // init
        void init(KINSysFn cb, N_Vector y) {
#if SUNDIALS_VERSION_MAJOR >= 3
            if (y_)
                throw std::runtime_error("y_ already allocated");
            y_ = N_VNew_Serial(NV_LENGTH_S(y));
            std::memcpy(NV_DATA_S(y_), NV_DATA_S(y), NV_LENGTH_S(y)*sizeof(realtype));
#else
            auto y_ = y;
#endif
            int status = KINInit(this->mem, cb, y_);
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
        int solve(SVectorView u, int strategy, SVector u_scale, SVector f_scale){
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
            int flag;
#if SUNDIALS_VERSION_MAJOR >= 3
            if (A_ == nullptr){
                if (A_)
                    throw std::runtime_error("matrix already set");
                A_ = SUNDenseMatrix(ny, ny);
                if (!A_)
                    throw std::runtime_error("SUNDenseMatrix failed.");
            }
            if (LS_ == nullptr){
                if (LS_)
                    throw std::runtime_error("linear solver already set");
                LS_ = SUNLapackDense(y_, A_);
                if (!LS_)
                    throw std::runtime_error("SUNDenseLinearSolver failed.");
            }
            flag = KINDlsSetLinearSolver(this->mem, LS_, A_);
            if (flag < 0)
                throw std::runtime_error("KINDlsSetLinearSolver failed.");
#else
            flag = KINLapackDense(this->mem, ny);
            if (flag != KINDLS_SUCCESS)
                throw std::runtime_error("KINLapackDense failed");
#endif
        }
        void set_dense_jac_fn(
#if SUNDIALS_VERSION_MAJOR >= 3
                              KINDlsJacFn
#else
                              KINDlsDenseJacFn
#endif

                              djac){
            int flag;
#if SUNDIALS_VERSION_MAJOR >= 3
            flag = KINDlsSetJacFn(this->mem, djac);
#else
            flag = KINDlsSetDenseJacFn(this->mem, djac);
#endif
            if (flag < 0)
                throw std::runtime_error("set_dense_jac_fn failed.");
        }

        // banded jacobian
        void set_linear_solver_to_banded(int ny, int mupper, int mlower){
            int flag;
#if SUNDIALS_VERSION_MAJOR >= 3
            if (A_ == nullptr){
                if (A_)
                    throw std::runtime_error("matrix already set");
                A_ = SUNBandMatrix(ny, mupper, mlower, std::min(ny-1, mlower+mupper));
                if (!A_)
                    throw std::runtime_error("SUNDenseMatrix failed.");
            }
            if (LS_ == nullptr){
                if (LS_)
                    throw std::runtime_error("linear solver already set");
                LS_ = SUNLapackBand(y_, A_);
                if (!LS_)
                    throw std::runtime_error("SUNDenseLinearSolver failed.");
            }
            flag = KINDlsSetLinearSolver(this->mem, LS_, A_);
            if (flag < 0)
                throw std::runtime_error("KINDlsSetLinearSolver failed.");
#else
            flag = KINLapackBand(this->mem, ny, mupper, mlower);
            if (flag != KINDLS_SUCCESS)
                throw std::runtime_error("KINLapackBand failed");
#endif
        }
        void set_band_jac_fn(
#if SUNDIALS_VERSION_MAJOR >= 3
                              KINDlsJacFn
#else
                              KINDlsBandJacFn
#endif
                              djac){
            int flag;
#if SUNDIALS_VERSION_MAJOR >= 3
            flag = KINDlsSetJacFn(this->mem, djac);
#else
            flag = KINDlsSetBandJacFn(this->mem, djac);
#endif
            if (flag < 0)
                throw std::runtime_error("set_band_jac_fn failed.");
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
    int jac_dense_cb(
#if SUNDIALS_VERSION_MAJOR < 3
                     long int N,
#endif
                     N_Vector u, N_Vector fu,
#if SUNDIALS_VERSION_MAJOR < 3
                     DlsMat Jac,
#else
                     SUNMatrix Jac,
#endif
                     void *user_data,
                     N_Vector tmp1, N_Vector tmp2){ // KINDlsDenseJacFn
        NeqSys * neqsys = (NeqSys*)user_data;
        neqsys->dense_jac_cmaj(NV_DATA_S(u), NV_DATA_S(fu),
#if SUNDIALS_VERSION_MAJOR < 3
                               DENSE_COL(Jac, 0), Jac->ldim
#else
                               SM_DATA_D(Jac), NV_LENGTH_S(u)
#endif
                               );
        return 0;
    }

    template <typename NeqSys>
    int jac_band_cb(
#if SUNDIALS_VERSION_MAJOR < 3
                    long int /* N */, long int mupper, long int mlower,
#endif
                    N_Vector u, N_Vector fu,
#if SUNDIALS_VERSION_MAJOR < 3
                     DlsMat Jac,
#else
                     SUNMatrix Jac,
#endif
                    void *user_data,
                    N_Vector tmp1, N_Vector tmp2){
        // callback of req. signature wrapping Neqsys method.
        NeqSys * neqsys = (NeqSys*)user_data;
#if SUNDIALS_VERSION_MAJOR < 3
        if (neqsys->mu != mupper)
            throw std::runtime_error("mupper mismatch");
        if (neqsys->ml != mlower)
            throw std::runtime_error("mlower mismatch");
        auto Jac_ = Jac;
#else
        auto Jac_ = SM_CONTENT_B(Jac);
#endif

        neqsys->banded_padded_jac_cmaj(NV_DATA_S(u), NV_DATA_S(fu), Jac_->data, Jac_->ldim);
        return 0;
    }

} // namespace cvodes_cxx
