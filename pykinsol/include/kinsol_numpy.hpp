#pragma once
#include <Python.h>
#include <numpy/arrayobject.h>
#include <chrono>
// #include <utility> // std::pair
#include <vector> // std::vector

#include "kinsol_cxx.hpp"
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */


namespace kinsol_numpy{

    using SVectorView = sundials_cxx::nvector_serial::VectorView;

    class PyKinsol {
    public:
        PyObject *py_func, *py_jac;
        const int nu;
        const int ml, mu;

        PyKinsol(PyObject * py_func, PyObject * py_jac, size_t nu, int ml=-1, int mu=-1) :
            py_func(py_func), py_jac(py_jac), nu(nu), ml(ml), mu(mu) {}

        PyObject * solve(PyObject *py_x0, double fnormtol, double scsteptol, long int mxiter,
                         PyObject *py_x_scale, PyObject *py_f_scale, PyObject * py_constraints){
            std::clock_t cputime0 = std::clock();
            auto solver = kinsol_cxx::Solver();
            solver.init(kinsol_cxx::f_cb<PyKinsol>, this->nu);
            solver.set_user_data(static_cast<void*>(this));
            if (ml == -1 && mu == -1){
                solver.set_linear_solver_to_dense(this->nu);
                solver.set_dense_jac_fn(kinsol_cxx::jac_dense_cb<PyKinsol>);
            } else {
                solver.set_linear_solver_to_banded(this->nu, this->mu, this->ml);
                solver.set_band_jac_fn(kinsol_cxx::jac_band_cb<PyKinsol>);
            }
            solver.set_num_max_iters(mxiter);
            solver.set_func_norm_tol(fnormtol);
            solver.set_scaled_steptol(scsteptol);
            solver.set_constraints(SVectorView(this->nu, static_cast<double*>(PyArray_GETPTR1(py_constraints, 0))));
            int flag = solver.solve(SVectorView(this->nu, static_cast<double*>(PyArray_GETPTR1(py_x0, 0))),
                                    KIN_LINESEARCH,
                                    SVectorView(this->nu, static_cast<double*>(PyArray_GETPTR1(py_x_scale, 0))),
                                    SVectorView(this->nu, static_cast<double*>(PyArray_GETPTR1(py_f_scale, 0))));
            PyObject *d = PyDict_New();
            // naming scheme from: scipy.optimize.OptimizeResult
            PyDict_SetItemString(d, "x", py_x0);
            PyDict_SetItemString(d, "success", (flag >= 0) ?
                                 (Py_INCREF(Py_True), Py_True) :
                                 (Py_INCREF(Py_False), Py_False));
            PyDict_SetItemString(d, "status", PyInt_FromLong(flag));
            PyDict_SetItemString(d, "nfev", PyInt_FromLong(solver.get_num_func_evals()));
            PyDict_SetItemString(d, "njev", PyInt_FromLong(solver.get_num_jac_evals()));
            PyDict_SetItemString(d, "nit", PyInt_FromLong(solver.get_num_nonlin_solv_iters()));
            PyDict_SetItemString(d, "time_cpu", PyFloat_FromDouble((std::clock() - cputime0) / (double)CLOCKS_PER_SEC));
            switch(flag){
            case KIN_SUCCESS:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_SUCCESS")); break;
            case KIN_INITIAL_GUESS_OK:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_INITIAL_GUESS_OK")); break;
            case KIN_STEP_LT_STPTOL:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_STEP_LT_STPTOL")); break;
            case KIN_LINESEARCH_NONCONV:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_LINESEARCH_NONCONV")); break;
            case KIN_MAXITER_REACHED:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_MAXITER_REACHED")); break;
            case KIN_MXNEWT_5X_EXCEEDED:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_MXNEWT_5X_EXCEEDED")); break;
            case KIN_LINESEARCH_BCFAIL:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_LINESEARCH_BCFAIL")); break;
            case KIN_LINSOLV_NO_RECOVERY:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_LINSOLV_NO_RECOVERY")); break;
            case KIN_LINIT_FAIL:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_LINIT_FAIL")); break;
            case KIN_LSETUP_FAIL:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_LSETUP_FAIL")); break;
            case KIN_LSOLVE_FAIL:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_LSOLVE_FAIL")); break;
            case KIN_MEM_NULL:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_MEM_NULL")); break;
            case KIN_NO_MALLOC:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_NO_MALLOC")); break;
            case KIN_ILL_INPUT:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "KIN_ILL_INPUT")); break;
            default:
                PyDict_SetItemString(d, "message", Py_BuildValue("s", "Unknown exit status.")); break;
            }
            return d;
        }

        void func(const double * const u, double * const fval){
            npy_intp dims[1] { static_cast<npy_intp>(this->nu) } ;
            PyObject * py_uarr = PyArray_SimpleNewFromData(
                1, dims, NPY_DOUBLE, static_cast<void*>(const_cast<double*>(u)));
            PyObject * py_fval = PyArray_SimpleNewFromData(
                1, dims, NPY_DOUBLE, static_cast<void*>(fval));
            PyObject * py_arglist = Py_BuildValue("(OO)", py_uarr, py_fval);
            PyObject * py_result = PyEval_CallObject(this->py_func, py_arglist);
            Py_DECREF(py_arglist);
            Py_DECREF(py_fval);
            Py_DECREF(py_uarr);
            if (py_result == nullptr){
                throw std::runtime_error("func() failed");
            } else if (py_result != Py_None){
                // py_result is not None
                throw std::runtime_error("func() failed (returned other than None)");
            }
            Py_DECREF(py_result);
        }
        void call_py_jac(const double * const u, const double * const fy,
                         PyObject * py_jmat){
            npy_intp udims[1] { static_cast<npy_intp>(this->nu) };
            PyObject * py_uarr = PyArray_SimpleNewFromData(1, udims, NPY_DOUBLE,
                const_cast<double *>(u));
            PyObject * py_fy = PyArray_SimpleNewFromData(1, udims, NPY_DOUBLE,
                                                         const_cast<double *>(fy));
            PyObject * py_arglist = Py_BuildValue("(OOO)", py_uarr, py_jmat, py_fy);
            PyObject * py_result = PyEval_CallObject(this->py_jac, py_arglist);
            Py_DECREF(py_arglist);
            Py_DECREF(py_fy);
            Py_DECREF(py_uarr);
            if (py_result == nullptr){
                throw std::runtime_error("jac() failed");
            } else if (py_result != Py_None){
                throw std::runtime_error("jac() did not return None");
            }
            Py_DECREF(py_result);
        }
        void dense_jac_cmaj(const double * const u, const double * const fu,
                            double * const jac, long int ldim){
            npy_intp Jdims[2] { static_cast<npy_intp>(this->nu), static_cast<npy_intp>(this->nu) };
            npy_intp strides[2] { sizeof(double), static_cast<npy_intp>(ldim*sizeof(double)) };
            PyObject * py_jmat = PyArray_New(
                &PyArray_Type, 2, Jdims, NPY_DOUBLE, strides,
                static_cast<void *>(const_cast<double *>(jac)), sizeof(double),
                NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_WRITEABLE, nullptr);
            call_py_jac(u, fu, py_jmat);
            Py_DECREF(py_jmat);
        }
        void banded_padded_jac_cmaj(const double * const u, const double * const fu,
                                    double * const jac, long int ldim){
            npy_intp Jdims[2] { 1 + this->ml + this->mu, static_cast<npy_intp>(this->nu) };
            npy_intp strides[2] { sizeof(double), static_cast<npy_intp>(ldim*sizeof(double)) };
            PyObject * py_jmat = PyArray_New(
                &PyArray_Type, 2, Jdims, NPY_DOUBLE, strides,
                static_cast<void *>(const_cast<double *>(jac + this->mu)), sizeof(double),
                NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_WRITEABLE, nullptr);
            call_py_jac(u, fu, py_jmat);
            Py_DECREF(py_jmat);
        }
    };

} // namespace kinsol_numpy
