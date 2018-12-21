#include "../libGravimetryInversion/Library.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <string>
#include <sstream>

namespace py = pybind11;

std::string ModelParameters_repr(const ModelParameters& mp) {
    std::ostringstream oss;
    oss << "ModelParameters(nu=" << mp.nu << ", misfit=" << mp.misfit << ", norm=" << mp.norm << ")";
    return oss.str();
}


    PYBIND11_MODULE(PyLibrary, m){
    m.doc() = "Pybind11 example plugin";

    // export inversion function
    m.def("inversion", &inversion, "Inversion with measurement errors", py::arg("filepath"), py::arg("norm_id"),
            py::arg("discretization_steps")=10000);
    // export inversion with errors function
    m.def("inversion_error", &inversion_error, "Inversion with measurement errors", py::arg("filepath"), py::arg("norm_id"),
            py::arg("discretization_steps")=10000, py::arg("nu")=-1);
    // export interpolation function
    m.def("interpolate", &interpolation, "Interpolation between points", py::arg("filepath"), py::arg("norm_id"),
            py::arg("lower"), py::arg("upper"), py::arg("discretization_steps")=10000);


    // export enum used for measurement error free inversion
    py::enum_<Norms>(m, "Norm")
            .value("L2Norm", Norms::L2Norm)
            .value("W12Norm", Norms::W12Norm)
            .value("SemiNorm", Norms::SemiNorm);
    // export enum used for norms that take measurement error into account
    py::enum_<ErrorNorms >(m, "ErrorNorm")
            .value("L2ErrorNorm", ErrorNorms ::L2ErrorNorm)
            .value("SemiErrorNorm", ErrorNorms ::SemiErrorNorm);
    // export enum used for interpolation norms
    py::enum_<InterpolationNorms>(m, "InterpolationNorm")
            .value("LinearInterpolationNorm", InterpolationNorms::LinearInterpolationNorm);
    // ModelParameter class saves parameters that were determined/used during inversion
    py::class_<ModelParameters>(m, "ModelParameters")
            .def_readonly("nu", &ModelParameters::nu)
            .def_readonly("misfit", &ModelParameters::misfit)
            .def_readonly("norm", &ModelParameters::norm)
            .def("__repr__", &ModelParameters_repr);
}

