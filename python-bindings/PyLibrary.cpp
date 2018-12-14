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

    m.def("inversion_error",&inversion_error, "Inversion with measurement errors", py::arg("filepath"), py::arg("norm_id"),
            py::arg("discretization_steps")=10000, py::arg("nu")=-1);


    py::enum_<Norms>(m, "Norm")
            .value("L2ErrorNorm", Norms::L2ErrorNorm)
            .value("SemiErrorNorm", Norms::SemiErrorNorm);

    py::class_<ModelParameters>(m, "ModelParameters")
            .def_readonly("nu", &ModelParameters::nu)
            .def_readonly("misfit", &ModelParameters::misfit)
            .def_readonly("norm", &ModelParameters::norm)
            .def("__repr__", &ModelParameters_repr);
}

