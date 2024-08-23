#include "irtools.h"
#include <iostream>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

// Wrapper for the c_computespec_core function
void py_computespec_core(int32_t nat, py::array_t<int32_t> at,
                         py::array_t<double> xyz, py::array_t<double> hess,
                         py::array_t<double> dipd, py::array_t<double> ams,
                         double fscal, py::array_t<double> freq,
                         py::array_t<double> ints) {

  // Ensure the arrays are contiguous and have the correct dimensions
  auto at_buf = at.request();
  auto xyz_buf = xyz.request();
  auto hess_buf = hess.request();
  auto dipd_buf = dipd.request();
  auto ams_buf = ams.request();
  auto freq_buf = freq.request();
  auto ints_buf = ints.request();

  if (at_buf.ndim != 1 || xyz_buf.ndim != 2 || xyz_buf.shape[1] != 3 ||
      hess_buf.ndim != 2 || dipd_buf.ndim != 2 || dipd_buf.shape[1] != 3 ||
      ams_buf.ndim != 1 || ams_buf.shape[0] != 118 || freq_buf.ndim != 1 ||
      ints_buf.ndim != 1) {
    throw std::runtime_error("Invalid array dimensions");
  }

  // Get pointers to the data
  int32_t *at_ptr = static_cast<int32_t *>(at_buf.ptr);
  double(*xyz_ptr)[3] = reinterpret_cast<double(*)[3]>(xyz_buf.ptr);
  double *hess_ptr = static_cast<double *>(hess_buf.ptr);
  double(*dipd_ptr)[3] = reinterpret_cast<double(*)[3]>(dipd_buf.ptr);
  double *ams_ptr = static_cast<double *>(ams_buf.ptr);
  double *freq_ptr = static_cast<double *>(freq_buf.ptr);
  double *ints_ptr = static_cast<double *>(ints_buf.ptr);

  // Call the Fortran routine
  c_computespec_core(nat, at_ptr, xyz_ptr, hess_ptr, dipd_ptr, ams_ptr, fscal,
                     freq_ptr, ints_ptr);
}


// Wrapper for the print_vib_spectrum_stdout function
void py_print_vib_spectrum_stdout(py::array_t<double> freq, py::array_t<double> intens) {
    // Check that the arrays are 1D and have the same size
    auto freq_buf = freq.request();
    auto intens_buf = intens.request();

    if (freq_buf.ndim != 1 || intens_buf.ndim != 1 || freq_buf.shape[0] != intens_buf.shape[0]) {
        throw std::runtime_error("Invalid array dimensions: 'freq' and 'intens' must be 1D arrays of the same size.");
    }

    int nat3 = freq_buf.shape[0];
    double *freq_ptr = static_cast<double *>(freq_buf.ptr);
    double *intens_ptr = static_cast<double *>(intens_buf.ptr);

    // Call the Fortran routine via the C wrapper
    c_print_vib_spectrum_stdout(nat3, freq_ptr, intens_ptr);
}


/////////////////////////////////////////////////////////////////////////////////////////
PYBIND11_MODULE(_irtools, m) {
  m.doc() = "Python bindings for IRtools";

  m.def("py_computespec_core", &py_computespec_core,"Compute spectrum core function",
        py::arg("nat"),py::arg("at"),py::arg("xyz"),py::arg("hessian"), 
        py::arg("dipole_gradient"),py::arg("amass"),py::arg("fscal"),
        py::arg("freq"),py::arg("intens"));
  m.def("py_print_vib_spectrum_stdout", &py_print_vib_spectrum_stdout, 
        "Print vibrational spectrum to stdout", py::arg("freq"),py::arg("intens"));


#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
