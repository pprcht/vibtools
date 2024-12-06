#include "vibtools.h"
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
void py_print_vib_spectrum_stdout(py::array_t<double> freq,
                                  py::array_t<double> intens) {
  // Check that the arrays are 1D and have the same size
  auto freq_buf = freq.request();
  auto intens_buf = intens.request();

  if (freq_buf.ndim != 1 || intens_buf.ndim != 1 ||
      freq_buf.shape[0] != intens_buf.shape[0]) {
    throw std::runtime_error("Invalid array dimensions: 'freq' and 'intens' "
                             "must be 1D arrays of the same size.");
  }

  int nat3 = freq_buf.shape[0];
  double *freq_ptr = static_cast<double *>(freq_buf.ptr);
  double *intens_ptr = static_cast<double *>(intens_buf.ptr);

  // Call the Fortran routine via the C wrapper
  c_print_vib_spectrum_stdout(nat3, freq_ptr, intens_ptr);
}

// Wrapper for the Lorentzian broadening routine
void py_lorentzian_broadening(int nmodes, py::array_t<double> freq,
                              py::array_t<double> intens, int npoints,
                              py::array_t<double> plt, double xmin,
                              double xmax, double dx, double fwhm) {
  // Request buffer information for the numpy arrays
  auto freq_buf = freq.request();
  auto intens_buf = intens.request();
  auto plt_buf = plt.request();

  // Single if statement for dimension checks
  if (freq_buf.ndim != 1 || intens_buf.ndim != 1 || plt_buf.ndim != 1 ||
      freq_buf.shape[0] != nmodes || intens_buf.shape[0] != nmodes ||
      plt_buf.shape[0] != npoints) {
    throw std::runtime_error(
        "Invalid array dimensions: "
        "freq and intens must be 1D arrays with length equal to nmodes, "
        "and plt must be a 1D array with length equal to npoints.");
  }

  // Convert numpy arrays to raw pointers
  double *freq_ptr = static_cast<double *>(freq_buf.ptr);
  double *intens_ptr = static_cast<double *>(intens_buf.ptr);
  double *plt_ptr = static_cast<double *>(plt_buf.ptr);

  // Call the Fortran subroutine via the C wrapper
  c_lorentzian_broadening(nmodes, freq_ptr, intens_ptr, npoints, plt_ptr,
                          xmin, xmax, dx, fwhm);
}


py::tuple py_compute_thermodynamics(int32_t nat, 
                                    py::array_t<int32_t> at,
                                    py::array_t<double> xyz,
                                    int32_t nfreq,
                                    py::array_t<double> freq,
                                    double T, double sthr, double ithr, int32_t rotnum) {

  // Request buffer information for input arrays
  auto at_buf = at.request();
  auto xyz_buf = xyz.request();
  auto freq_buf = freq.request();

  // Dimension checks for input arrays
  if (at_buf.ndim != 1 || at_buf.shape[0] != nat) {
    throw std::runtime_error("Invalid dimensions for 'at'");
  }
  if (xyz_buf.ndim != 2 || xyz_buf.shape[0] != nat || xyz_buf.shape[1] != 3) {
    throw std::runtime_error("Invalid dimensions for 'xyz'");
  }
  if (freq_buf.ndim != 1 || freq_buf.shape[0] != nfreq) {
    throw std::runtime_error("Invalid dimensions for 'freq'");
  }

  // Get pointers to input data
  int32_t *at_ptr = static_cast<int32_t *>(at_buf.ptr);
  double (*xyz_ptr)[3] = reinterpret_cast<double (*)[3]>(xyz_buf.ptr);
  double *freq_ptr = static_cast<double *>(freq_buf.ptr);

  // Prepare output scalars
  double zpve_val, et_val, ht_val, ts_val, g_val;

  // Call the Fortran routine
  c_compute_thermodynamics(nat, at_ptr, xyz_ptr, nfreq, freq_ptr,
                           T, sthr, ithr, rotnum,
                           &zpve_val, &et_val, &ht_val, &ts_val, &g_val);

  // Return as a Python tuple of scalars
  return py::make_tuple(zpve_val, et_val, ht_val, ts_val, g_val);
}




/////////////////////////////////////////////////////////////////////////////////////////
PYBIND11_MODULE(_vibtools, m) {
  m.doc() = "Python bindings for vibtools";

  m.def("py_computespec_core", &py_computespec_core,
        "Compute spectrum core function", py::arg("nat"), py::arg("at"),
        py::arg("xyz"), py::arg("hessian"), py::arg("dipole_gradient"),
        py::arg("amass"), py::arg("fscal"), py::arg("freq"), py::arg("intens"));
  m.def("py_print_vib_spectrum_stdout", &py_print_vib_spectrum_stdout,
        "Print vibrational spectrum to stdout", py::arg("freq"),
        py::arg("intens"));
  m.def("py_lorentzian_broadening", &py_lorentzian_broadening,
        "Apply Lorentzian broadening to a spectrum", py::arg("nmodes"),
        py::arg("freq"), py::arg("intens"), py::arg("npoints"),
        py::arg("plt"), py::arg("xmin"), py::arg("xmax"), py::arg("dx"),
        py::arg("fwhm"));
  m.def("py_compute_thermodynamics", &py_compute_thermodynamics,
      "Compute thermodynamic quantities, returns a tuple",
      py::arg("nat"), py::arg("at"), py::arg("xyz"), py::arg("nfreq"),
      py::arg("freq"), py::arg("T"), py::arg("sthr"), py::arg("ithr"),
      py::arg("rotnum"));

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
