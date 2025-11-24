#include <algorithm>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

static PyObject *do_osc(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
  {"do_osc", do_osc, METH_VARARGS, ""}, 
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef cModPyDem =
{
  PyModuleDef_HEAD_INIT,
  "oschelper",
  "",
  -1,
  module_methods
};


PyMODINIT_FUNC
PyInit_oschelper(void) {
   import_array();
   return PyModule_Create(&cModPyDem);
}

const static int NNU = 4;

static PyObject *do_osc(PyObject *self, PyObject *args) {
  // inputs
  PyObject *energy_obj;
  PyObject *steps_obj;
  PyObject *density_obj;
  PyObject *U_obj;
  PyObject *MP_obj;
  double dm21;
  double dm31;
  double dm41;
  int nuind;
  double dL;

  // output 
  PyObject *ret = NULL;

  if (!PyArg_ParseTuple(args, "OOOOOdddid", &energy_obj, &steps_obj, &density_obj, &U_obj, &MP_obj, &dm21, &dm31, &dm41, &nuind, &dL)) { 
    PyErr_SetString(PyExc_TypeError, "BAD INPUT TYPES");
    return NULL;
  }
  if (nuind < 0 || nuind > 2) {
    std::stringstream error;
    error << "ERROR BAD NUIND (" << nuind << "). Must be between 0 and 2.";
    PyErr_SetString(PyExc_ValueError, error.str().c_str());
    return NULL;
  }

  PyArrayObject *energy = (PyArrayObject *)PyArray_FROM_OTF(energy_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyArrayObject *steps = (PyArrayObject *)PyArray_FROM_OTF(steps_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyArrayObject *density = (PyArrayObject *)PyArray_FROM_OTF(density_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyArrayObject *U = (PyArrayObject *)PyArray_FROM_OTF(U_obj, NPY_CDOUBLE, NPY_IN_ARRAY);
  PyArrayObject *MP = (PyArrayObject *)PyArray_FROM_OTF(MP_obj, NPY_CDOUBLE, NPY_IN_ARRAY);
  if (energy == NULL || steps == NULL || density == NULL || U == NULL || MP == NULL) {
    Py_XDECREF(energy);
    Py_XDECREF(steps);
    Py_XDECREF(density);
    Py_XDECREF(U);
    Py_XDECREF(MP);
    PyErr_SetString(PyExc_RuntimeError, "ERROR: FAILED TO GET PYARRAY");
    return NULL;
  }

  // Check shapes
  int energy_ndim = (int) PyArray_NDIM(energy);
  int steps_ndim = (int) PyArray_NDIM(steps);
  int density_ndim = (int) PyArray_NDIM(density);
  int U_ndim = (int) PyArray_NDIM(U);
  int MP_ndim = (int) PyArray_NDIM(MP);

  if (energy_ndim != 1 || steps_ndim != 1 || density_ndim != 1) {
    Py_XDECREF(energy);
    Py_XDECREF(steps);
    Py_XDECREF(density);
    Py_XDECREF(U);
    Py_XDECREF(MP);

    std::stringstream error;
    error << "ERROR DIMS MUST BE 1. ENERGY: " << energy_ndim << " STEPS: " << steps_ndim << " DENSITY: " << density_ndim;
    PyErr_SetString(PyExc_ValueError, error.str().c_str());

    return NULL;
  }
  if (U_ndim != 2 || MP_ndim != 2) {
    Py_XDECREF(energy);
    Py_XDECREF(steps);
    Py_XDECREF(density);
    Py_XDECREF(U);
    Py_XDECREF(MP);
    std::stringstream error;
    error << "ERROR DIMS MUST BE 2. U: " << U_ndim << " MP: " << MP_ndim;
    PyErr_SetString(PyExc_ValueError, error.str().c_str());
    return NULL;
  }

  // Check sizes
  int energy_size = PyArray_SHAPE(energy)[0];
  int density_size = PyArray_SHAPE(density)[0];
  int steps_size = PyArray_SHAPE(steps)[0];
  if (steps_size != density_size+1) {
    Py_XDECREF(energy);
    Py_XDECREF(steps);
    Py_XDECREF(density);
    Py_XDECREF(U);
    Py_XDECREF(MP);

    std::stringstream error;
    error << "ERROR BAD NUMBER OF STEPS. STEPS: " << steps_size << " DENSITY: " << density_size;
    PyErr_SetString(PyExc_ValueError, error.str().c_str());

    return NULL;
  }

  npy_intp *U_dim = PyArray_SHAPE(U);
  int U_N = U_dim[0];
  int U_M = U_dim[1];
  npy_intp *MP_dim = PyArray_SHAPE(MP);
  int MP_N = MP_dim[0];
  int MP_M = MP_dim[1];

  if (U_N != NNU || U_M != NNU || MP_N != NNU || MP_M != NNU) {
    Py_XDECREF(energy);
    Py_XDECREF(steps);
    Py_XDECREF(density);
    Py_XDECREF(U);
    Py_XDECREF(MP);

    std::stringstream error;
    error << "ERROR BAD MATRIX SIZE, MUST BE " << NNU << ". U ROW " << U_N << " U COL " << U_M << " MP ROW " << MP_N << " MP COL " << MP_M;
    PyErr_SetString(PyExc_ValueError, error.str().c_str());

    return NULL;
  }

  // Get local copy of Hamiltonian
  std::complex<double> Uosc[NNU][NNU];
  for (unsigned i_A = 0; i_A < NNU; i_A++) {
    for (unsigned j_A = 0; j_A < NNU; j_A++) {
      Uosc[i_A][j_A] = *((std::complex<double>*)PyArray_GETPTR2(U, i_A, j_A));
    }
  }

  // zero-initialize output array
  npy_intp nu_dim[2] = {NNU, energy_size};
  PyArrayObject *nu = (PyArrayObject *)PyArray_SimpleNew(2, nu_dim, NPY_CDOUBLE);
  if (nu == NULL) {
    Py_XDECREF(energy);
    Py_XDECREF(steps);
    Py_XDECREF(density);
    Py_XDECREF(U);
    Py_XDECREF(MP);
    PyErr_SetString(PyExc_RuntimeError, "ERROR CREATING OUTPUT ARRAY");
    return NULL;
  }

  // Initialize neutrino state
  for (int i = 0; i < energy_size; i++) {
    *((std::complex<double> *)PyArray_GETPTR2(nu, 0, i)) = 0;
    *((std::complex<double> *)PyArray_GETPTR2(nu, 1, i)) = 0;
    *((std::complex<double> *)PyArray_GETPTR2(nu, 2, i)) = 0;
    *((std::complex<double> *)PyArray_GETPTR2(nu, 3, i)) = 0;

    *((std::complex<double> *)PyArray_GETPTR2(nu, nuind, i)) = 1;
  }

  // Time Evolution Loop!
  for (int i = 0; i < density_size; i++) { // Time / Length evolution
    double NA = 6e23; // Avogadro number, mol/g
    double Ye = 0.5; // Electron fraction, neutron fraction, proton fraction
    double rho = *((double*)PyArray_GETPTR1(density, i)); // density, g/cm3
    double hbarc_GeV_cm = 197.3e-16; // GeV*cm 
    double Gfermi = 1.166e-5; // GeV^-2
    double MP_v = Ye*rho*Gfermi*NA*(hbarc_GeV_cm*hbarc_GeV_cm*hbarc_GeV_cm)*sqrt(2)*1e9*1e9; // eV^2 * (E / GeV)

    double hbarc = 6.582119e-16*3e8; // ev*m

    for (int j = 0; j < energy_size; j++) { // Per neutrino energy

      double E = *((double*)PyArray_GETPTR1(energy, j));

      std::complex<double> A[NNU][NNU]; // Build matter potential
      for (unsigned i_A = 0; i_A < NNU; i_A++) {
        for (unsigned j_A = 0; j_A < NNU; j_A++) {
	  A[i_A][j_A] = 2*E*MP_v*(*((std::complex<double>*)PyArray_GETPTR2(MP, i_A, j_A)));
        }
      }

      double Losc = dL/hbarc; // m -> 1/eV
      double Eosc = E*1e9; // GeV -> eV

      std::complex<double> thisnu[NNU] = {*((std::complex<double>*)PyArray_GETPTR2(nu, 0, j)), 
	                          *((std::complex<double>*)PyArray_GETPTR2(nu, 1, j)),
                                  *((std::complex<double>*)PyArray_GETPTR2(nu, 2, j)),
                                  *((std::complex<double>*)PyArray_GETPTR2(nu, 3, j))};
      std::complex<double> J(0,1);
      std::complex<double> dnu[NNU] = {0, 0, 0, 0};
      dnu[0] = (-J*Losc/(2*Eosc))*(dm31*Uosc[0][2]*std::conj(Uosc[0][2]) + dm21*Uosc[0][1]*std::conj(Uosc[0][1]) + dm41*Uosc[0][3]*std::conj(Uosc[0][3]) + A[0][0])*thisnu[0] +\
               (-J*Losc/(2*Eosc))*(dm31*Uosc[0][2]*std::conj(Uosc[1][2]) + dm21*Uosc[0][1]*std::conj(Uosc[1][1]) + dm41*Uosc[0][3]*std::conj(Uosc[1][3]) + A[0][1])*thisnu[1] +\
               (-J*Losc/(2*Eosc))*(dm31*Uosc[0][2]*std::conj(Uosc[2][2]) + dm21*Uosc[0][1]*std::conj(Uosc[2][1]) + dm41*Uosc[0][3]*std::conj(Uosc[2][3]) + A[0][2])*thisnu[2] +\
               (-J*Losc/(2*Eosc))*(dm31*Uosc[0][2]*std::conj(Uosc[3][2]) + dm21*Uosc[0][1]*std::conj(Uosc[3][1]) + dm41*Uosc[0][3]*std::conj(Uosc[3][3]) + A[0][3])*thisnu[3];
        
      dnu[1] = (-J*Losc/(2*Eosc))*(dm31*Uosc[1][2]*std::conj(Uosc[0][2]) + dm21*Uosc[1][1]*std::conj(Uosc[0][1]) + dm41*Uosc[1][3]*std::conj(Uosc[0][3]) + A[1][0])*thisnu[0] +\
               (-J*Losc/(2*Eosc))*(dm31*Uosc[1][2]*std::conj(Uosc[1][2]) + dm21*Uosc[1][1]*std::conj(Uosc[1][1]) + dm41*Uosc[1][3]*std::conj(Uosc[1][3]) + A[1][1])*thisnu[1] +\
               (-J*Losc/(2*Eosc))*(dm31*Uosc[1][2]*std::conj(Uosc[2][2]) + dm21*Uosc[1][1]*std::conj(Uosc[2][1]) + dm41*Uosc[1][3]*std::conj(Uosc[2][3]) + A[1][2])*thisnu[2] +\
               (-J*Losc/(2*Eosc))*(dm31*Uosc[1][2]*std::conj(Uosc[3][2]) + dm21*Uosc[1][1]*std::conj(Uosc[3][1]) + dm41*Uosc[1][3]*std::conj(Uosc[3][3]) + A[1][3])*thisnu[3];

      dnu[2] = (-J*Losc/(2*Eosc))*(dm31*Uosc[2][2]*std::conj(Uosc[0][2]) + dm21*Uosc[2][1]*std::conj(Uosc[0][1]) + dm41*Uosc[2][3]*std::conj(Uosc[0][3]) + A[2][0])*thisnu[0] +\
               (-J*Losc/(2*Eosc))*(dm31*Uosc[2][2]*std::conj(Uosc[1][2]) + dm21*Uosc[2][1]*std::conj(Uosc[1][1]) + dm41*Uosc[2][3]*std::conj(Uosc[1][3]) + A[2][1])*thisnu[1] +\
               (-J*Losc/(2*Eosc))*(dm31*Uosc[2][2]*std::conj(Uosc[2][2]) + dm21*Uosc[2][1]*std::conj(Uosc[2][1]) + dm41*Uosc[2][3]*std::conj(Uosc[2][3]) + A[2][2])*thisnu[2] +\
               (-J*Losc/(2*Eosc))*(dm31*Uosc[2][2]*std::conj(Uosc[3][2]) + dm21*Uosc[2][1]*std::conj(Uosc[3][1]) + dm41*Uosc[2][3]*std::conj(Uosc[3][3]) + A[2][3])*thisnu[3];
        
      dnu[3] = (-J*Losc/(2*Eosc))*(dm31*Uosc[3][2]*std::conj(Uosc[0][2]) + dm21*Uosc[3][1]*std::conj(Uosc[0][1]) + dm41*Uosc[3][3]*std::conj(Uosc[0][3]) + A[3][0])*thisnu[0] +\
               (-J*Losc/(2*Eosc))*(dm31*Uosc[3][2]*std::conj(Uosc[1][2]) + dm21*Uosc[3][1]*std::conj(Uosc[1][1]) + dm41*Uosc[3][3]*std::conj(Uosc[1][3]) + A[3][1])*thisnu[1] +\
               (-J*Losc/(2*Eosc))*(dm31*Uosc[3][2]*std::conj(Uosc[2][2]) + dm21*Uosc[3][1]*std::conj(Uosc[2][1]) + dm41*Uosc[3][3]*std::conj(Uosc[2][3]) + A[3][2])*thisnu[2] +\
               (-J*Losc/(2*Eosc))*(dm31*Uosc[3][2]*std::conj(Uosc[3][2]) + dm21*Uosc[3][1]*std::conj(Uosc[3][1]) + dm41*Uosc[3][3]*std::conj(Uosc[3][3]) + A[3][3])*thisnu[3];

      *((std::complex<double>*)PyArray_GETPTR2(nu, 0, j)) += dnu[0];
      *((std::complex<double>*)PyArray_GETPTR2(nu, 1, j)) += dnu[1];
      *((std::complex<double>*)PyArray_GETPTR2(nu, 2, j)) += dnu[2];
      *((std::complex<double>*)PyArray_GETPTR2(nu, 3, j)) += dnu[3];

    }
  }


  Py_XDECREF(energy);
  Py_XDECREF(steps);
  Py_XDECREF(density);
  Py_XDECREF(U);
  Py_XDECREF(MP);

  ret = Py_BuildValue("O", nu);
  return ret;

}


