from setuptools import setup, Extension
import numpy

setup(
  name="oschelper",
  ext_modules=[Extension("oschelper", ["src/oschelper.C"])],
  include_dirs=[numpy.get_include()],
  version="0.1",
  description="Helper functions for Neutrino Oscillation Interpolation",
  author="Gray Putnam",
  # packages=["oschelper"],
  install_requires=["numpy"],
  zip_safe=False
)
