# COMPILE CYTHON

# keep html
# python3.6 setup.py build_ext --inplace

# Keep only the .c and .so
# python3.6 setup.py build_ext --inplace clean --all


from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

#Change "mult_fast" & ["mult_fast.pyx"] depends on the file name
extensions = [
    Extension(
        "fast_loop",
        ["fast_loop.pyx"],
    ),
]

# If you want html (Python-interpreter interactions with C): annotate=True
setup(
    name = "fast_loop.pyx",
    packages = find_packages(),
    ext_modules = cythonize(extensions, annotate=True)
)
