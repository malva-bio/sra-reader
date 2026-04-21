from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension(
        "sra_reader._parser",
        ["sra_reader/_parser.pyx"],
        extra_compile_args=["-O3", "-march=native"],
    ),
]

setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
            "initializedcheck": False,
        },
    ),
)
