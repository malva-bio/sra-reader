"""
Compatibility shim — delegates to build.py.

Use ``python build.py build_ext --inplace`` for in-place Cython compilation,
or ``pip install .`` / ``poetry install`` for a full install.
"""

from setuptools import setup
from build import get_extensions, BuildExt

setup(
    name="sra-reader",
    packages=["sra_reader"],
    ext_modules=get_extensions(),
    cmdclass=dict(build_ext=BuildExt),
)
