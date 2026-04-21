"""
Cython build script for sra-reader.

Invoked automatically by Poetry (via [tool.poetry.build] script = "build.py").
Can also be run directly:  python build.py build_ext --inplace
"""

import os
import sys
import numpy
from distutils.command.build_ext import build_ext
from setuptools import Extension
from Cython.Build import cythonize


# Set SRA_READER_DEBUG_BUILD=1 to build with profiling hooks and debug symbols.
DEBUG_BUILD = os.environ.get('SRA_READER_DEBUG_BUILD', '0') == '1'


class BuildExt(build_ext):
    """Silently skip failed extensions rather than aborting the whole build."""
    def build_extensions(self):
        try:
            super().build_extensions()
        except Exception:
            pass


def get_extensions():
    if DEBUG_BUILD:
        compiler_directives = {
            'language_level': 3,
            'boundscheck': False,
            'wraparound': False,
            'initializedcheck': False,
            'cdivision': True,
            'linetrace': True,
            'profile': True,
        }
        compile_args = ["-O0", "-g"]
        link_args = []
        macros = [('CYTHON_TRACE', '1'), ('CYTHON_TRACE_NOGIL', '1')]
    else:
        compiler_directives = {
            'language_level': 3,
            'boundscheck': False,
            'wraparound': False,
            'initializedcheck': False,
            'cdivision': True,
            'linetrace': False,
            'profile': False,
            'embedsignature': False,
            'emit_code_comments': False,
        }
        compile_args = [
            "-O3",
            "-ffast-math",
            "-DNDEBUG",
            "-fvisibility=hidden",
            "-ffunction-sections",
            "-fdata-sections",
        ]
        macros = [('NDEBUG', '1')]

        if sys.platform == 'linux':
            compile_args.append("-fno-semantic-interposition")
            link_args = [
                "-Wl,--strip-all",
                "-Wl,--gc-sections",
                "-Wl,--exclude-libs,ALL",
                "-Wl,-z,relro",
                "-Wl,-z,now",
            ]
        else:
            link_args = ["-Wl,-dead_strip"]

    extensions = [
        Extension(
            "sra_reader._parser",
            ["sra_reader/_parser.pyx"],
            include_dirs=[numpy.get_include()],
            extra_compile_args=compile_args,
            extra_link_args=link_args,
            define_macros=macros,
        ),
    ]

    return cythonize(extensions,
                     compiler_directives=compiler_directives,
                     force=True)


def build(setup_kwargs):
    """Called by Poetry's build backend."""
    setup_kwargs.update(
        dict(
            cmdclass=dict(build_ext=BuildExt),
            ext_modules=get_extensions(),
        )
    )


if __name__ == "__main__":
    from setuptools import setup
    setup(
        name="sra-reader",
        ext_modules=get_extensions(),
        cmdclass=dict(build_ext=BuildExt),
    )
