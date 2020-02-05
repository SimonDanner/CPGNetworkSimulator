from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import sys
import setuptools

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.
    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')

MAJOR_VERSION=1
MINOR_VERSION=0
print(str(MAJOR_VERSION) + '.' + str(MINOR_VERSION))
module1 = Extension('CPGNetworkSimulator._CPGNetworkSimulator',
                    extra_compile_args=["-std=c++14", "-Ofast", "-march=native"],
                    define_macros = [('MAJOR_VERSION', str(MAJOR_VERSION)),
                                     ('MINOR_VERSION', str(MINOR_VERSION))],
                    include_dirs = ['/usr/local/include',
                                   './include',
                                   get_pybind_include(),
                                   get_pybind_include(user=True)],
                    libraries = ['yaml-cpp'],
                    library_dirs = ['/usr/local/lib'], 
                    sources = ['./src/typedefs.cpp','./src/Network.cpp','./src/CPGNetworkSimulator.cpp','./src/CPGNetworkSimulatorPY.cpp',],
                    language = 'c++' )

class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)

setup (name = 'CPGNetworkSimulator',
        version = str(MAJOR_VERSION) + '.' + str(MINOR_VERSION),
        description = 'Simulator for neural network models of CPGs using activity-dependent population models',
        author = 'Simon Danner',
        author_email = 'simon.danner@gmail.com',
        #platforms = ['Mac OSX', 'POSIX'],
        url = 'https://github.com/SimonDanner/CPGNetworkSimulator',
        long_description = 'Simulator for neural network models using activity-dependent population models',
        ext_modules = [module1],
        install_requires=['pybind11>=2.2','numpy','scipy','matplotlib'],
        cmdclass={'build_ext': BuildExt},
        zip_safe=False,
        packages=find_packages()
    )