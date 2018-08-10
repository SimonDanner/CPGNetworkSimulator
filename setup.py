from distutils.core import setup, Extension

MAJOR_VERSION=1
MINOR_VERSION=1
print(str(MAJOR_VERSION) + '.' + str(MINOR_VERSION))
module1 = Extension('CPGNetworkSimulator',
                    extra_compile_args=["-std=c++14", "-Ofast", "-march=native"],
                    define_macros = [('MAJOR_VERSION', str(MAJOR_VERSION)),
                                     ('MINOR_VERSION', str(MINOR_VERSION))],
                    include_dirs = ['/usr/local/include','./include'],
                    libraries = ['yaml-cpp'],
                    library_dirs = ['/usr/local/lib'], 
                    sources = ['./src/typedefs.cpp','./src/Network.cpp','./src/CPGNetworkSimulator.cpp','./src/CPGNetworkSimulatorPY.cpp',],
                    language = 'c++' )

setup (name = 'CPGNetworkSimulator',
       version = str(MAJOR_VERSION) + '.' + str(MINOR_VERSION),
       description = 'Simulator for neural network models using activity-dependent population models',
       author = 'Simon Danner',
       author_email = 'simon.danner@gmail.com',
       #platforms = ['Mac OSX', 'POSIX'],
       url = '',
       long_description = 'Simulator for neural network models using activity-dependent population models',
       ext_modules = [module1])