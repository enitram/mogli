########################################################################################################################
#     mogli - molecular graph library                                                                                  #
#                                                                                                                      #
#     Copyright (C) 2016-2019  Martin S. Engler                                                                        #
#                                                                                                                      #
#     This program is free software: you can redistribute it and/or modify                                             #
#     it under the terms of the GNU Lesser General Public License as published                                         #
#     by the Free Software Foundation, either version 3 of the License, or                                             #
#     (at your option) any later version.                                                                              #
#                                                                                                                      #
#     This program is distributed in the hope that it will be useful,                                                  #
#     but WITHOUT ANY WARRANTY; without even the implied warranty of                                                   #
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                     #
#     GNU General Public License for more details.                                                                     #
#                                                                                                                      #
#     You should have received a copy of the GNU Lesser General Public License                                         #
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.                                           #
########################################################################################################################

import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      '-DBUILD_TESTS=OFF']

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.', '--target', 'mogli'] + build_args, cwd=self.build_temp)


setup(
    name='mogli',
    version='1.0',
    author='Martin S. Engler',
    author_email='martin.engler@hhu.de',
    description='The molecular graph library.',
    license='GNU Lesser General Public License v3 or later (LGPLv3+)',
    ext_modules=[CMakeExtension('mogli')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'],
)