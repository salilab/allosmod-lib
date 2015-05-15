import sys

EnsureSConsVersion(0, 98)

vars = Variables('config.py', ARGUMENTS)
vars.Add(PathVariable('prefix', 'Top-level installation directory', '/usr',
                      PathVariable.PathAccept))
vars.Add(PathVariable('destdir',
                      'String to prepend to every installed filename',
                      '', PathVariable.PathAccept))
vars.Add(PathVariable('libdir', 'Shared library installation directory',
                      '${prefix}/lib', PathVariable.PathAccept))
vars.Add(PathVariable('bindir', 'Binary installation directory',
                      '${prefix}/bin', PathVariable.PathAccept))
vars.Add(PathVariable('pythondir', 'Python module installation directory',
                      '${libdir}/python%d.%d/site-packages' \
                      % sys.version_info[0:2], PathVariable.PathAccept))
vars.Add(PathVariable('html_coverage',
                      'Directory to output HTML coverage reports into '
                      '(requires the Python coverage module).',
                      None, PathVariable.PathIsDirCreate))

env = Environment(variables=vars)

Help(vars.GenerateHelpText(env))

Export('env')
SConscript('bin/SConscript')
SConscript('test/SConscript')
SConscript('lib/allosmod/SConscript')
