import sys
import os

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
vars.Add(PathVariable('datadir', 'Data file installation directory',
                      '${prefix}/share', PathVariable.PathAccept))
vars.Add(PathVariable('pythondir', 'Python module installation directory',
                      '${libdir}/python%d.%d/site-packages' \
                      % sys.version_info[0:2], PathVariable.PathAccept))
vars.Add(PathVariable('local_scratch',
                      'Local disk to use as temporary storage on each '
                      'cluster node', '/tmp', PathVariable.PathAccept))
vars.Add(PathVariable('global_scratch',
                      'Disk where job results are deposited if SCRAPP is '
                      'set in input.dat; must be on network storage '
                      '(visible to all nodes)', '/scrapp',
                      PathVariable.PathAccept))
vars.Add(PathVariable('html_coverage',
                      'Directory to output HTML coverage reports into '
                      '(requires the Python coverage module).',
                      None, PathVariable.PathIsDirCreate))

env = Environment(variables=vars)
# Inherit setup from environment (so we can do "module load modeller")
if 'PYTHONPATH' in os.environ:
    env['ENV']['PYTHONPATH'] = os.environ['PYTHONPATH']
if 'LD_LIBRARY_PATH' in os.environ:
    env['ENV']['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH']
env['ENV']['PATH'] = os.environ['PATH']

Help(vars.GenerateHelpText(env))

if not env.GetOption('clean') and not env.GetOption('help'):
    try:
        import modeller
    except ImportError, detail:
        env.Exit("Cannot import 'modeller' Python package: %s. "
                 "You may need to set PYTHONPATH, LD_LIBRARY_PATH, "
                 "and/or PATH environment variables." % detail)
    modeller_version = modeller.info.version_info
    if isinstance(modeller_version, tuple):
        modeller_version = "%s.%s" % modeller_version
    modeller_binary = "mod" + modeller_version
    modeller_path = env.WhereIs(modeller_binary)
    if not modeller_path:
        env.Exit("Could not find Modeller '%s' binary. You may need to set "
                 "the PATH environment variable." % modeller_binary)
    env.ParseConfig("%s --cflags --libs" % modeller_path)
    env.ParseConfig("pkg-config --cflags glib-2.0")
    import distutils.sysconfig
    env.Append(CPPPATH=[distutils.sysconfig.get_python_inc()])

Export('env')
SConscript('bin/SConscript')
pyso = SConscript('src/SConscript')
test = SConscript('test/SConscript')
SConscript('data/SConscript')
SConscript('lib/allosmod/SConscript')

# Tests need the C extension module built:
env.Depends(test, pyso)
