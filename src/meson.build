project('soap_turbo_interface', 'c', 'fortran')

py = import('python').find_installation()
fc = meson.get_compiler('fortran')

# --------------------------------------------------------------------
# Build the static Fortran library (replaces `make`)
# --------------------------------------------------------------------
fortran_sources = files(
  'soap_turbo/src/mod_types.f90',
  'soap_turbo/src/soap_turbo_functions.f90',
  'soap_turbo/src/soap_turbo_compress.f90',
  'soap_turbo/src/soap_turbo_radial.f90',
  'soap_turbo/src/soap_turbo_angular.f90',
  'soap_turbo/src/soap_turbo.f90',
)

libsoap_turbo = static_library(
  'soap_turbo_fortran',
  sources: fortran_sources,
  fortran_args: ['-O3', '-fPIC'],
  install: false,
)

# --------------------------------------------------------------------
# Query the active Python for paths at configure time
# --------------------------------------------------------------------
ext_suffix = run_command(
  py, '-c', 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"), end="")',
  check: true,
).stdout()

py_include = run_command(
  py, '-c', 'import sysconfig; print(sysconfig.get_path("include"), end="")',
  check: true,
).stdout()

numpy_include = run_command(
  py, '-c', 'import numpy; print(numpy.get_include(), end="")',
  check: true,
).stdout()

f2py_include = run_command(
  py, '-c', 'import numpy.f2py, os; print(os.path.join(os.path.dirname(numpy.f2py.__file__), "src"), end="")',
  check: true,
).stdout()

py_install_dir = run_command(
  py, '-c', 'import site; print(site.getsitepackages()[0], end="")',
  check: true,
).stdout()

f2py_f90wrap = find_program('f2py-f90wrap')

# --------------------------------------------------------------------
# Step A: f90wrap — generates six Fortran wrapper .f90 files
# --------------------------------------------------------------------
f90wrap_fortran = custom_target(
  'f90wrap_fortran',
  input : files(
    'soap_turbo/src/mod_types.f90',
    'soap_turbo/src/soap_turbo_functions.f90',
    'soap_turbo/src/soap_turbo_compress.f90',
    'soap_turbo/src/soap_turbo_radial.f90',
    'soap_turbo/src/soap_turbo_angular.f90',
    'soap_turbo/src/soap_turbo.f90',
  ),
  output : [
    'f90wrap_mod_types.f90',
    'f90wrap_soap_turbo_functions.f90',
    'f90wrap_soap_turbo_compress.f90',
    'f90wrap_soap_turbo_radial.f90',
    'f90wrap_soap_turbo_angular.f90',
    'f90wrap_soap_turbo.f90',
  ],
  command : [
    py, '-m', 'f90wrap',
    '@INPUT@',
    '-m', 'soap_turbo',
    '-k', meson.project_source_root() / 'kind_map.json',
  ],
  install : false,
)

# Run f90wrap again to get soap_turbo.py, then patch it in the same command:
# replace `import _soap_turbo` with `from . import _soap_turbo` so it works
# as a submodule inside the soap_turbo package.
soap_turbo_py_patched = custom_target(
  'f90wrap_py_patched',
  input : files(
    'soap_turbo/src/mod_types.f90',
    'soap_turbo/src/soap_turbo_functions.f90',
    'soap_turbo/src/soap_turbo_compress.f90',
    'soap_turbo/src/soap_turbo_radial.f90',
    'soap_turbo/src/soap_turbo_angular.f90',
    'soap_turbo/src/soap_turbo.f90',
  ),
  output : ['soap_turbo_ext.py'],
  command : [
    py, '-c',
    'import subprocess, sys, os; subprocess.run([sys.executable, "-m", "f90wrap"] + sys.argv[3:] + ["-m", "soap_turbo", "-k", sys.argv[1]], check=True, cwd=sys.argv[2]); src=open(os.path.join(sys.argv[2],"soap_turbo.py")).read(); open(os.path.join(sys.argv[2],"soap_turbo_ext.py"),"w").write(src.replace("import _soap_turbo","from . import _soap_turbo"))',
    meson.project_source_root() / 'kind_map.json',
    '@OUTDIR@',
    '@INPUT@',
  ],
  install : true,
  install_dir : py_install_dir / 'soap_turbo',
)

# --------------------------------------------------------------------
# Step B: f2py-f90wrap — generates the C glue with PyInit__soap_turbo
# --------------------------------------------------------------------
f2py_c = custom_target(
  'f2py_c',
  input : [
    f90wrap_fortran[0], f90wrap_fortran[1], f90wrap_fortran[2],
    f90wrap_fortran[3], f90wrap_fortran[4], f90wrap_fortran[5],
  ],
  output : ['_soap_turbomodule.c'],
  command : [
    f2py_f90wrap,
    '-m', '_soap_turbo',
    '--lower',
    '@INPUT@',
    '--build-dir', '@OUTDIR@',
  ],
  install : false,
)

# --------------------------------------------------------------------
# Step C: compile everything into _soap_turbo.so
# --------------------------------------------------------------------
custom_target(
  '_soap_turbo_ext',
  input  : [
    f90wrap_fortran[0], f90wrap_fortran[1], f90wrap_fortran[2],
    f90wrap_fortran[3], f90wrap_fortran[4], f90wrap_fortran[5],
    f2py_c[0],
    f2py_include / 'fortranobject.c',
  ],
  output : ['_soap_turbo' + ext_suffix],
  depends: [libsoap_turbo],
  command : [
    'gfortran',
    '-shared', '-fPIC',
    '-I' + meson.current_build_dir() / 'libsoap_turbo_fortran.a.p',
    '-I' + py_include,
    '-I' + numpy_include,
    '-I' + f2py_include,
    '@INPUT@',
    libsoap_turbo.full_path(),
    '-L/usr/lib/x86_64-linux-gnu',
    '-llapack', '-lblas', '-lgfortran', '-lm',
    '-o', '@OUTPUT@',
  ],
  install : true,
  install_dir : py_install_dir / 'soap_turbo',
)

# --------------------------------------------------------------------
# Install the Python package files
# --------------------------------------------------------------------
install_data(
  'descriptor.py',
  install_dir : py_install_dir / 'soap_turbo',
)

install_data(
  '__init__.py',
  install_dir : py_install_dir / 'soap_turbo',
)
