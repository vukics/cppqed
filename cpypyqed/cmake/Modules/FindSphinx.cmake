# Arch linux hack
if(EXISTS /etc/arch-release)
  set(sphinx sphinx-build2)
else()
  set(sphinx sphinx-build)
endif()

find_program(SPHINX_EXECUTABLE NAMES ${sphinx}
    HINTS
    $ENV{SPHINX_DIR}
    PATH_SUFFIXES bin
    DOC "Sphinx documentation generator"
)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Sphinx DEFAULT_MSG
    SPHINX_EXECUTABLE
)

mark_as_advanced(SPHINX_EXECUTABLE)
