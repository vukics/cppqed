include(PythonSetup)
include(${core_BINARY_DIR}/CPPQEDConfig.cmake)
CPPQED_SETUP()

include_directories(${PYTHON_INCLUDE_DIRS})

SET(MANDATORY_LIBRARIES ${CPPQEDLIB} ${Boost_PYTHON_LIBRARY})

set(CPYPYQED_INSTALL_PATH ${PYTHON_MODULES_PATH}/cpypyqed)

set(PYTHON_MAX_RANK 5 CACHE STRING "Pre-instantiate template classes up to this rank in Python bindings.")
set_property(CACHE PYTHON_MAX_RANK PROPERTY STRINGS 1 2 3 4 5 6 7 8 9)

# Macro used for each module to select the correct extension according to the OS
MACRO(PYTHON_CPPMODULE arg1)
  SET_TARGET_PROPERTIES(_${arg1} PROPERTIES PREFIX "")
  IF(PYDEXTENSION)
    SET_TARGET_PROPERTIES(_${arg1} PROPERTIES SUFFIX ".pyd")
  ELSEIF(SOEXTENSION)
    SET_TARGET_PROPERTIES(_${arg1} PROPERTIES SUFFIX ".so")
  ENDIF(PYDEXTENSION)
  TARGET_LINK_LIBRARIES(_${arg1} ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARIES})
  #STRING(REGEX REPLACE "^_" "" DIR_ARG ${arg1})
  INSTALL(TARGETS _${arg1}
                RUNTIME DESTINATION ${CPYPYQED_INSTALL_PATH}/${arg1}
                LIBRARY DESTINATION ${CPYPYQED_INSTALL_PATH}/${arg1}
  )
ENDMACRO(PYTHON_CPPMODULE)

# Macro used for sub-package creation
MACRO(PYTHON_SUBPACKAGE PACKAGE)
  SET(PACKAGENAME _${PACKAGE})
  AUX_SOURCE_DIRECTORY(. PY_SOURCE_FILES)

  foreach(s ${PY_SOURCE_FILES})
    file( STRINGS ${s} EXPORT_HELPER_LINE REGEX "[ ]*void[ ]+export_" )
    #get_filename_component(ITEM ${s} NAME_WE)
    #set(SEQUENCE ${SEQUENCE}\(${ITEM}\))
    string(REGEX REPLACE ".*(export_.*)\\(.*" "\\1" EXPORT_HELPER "${EXPORT_HELPER_LINE}"  )
    set(EXPORT_HELPERS ${EXPORT_HELPERS} ${EXPORT_HELPER})
  endforeach(s)
  if(EXPORT_HELPERS)
    list(SORT EXPORT_HELPERS)
  endif(EXPORT_HELPERS)
  foreach(s ${EXPORT_HELPERS})
    set(EXPORT_HELPER_SEQUENCE ${EXPORT_HELPER_SEQUENCE}\(${s}\))
  endforeach(s)
  

  configure_file(../details/config.h.in config.h @ONLY)
  configure_file(../details/exportfunctions.h.in exportfunctions.h @ONLY)
  include_directories(${CMAKE_CURRENT_BINARY_DIR})

  configure_file(../details/build_module.cc.in build${PACKAGENAME}.cc @ONLY)

  add_library(${PACKAGENAME} SHARED build${PACKAGENAME}.cc ${PY_SOURCE_FILES} ${ARGN})
  PYTHON_CPPMODULE(${PACKAGE})

  configure_file(../details/__init__.py.in __init__.py @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/__init__.py DESTINATION ${CPYPYQED_INSTALL_PATH}/${PACKAGE})
ENDMACRO(PYTHON_SUBPACKAGE)

include_directories(include)

add_custom_target(cpypyqed)

# build python io module
if(ENABLE_PYTHON_IO)
  add_subdirectory(io)
  add_dependencies(cpypyqed _io)
endif()

foreach(file __init__.py)
  configure_file(${file} ${file} COPYONLY)
  install(FILES ${file} DESTINATION ${CPYPYQED_INSTALL_PATH})
endforeach(file)