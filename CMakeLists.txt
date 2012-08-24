
set(ELEMENTS_SOURCE_DIRS utils frees interactions composites)

foreach(d ${ELEMENTS_SOURCE_DIRS})
  aux_source_directory(${d} ELEMENTS_SOURCE_FILES)
  set(ELEMENTS_INCLUDE_DIRS ${ELEMENTS_INCLUDE_DIRS} elements/${d})
endforeach(d)
include_directories(${ELEMENTS_SOURCE_DIRS})
# define ELEMENTS_INCLUDE_DIRS in parent scope, so that it is accessible to scripts
set(ELEMENTS_INCLUDE_DIRS ${ELEMENTS_INCLUDE_DIRS} PARENT_SCOPE)

add_library(el SHARED ${ELEMENTS_SOURCE_FILES})
target_link_libraries(el C++QED)
set_target_properties( el PROPERTIES INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib )

install(TARGETS el
        EXPORT C++QEDLibraryDepends
	LIBRARY DESTINATION ${INSTALL_LIB_DIR}
	COMPONENT shlib)
