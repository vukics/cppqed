find_package(GSL REQUIRED)

include_directories(include ${GSL_INCLUDE_DIRS})

AUX_SOURCE_DIRECTORY(src UTILS_SOURCE_FILES)

# checks if libutils should be linked with FLENS
IF( FLENS )
  set(UTILS_OPTLIBS ${UTILS_OPTLIBS} ${flens_LIBRARIES})
ELSE( FLENS )
  # linkage of this will fail without FLENS support
  LIST(REMOVE_ITEM UTILS_SOURCE_FILES src/DrivenDampedHarmonicOscillator.cc)
ENDIF( FLENS )

set(UTILS_LIBS ${UTILS_OPTLIBS} ${GSL_LIBRARIES} ${blitz++_LIBRARIES})

add_library(utils SHARED ${UTILS_SOURCE_FILES})
target_link_libraries(utils ${UTILS_LIBS})
set_target_properties( utils PROPERTIES INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib )

install(TARGETS utils
        EXPORT C++QEDLibraryDepends
	LIBRARY DESTINATION ${INSTALL_LIB_DIR}
	COMPONENT shlib)
