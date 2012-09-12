find_package(GSL REQUIRED)

include_directories(include ${GSL_INCLUDE_DIRS})

#AUX_SOURCE_DIRECTORY(src UTILS_SOURCE_FILES)
file(GLOB UTILS_SOURCE_FILES RELATIVE .. src/*.cc)

# checks if libutils should be linked with FLENS
IF( FLENS )
  set(UTILS_OPTLIBS ${UTILS_OPTLIBS} ${flens_LIBRARIES})
ELSE( FLENS )
  # linkage of this will fail without FLENS support
  LIST(REMOVE_ITEM UTILS_SOURCE_FILES utils/src/DrivenDampedHarmonicOscillator.cc)
ENDIF( FLENS )

set(UTILS_LIBS ${GSL_LIBRARIES} ${blitz_LIBRARIES})
set(UTILS_LIBS ${UTILS_LIBS} PARENT_SCOPE)
set(UTILS_SOURCE_FILES ${UTILS_SOURCE_FILES} PARENT_SCOPE)

#add_library(C++QEDutils SHARED ${UTILS_SOURCE_FILES})
#target_link_libraries(C++QEDutils ${UTILS_LIBS})
#set_target_properties(C++QEDutils
#		      PROPERTIES
#		      INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib
#                      VERSION ${CPPQED_VERSION}
#                      SOVERSION 1)
#install(TARGETS C++QEDutils
#        EXPORT CPPQEDLibraryDepends
#	LIBRARY DESTINATION ${INSTALL_LIB_DIR}
#	COMPONENT shlib)
