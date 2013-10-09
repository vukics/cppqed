
set(ELEMENTS_SOURCE_DIRS utils frees interactions)

foreach(d ${ELEMENTS_SOURCE_DIRS})
  #aux_source_directory(${d} ELEMENTS_SOURCE_FILES)
  file(GLOB SOURCE_FILES RELATIVE ${CMAKE_SOURCE_DIR} ${d}/*.cc)
  set(ELEMENTS_SOURCE_FILES ${ELEMENTS_SOURCE_FILES} ${SOURCE_FILES})
  set(ELEMENTS_INCLUDE_DIRS ${ELEMENTS_INCLUDE_DIRS} elements/${d})
endforeach(d)
include_directories(${ELEMENTS_SOURCE_DIRS})
# define ELEMENTS_INCLUDE_DIRS in parent scope, so that it is accessible to scripts
set(ELEMENTS_INCLUDE_DIRS ${ELEMENTS_INCLUDE_DIRS} PARENT_SCOPE)
set(ELEMENTS_SOURCE_FILES ${ELEMENTS_SOURCE_FILES} PARENT_SCOPE)

#add_library(C++QEDel SHARED ${ELEMENTS_SOURCE_FILES})
#target_link_libraries(C++QEDel C++QEDcore)
#set_target_properties(C++QEDel
#		      PROPERTIES INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib
#                      VERSION ${CPPQED_VERSION}
#                      SOVERSION 1)

#install(TARGETS C++QEDel
#        EXPORT CPPQEDLibraryDepends
#	LIBRARY DESTINATION ${INSTALL_LIB_DIR}
#	COMPONENT shlib)
