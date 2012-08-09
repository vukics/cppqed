
set(ELEMENTS_SOURCE_DIRS utils frees interactions composites)

foreach(d ${ELEMENTS_SOURCE_DIRS})
  aux_source_directory(${d} ELEMENTS_SOURCE_FILES)
  set(ELEMENTS_INCLUDE_DIRS ${ELEMENTS_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/${d})
endforeach(d)
include_directories(${ELEMENTS_SOURCE_DIRS})
# define ELEMENTS_INCLUDE_DIRS in parent scope, so that it is accessible to scripts
set(ELEMENTS_INCLUDE_DIRS ${ELEMENTS_INCLUDE_DIRS} PARENT_SCOPE)

add_library(el SHARED ${ELEMENTS_SOURCE_FILES})
target_link_libraries(el C++QED)
