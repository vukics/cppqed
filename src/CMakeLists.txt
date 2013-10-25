
set(ELEMENTS_SOURCE_DIRS utils frees interactions)

foreach(d ${ELEMENTS_SOURCE_DIRS})
  add_subdirectory(${d})
endforeach(d)

add_library(C++QEDelements SHARED $<TARGET_OBJECTS:utils_objs> 
                                  $<TARGET_OBJECTS:frees_objs>
                                  $<TARGET_OBJECTS:interactions_objs>)
target_link_libraries(C++QEDelements ${CPPQED_LIBRARIES})