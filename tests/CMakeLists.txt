# Copyright András Vukics 2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

message(STATUS "Configuring tests")

# Regression: jc_qjmc against the output of the code before the switch to cppqed-multidiagonal.
# Regenerate (only for an intended change of results): compareTrajectory.py <jc_qjmc> jc_qjmc.ref --regenerate
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME jc_qjmc_regression
           COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/regression/compareTrajectory.py
                   $<TARGET_FILE:jc_qjmc> ${CMAKE_CURRENT_SOURCE_DIR}/regression/jc_qjmc.ref)
else()
  message(WARNING "Python 3 not found: regression tests disabled")
endif()
