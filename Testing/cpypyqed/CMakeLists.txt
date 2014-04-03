# Copyright Raimar Sandner 2012–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

testdir(cpypyqed 1particle1mode 2particles1mode QbitMode_C++QED)

set(ONDEMAND ${PROJECT_BINARY_DIR}/ondemand_temp)

set(CPYPYQED_CONFIG ${ONDEMAND}/config.txt)
configure_file(config.txt.in ${CPYPYQED_CONFIG} @ONLY)

declaretest(cpypyqed_2particles1mode_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=PythonContinuer --script=${cpypyqed_SOURCE_DIR}/scripts/2particles1mode.py
            --cpypyqed_builddir=${ONDEMAND}/${TESTNAME} --cpypyqed_config=${CPYPYQED_CONFIG})
declaretest(cppqed_2particles1mode_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Continuer --script=$<TARGET_FILE:2particles1mode>)
declaretest(verify_2particles1mode_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=Verifier)
set_tests_properties(${TESTNAME} PROPERTIES DEPENDS "cppqed_2particles1mode_dc;cpypyqed_2particles1mode_dc")

declaretest(cpypyqed_QbitMode_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=PythonRunner --script=${cpypyqed_SOURCE_DIR}/scripts/QbitMode.py
            --cpypyqed_builddir=${ONDEMAND}/${TESTNAME} --cpypyqed_config=${CPYPYQED_CONFIG})
declaretest(cppqed_QbitMode_dc)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:QbitMode_C++QED>)
set_tests_properties(${TESTNAME} PROPERTIES DEPENDS "cpypyqed_QbitMode_dc")

# Python 2.6 does not have argparse, we have to use the version which works
# around this: 1particle1mode_without_argparse.py

if(PYTHON_VERSION_STRING VERSION_LESS 2.7)
  set(1particle1mode 1particle1mode_without_argparse.py)
else()
  set(1particle1mode 1particle1mode.py)
endif()


declaretest(cpypyqed_1particle1mode_conf1)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=PythonRunner --script=${cpypyqed_SOURCE_DIR}/scripts/${1particle1mode}
            --cpypyqed_builddir=${ONDEMAND}/${TESTNAME} --cpypyqed_config=${CPYPYQED_CONFIG})
declaretest(cppqed_1particle1mode_conf1)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:1particle1mode>)
set_tests_properties(${TESTNAME} PROPERTIES DEPENDS "cpypyqed_1particle1mode_conf1")

declaretest(cpypyqed_1particle1mode_conf2)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=PythonRunner --script=${cpypyqed_SOURCE_DIR}/scripts/${1particle1mode}
            --cpypyqed_builddir=${ONDEMAND}/${TESTNAME} --cpypyqed_config=${CPYPYQED_CONFIG})
declaretest(cppqed_1particle1mode_conf2)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:1particle1mode>)
set_tests_properties(${TESTNAME} PROPERTIES DEPENDS "cpypyqed_1particle1mode_conf2")

declaretest(cpypyqed_1particle1mode_conf3)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=PythonRunner --script=${cpypyqed_SOURCE_DIR}/scripts/${1particle1mode}
            --cpypyqed_builddir=${ONDEMAND}/${TESTNAME} --cpypyqed_config=${CPYPYQED_CONFIG})
declaretest(cppqed_1particle1mode_conf3)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:1particle1mode>)
set_tests_properties(${TESTNAME} PROPERTIES DEPENDS "cpypyqed_1particle1mode_conf3")

declaretest(cpypyqed_1particle1mode_conf4)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=PythonRunner --script=${cpypyqed_SOURCE_DIR}/scripts/${1particle1mode}
            --cpypyqed_builddir=${ONDEMAND}/${TESTNAME} --cpypyqed_config=${CPYPYQED_CONFIG})
declaretest(cppqed_1particle1mode_conf4)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=VerifiedRunner --script=$<TARGET_FILE:1particle1mode>)
set_tests_properties(${TESTNAME} PROPERTIES DEPENDS "cpypyqed_1particle1mode_conf4")
