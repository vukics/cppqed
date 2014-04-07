# Copyright Raimar Sandner 2012–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

testdir(compile)

CPPQED_SETUP()

add_library(TMP_ToolsCompile OBJECT EXCLUDE_FROM_ALL TMP_ToolsCompile.cc)
add_library(BlitzArraySliceIteratorTMP_Compile OBJECT EXCLUDE_FROM_ALL BlitzArraySliceIteratorTMP_Compile.cc)

declaretest(test_CustomElementsExample)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=CompileTarget --script=C++QEDelements_custom_example-${CPPQED_ID})
declaretest(test_CustomScriptsExample)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=CompileTarget --script=scripts_custom_example_all)
set_tests_properties(${TESTNAME} PROPERTIES DEPENDS test_CustomElementsExample)
declaretest(test_scripts_all_target)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=CompileTarget --script=scripts_all)
declaretest(test_core_examples_main_targets)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=CompileTarget --script=core_examples_main_targets)
declaretest(test_TMP_ToolsCompile)
add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
            --test=${TESTNAME} --testclass=CompileTarget --script=TMP_ToolsCompile)
#declaretest(test_BlitzArraySliceIteratorTMP_Compile)
#add_test( NAME ${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} ${TESTSCRIPT}
#            --test=${TESTNAME} --testclass=CompileTarget --script=BlitzArraySliceIteratorTMP_Compile)