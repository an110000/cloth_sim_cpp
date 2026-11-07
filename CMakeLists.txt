cmake_minimum_required(VERSION 3.21)
project(cloth_sim_cpp VERSION 0.6.0 LANGUAGES CXX)

option(CLOTH_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
include(cmake/ClothEigen.cmake)

add_library(cloth_core src/mesh.cpp src/material.cpp src/sparse.cpp
    src/mass.cpp src/assembly.cpp src/dofs.cpp src/linear_solver.cpp src/implicit.cpp
    src/constraints.cpp src/xpbd.cpp src/contact.cpp src/scene.cpp src/experiment.cpp)
target_link_libraries(cloth_core PRIVATE cloth_eigen_headers)
target_include_directories(cloth_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_features(cloth_core PUBLIC cxx_std_17)
set_target_properties(cloth_core PROPERTIES CXX_EXTENSIONS OFF)

function(cloth_compile_options target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive- /utf-8)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    if(CLOTH_ENABLE_SANITIZERS)
        if(WIN32 OR NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            message(FATAL_ERROR "This sanitizer preset requires GCC or Clang on Linux/macOS")
        endif()
        target_compile_options(${target} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE -fsanitize=address,undefined)
    endif()
endfunction()
cloth_compile_options(cloth_core)

add_executable(cloth_inspect apps/cloth_inspect.cpp)
target_link_libraries(cloth_inspect PRIVATE cloth_core)
cloth_compile_options(cloth_inspect)

add_executable(cloth_assembly_inspect apps/cloth_assembly_inspect.cpp)
target_link_libraries(cloth_assembly_inspect PRIVATE cloth_core)
cloth_compile_options(cloth_assembly_inspect)

add_executable(cloth_simulate apps/cloth_simulate.cpp)
target_link_libraries(cloth_simulate PRIVATE cloth_core)
cloth_compile_options(cloth_simulate)

add_executable(cloth_xpbd_simulate apps/cloth_xpbd_simulate.cpp)
target_link_libraries(cloth_xpbd_simulate PRIVATE cloth_core)
cloth_compile_options(cloth_xpbd_simulate)

add_executable(cloth_scene_run apps/cloth_scene_run.cpp)
target_link_libraries(cloth_scene_run PRIVATE cloth_core)
cloth_compile_options(cloth_scene_run)

include(CTest)
if(BUILD_TESTING)
    add_executable(test_mesh tests/test_mesh.cpp)
    target_link_libraries(test_mesh PRIVATE cloth_core)
    cloth_compile_options(test_mesh)
    add_test(NAME mesh_and_mass COMMAND test_mesh)

    add_executable(test_material tests/test_material.cpp)
    target_link_libraries(test_material PRIVATE cloth_core)
    cloth_compile_options(test_material)
    add_test(NAME material_properties COMMAND test_material)

    add_executable(test_python_reference tests/test_python_reference.cpp)
    target_link_libraries(test_python_reference PRIVATE cloth_core)
    cloth_compile_options(test_python_reference)
    add_test(NAME python_reference COMMAND test_python_reference
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/python_reference.txt")
    set_tests_properties(mesh_and_mass material_properties python_reference PROPERTIES TIMEOUT 30)

    foreach(test_target IN ITEMS test_hessian test_sparse_mass test_assembly test_linear_solver test_implicit
            test_constraints test_xpbd test_energy_xpbd test_contact test_implicit_contact test_xpbd_contact
            test_scene test_experiment)
        add_executable(${test_target} tests/${test_target}.cpp)
        target_link_libraries(${test_target} PRIVATE cloth_core)
        cloth_compile_options(${test_target})
        add_test(NAME ${test_target} COMMAND ${test_target})
        set_tests_properties(${test_target} PROPERTIES TIMEOUT 30)
    endforeach()
    add_executable(test_task02_reference tests/test_task02_reference.cpp)
    target_link_libraries(test_task02_reference PRIVATE cloth_core)
    cloth_compile_options(test_task02_reference)
    add_test(NAME task02_python_reference COMMAND test_task02_reference
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/task02_reference.txt")
    set_tests_properties(task02_python_reference PROPERTIES TIMEOUT 30)
    add_executable(test_task03_reference tests/test_task03_reference.cpp)
    target_link_libraries(test_task03_reference PRIVATE cloth_core)
    cloth_compile_options(test_task03_reference)
    add_test(NAME task03_python_trajectories COMMAND test_task03_reference
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/task03_reference.txt")
    set_tests_properties(task03_python_trajectories PROPERTIES TIMEOUT 30)
    add_test(NAME cli_success_and_rollback COMMAND "${CMAKE_COMMAND}"
        "-DSIMULATOR=$<TARGET_FILE:cloth_simulate>"
        "-DOUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR}/cli-test"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cli.cmake")
    set_tests_properties(cli_success_and_rollback PROPERTIES TIMEOUT 60)
    add_executable(test_task04_reference tests/test_task04_reference.cpp)
    target_link_libraries(test_task04_reference PRIVATE cloth_core)
    cloth_compile_options(test_task04_reference)
    add_test(NAME task04_python_reference COMMAND test_task04_reference
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/task04_reference.txt")
    set_tests_properties(task04_python_reference PROPERTIES TIMEOUT 30)
    add_test(NAME xpbd_cli_budget_and_strict COMMAND "${CMAKE_COMMAND}"
        "-DSIMULATOR=$<TARGET_FILE:cloth_xpbd_simulate>"
        "-DOUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR}/xpbd-cli-test"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_xpbd_cli.cmake")
    set_tests_properties(xpbd_cli_budget_and_strict PROPERTIES TIMEOUT 60)
    add_executable(test_task05a_reference tests/test_task05a_reference.cpp)
    target_link_libraries(test_task05a_reference PRIVATE cloth_core)
    cloth_compile_options(test_task05a_reference)
    add_test(NAME task05a_python_reference COMMAND test_task05a_reference
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/task05a_reference.txt")
    set_tests_properties(task05a_python_reference PROPERTIES TIMEOUT 30)
    add_test(NAME energy_xpbd_cli_contract COMMAND "${CMAKE_COMMAND}"
        "-DSIMULATOR=$<TARGET_FILE:cloth_xpbd_simulate>"
        "-DOUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR}/energy-xpbd-cli-test"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_energy_xpbd_cli.cmake")
    set_tests_properties(energy_xpbd_cli_contract PROPERTIES TIMEOUT 60)
    add_executable(test_task05b_reference tests/test_task05b_reference.cpp)
    target_link_libraries(test_task05b_reference PRIVATE cloth_core)
    cloth_compile_options(test_task05b_reference)
    add_test(NAME task05b_python_reference COMMAND test_task05b_reference
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/task05b_reference.txt")
    set_tests_properties(task05b_python_reference PROPERTIES TIMEOUT 30)
    add_test(NAME contact_cli_contract COMMAND "${CMAKE_COMMAND}"
        "-DIMPLICIT=$<TARGET_FILE:cloth_simulate>" "-DXPBD=$<TARGET_FILE:cloth_xpbd_simulate>"
        "-DOUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR}/contact-cli-test"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_contact_cli.cmake")
    set_tests_properties(contact_cli_contract PROPERTIES TIMEOUT 60)
    add_test(NAME task06_scene_cli_contract COMMAND "${CMAKE_COMMAND}"
        "-DSIMULATOR=$<TARGET_FILE:cloth_scene_run>"
        "-DSCENE=${CMAKE_CURRENT_SOURCE_DIR}/tests/scenes/task06_teaching.scene"
        "-DOUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR}/task06-scene-cli-test"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_scene_cli.cmake")
    set_tests_properties(task06_scene_cli_contract PROPERTIES TIMEOUT 60)
    add_test(NAME task06_implicit_scene_cli_contract COMMAND "${CMAKE_COMMAND}"
        "-DSIMULATOR=$<TARGET_FILE:cloth_scene_run>"
        "-DSCENE=${CMAKE_CURRENT_SOURCE_DIR}/tests/scenes/task06_implicit.scene"
        "-DOUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR}/task06-implicit-scene-cli-test"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_scene_cli.cmake")
    set_tests_properties(task06_implicit_scene_cli_contract PROPERTIES TIMEOUT 60)
    add_test(NAME task06_energy_scene_cli_contract COMMAND "${CMAKE_COMMAND}"
        "-DSIMULATOR=$<TARGET_FILE:cloth_scene_run>"
        "-DSCENE=${CMAKE_CURRENT_SOURCE_DIR}/tests/scenes/task06_energy.scene"
        "-DOUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR}/task06-energy-scene-cli-test"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_scene_cli.cmake")
    set_tests_properties(task06_energy_scene_cli_contract PROPERTIES TIMEOUT 60)
    add_test(NAME task07_custom_mesh_scene_cli_contract COMMAND "${CMAKE_COMMAND}"
        "-DSIMULATOR=$<TARGET_FILE:cloth_scene_run>"
        "-DSCENE=${CMAKE_CURRENT_SOURCE_DIR}/tests/scenes/task07_custom_mesh.scene"
        "-DOUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR}/task07-custom-mesh-scene-cli-test"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_scene_cli.cmake")
    set_tests_properties(task07_custom_mesh_scene_cli_contract PROPERTIES TIMEOUT 60)
endif()
