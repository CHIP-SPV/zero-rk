
#N.B. CMAKE_CURRENT_BINARY_DIR needed to include these files
configure_file(zerork_version_info.h.in zerork_version_info.h)
configure_file(zerork_conf.h.in zerork_conf.h)

add_library(zerork element.cpp species.cpp mechanism.cpp utilities.cpp
            nasa_poly.cpp info_net.cpp rate_const.cpp perf_net.cpp
            fast_exps.cpp plog_reaction.cpp
            non_integer_reaction_network.cpp constants_api.cpp 
            elemental_composition.cpp impls/elemental_composition_impl.cpp)

#TODO: Cleanup include directories (e.g. "zerork/rate_const.h" vs. "rate_const.h")
target_include_directories(zerork PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../>
                                         $<INSTALL_INTERFACE:include>)

if(HAVE_ALIGNED_ALLOC)
  target_compile_definitions(zerork PRIVATE "HAVE_ALIGNED_ALLOC")
endif()

set(public_headers atomicMassDB.h constants.h constants_api.h
   element.h elemental_composition.h external_funcs.h fast_exps.h
   info_net.h mechanism.h nasa_poly.h non_integer_reaction_network.h
   perf_net.h plog_reaction.h rate_const.h species.h utilities.h)
target_link_libraries(zerork PUBLIC ckconverter)
if(NOT WIN32)
target_link_libraries(zerork PUBLIC dl m)
endif()
target_link_libraries(zerork PRIVATE zerork_vectormath)

set_target_properties(zerork PROPERTIES
 PUBLIC_HEADER  "${public_headers}")
install(TARGETS zerork
 LIBRARY DESTINATION lib
 ARCHIVE DESTINATION lib
 PUBLIC_HEADER DESTINATION include/zerork)

if(ENABLE_GPU)
add_library(zerork_cuda element.cpp species.cpp mechanism.cpp utilities.cpp
            nasa_poly.cpp info_net.cpp rate_const.cpp perf_net.cpp
            fast_exps.cpp plog_reaction.cpp
            non_integer_reaction_network.cpp constants_api.cpp
            elemental_composition.cpp impls/elemental_composition_impl.cpp
            zerork_cuda_defs.cpp nasa_poly_cuda.cpp nasa_poly_kernels.cpp
            rate_const_cuda.cpp rate_const_kernels.cpp
            perf_net_cuda.cpp perf_net_kernels.cpp
            mechanism_cuda.cpp mechanism_kernels.cpp
            scatter_add_kernels.cpp misc_kernels.cpp)

target_include_directories(zerork_cuda PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
                                              $<INSTALL_INTERFACE:include>
                                              $<BUILD_INTERFACE:/space/pvelesko/install/HIP/rocThrust/main/include>
                                       PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                                               /space/pvelesko/install/HIP/rocPRIM/2025.04.22/include)

set(public_headers_cuda atomicMassDB.h constants.h
   constants_api.h zerork_cuda_defs.h element.h
   elemental_composition.h external_funcs.h
   fast_exps.h info_net.h mechanism.h mechanism_cuda.h
   mechanism_kernels.h misc_kernels.h nasa_poly.h
   nasa_poly_cuda.h nasa_poly_kernels.h
   non_integer_reaction_network.h perf_net.h
   perf_net_cuda.h perf_net_kernels.h plog_reaction.h
   rate_const.h rate_const_cuda.h rate_const_kernels.h
   scatter_add_kernels.h species.h utilities.h)
target_link_libraries(zerork_cuda PUBLIC ckconverter)
if(NOT WIN32)
target_link_libraries(zerork_cuda PUBLIC dl m)
endif()
target_link_libraries(zerork_cuda PRIVATE zerork_vectormath)
set_target_properties(zerork_cuda PROPERTIES
 PUBLIC_HEADER  "${public_headers_cuda}")
install(TARGETS zerork_cuda
 LIBRARY DESTINATION lib
 PUBLIC_HEADER DESTINATION include/zerork)
endif()

