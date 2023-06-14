cmake_minimum_required(VERSION 3.26)

#コンパイラー設定
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_C_COMPILER "gcc")

project(mc_sim)


#---設定
#compile_commands.jsonを生成する
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

#Debug/Releaseのどちらのビルドかコンパイラにわたす
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DBUILD_TYPE=\\\"Debug\\\"")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DBUILD_TYPE=\\\"Release\\\"")


#spdlog
set(spdlog_DIR "spdlog/build")
find_package(spdlog REQUIRED)


#integral
add_library(integral STATIC module/Integral.cpp)

target_include_directories(integral PUBLIC module)
target_compile_options(integral PRIVATE -Wall)
target_compile_features(integral PRIVATE cxx_std_20)


#logger_obj
add_library(logger_obj STATIC module/logger_obj.cpp)

target_include_directories(logger_obj PUBLIC module)
target_compile_options(logger_obj PRIVATE -Wall)
target_compile_options(logger_obj PRIVATE "$<$<CONFIG:Debug>:-Og>")
target_compile_options(logger_obj PRIVATE "$<$<CONFIG:Debug>:-DMC_SIM_DEBUG>")
target_compile_features(logger_obj PRIVATE cxx_std_20)

target_link_libraries(logger_obj PRIVATE spdlog::spdlog)
target_include_directories(logger_obj PUBLIC {spdlog_DIR})


#scatconst
add_library(scatconst STATIC module/scatconst.cpp)

target_include_directories(scatconst PUBLIC module)

target_compile_options(scatconst PRIVATE -Wall)
target_compile_features(scatconst PRIVATE cxx_std_20)


#curve
add_library(curve STATIC module/curve.cpp)

target_include_directories(curve PUBLIC module)
target_compile_options(curve PRIVATE -Wall)
target_compile_features(curve PRIVATE cxx_std_20)


#band_obj
add_library(band_obj STATIC module/band_obj.cpp)

target_include_directories(band_obj PUBLIC module)
target_compile_options(band_obj PRIVATE -Wall)
target_compile_features(band_obj PRIVATE cxx_std_20)

target_link_libraries(band_obj PRIVATE curve)
target_link_libraries(band_obj PRIVATE scatconst)


#physconst
add_library(physconst STATIC module/physconst.cpp)

target_include_directories(physconst PUBLIC module)
target_compile_options(physconst PRIVATE -Wall)
target_compile_features(physconst PRIVATE cxx_std_20)


#mcparticle
add_library(mcparticle STATIC module/mcparticle.cpp)

target_include_directories(mcparticle PUBLIC module)
target_compile_options(mcparticle PRIVATE -Wall)
target_compile_features(mcparticle PRIVATE cxx_std_20)

target_link_libraries(mcparticle PRIVATE physconst)


#massconst
add_library(massconst STATIC module/massconst.cpp)

target_include_directories(massconst PUBLIC module)
target_compile_options(massconst PRIVATE -Wall)
target_compile_features(massconst PRIVATE cxx_std_20)

target_link_libraries(massconst PRIVATE integral)


#simulation
add_library(simulation STATIC module/simulation.cpp)

target_include_directories(simulation PUBLIC module)
target_compile_options(simulation PRIVATE -Wall)
target_compile_options(simulation PRIVATE -fopenmp)
target_compile_options(simulation PRIVATE -MT)
target_compile_features(simulation PRIVATE cxx_std_20)

target_link_libraries(simulation PRIVATE mcparticle)




#メイン
add_executable(main test.cpp)

target_compile_options(main PRIVATE -Wall)
target_compile_features(main PRIVATE cxx_std_20)



add_executable(calculator calculator.cpp)

target_compile_options(calculator PRIVATE -Wall)
target_compile_features(calculator PRIVATE cxx_std_20)

target_link_libraries(calculator PRIVATE massconst)
target_link_libraries(calculator PRIVATE logger_obj)
target_link_libraries(calculator PRIVATE band_obj)
target_link_libraries(calculator PRIVATE curve)
target_link_libraries(calculator PRIVATE spdlog::spdlog)
target_link_libraries(calculator PRIVATE scatconst)

target_include_directories(calculator PRIVATE ${spdlog_DIR})


# target_link_libraries(main PRIVATE massconst)
# target_link_libraries(main PRIVATE simulation)
# target_link_libraries(main PRIVATE mcparticles)
# target_link_libraries(main PRIVATE physconst)
# target_link_libraries(main PRIVATE logger_obj)

