cmake_minimum_required(VERSION 3.26)

project(cmaketest)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER gcc)


#---設定
#compile_commands.jsonを生成する
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

#warning表示
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

#c++17を使用
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

#OpenMP並列化を使用
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

#Debug/Releaseのどちらのビルドかコンパイラにわたす
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DBUILD_TYPE=\\\"Debug\\\"")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DBUILD_TYPE=\\\"Release\\\"")


#---コンパイル
#メイン
add_executable(main test.cpp)

#spdlog
set(spdlog_DIR "spdlog/build")
find_package(spdlog REQUIRED)
include_directories(${spdlog_SOURCE_DIR})

include_directories(module)

add_library(simulation STATIC module/Simulation.cpp)
add_library(mcparticles STATIC module/mcparticles.cpp)
add_library(physconst STATIC module/physconst.cpp)
add_library(massconst STATIC module/massconst.cpp)
add_library(integral STATIC module/Integral.cpp)

target_link_libraries(main simulation)
target_link_libraries(main mcparticles)
target_link_libraries(main physconst)
target_link_libraries(main massconst)
target_link_libraries(main integral)
target_link_libraries(main spdlog::spdlog)
