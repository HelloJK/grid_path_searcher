cmake_minimum_required(VERSION 2.8.3)
project(grid_path_searcher)

find_package(catkin REQUIRED COMPONENTS 
    roscpp 
    std_msgs
    nav_msgs
    visualization_msgs 
)

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
#find_package(OMPL REQUIRED)

set(Eigen3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

catkin_package(
  INCLUDE_DIRS include
)

include_directories(
    include 
    include/ompl_rrt/
    SYSTEM 
    third_party
    ${catkin_INCLUDE_DIRS} 
    ${Eigen3_INCLUDE_DIRS} 
    ${PCL_INCLUDE_DIRS}
    # ${OMPL_INCLUDE_DIRS}
)

link_libraries(${OMPL_LIBRARY_DIRS})

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} -O3 -Wall") # -Wextra -Werror

add_executable( demo_node 
    src/demo_node.cpp
    src/Astar_searcher.cpp
    src/jps_search/JPS_utils.cpp
    src/jps_search/JPS_searcher.cpp
    )

target_link_libraries(demo_node 
    ${catkin_LIBRARIES}
    ${PCL_LIBRARIES} 
)

add_executable ( random_complex 
    src/random_complex_generator.cpp )

target_link_libraries( random_complex
    ${catkin_LIBRARIES}
    ${PCL_LIBRARIES} )  

#add_executable (ompl_rrt
#    src/ompl_rrt_node.cpp
#    src/ompl_rrt/ompl_rrt_searcher.cpp)

#target_link_libraries(ompl_rrt
#   ${catkin_LIBRARIES}
#   ${PCL_LIBRARIES}
#   ${OMPL_LIBRARIES} )