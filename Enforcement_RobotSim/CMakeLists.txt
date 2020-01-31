cmake_minimum_required(VERSION 3.15)
project(robotsim_master C)

set(CMAKE_C_STANDARD 11)

include_directories(src)
include_directories(src/utilities)
include_directories(src/vizualizer)

add_executable(robotsim_master
        src/utilities/arenas.c
        src/utilities/arenas.h
        src/utilities/datastructs.c
        src/utilities/datastructs.h
        src/utilities/debug_utils.c
        src/utilities/debug_utils.h
        src/utilities/dynamics_utils.c
        src/utilities/dynamics_utils.h
        src/utilities/file_utils.c
        src/utilities/file_utils.h
        "src/utilities/interactions (复件).c"
        src/utilities/interactions.c
        src/utilities/interactions.h
        src/utilities/math_utils.c
        src/utilities/math_utils.h
        src/utilities/obstacles.c
        src/utilities/obstacles.h
        src/utilities/output_utils.c
        src/utilities/output_utils.h
        src/utilities/param_utils.c
        src/utilities/param_utils.h
        src/utilities/pngout_utils.c
        src/utilities/pngout_utils.h
        src/vizualizer/objects_2d.c
        src/vizualizer/objects_2d.h
        src/vizualizer/objects_3d.c
        src/vizualizer/objects_3d.h
        src/algo.h
        src/algo_gui.h
        src/algo_spp_evol.c
        src/algo_spp_evol.h
        src/algo_spp_evol_gui.c
        src/algo_spp_evol_stat.c
        src/algo_stat.h
        src/colors.c
        src/colors.h
        src/dynspecviz.c
        src/dynspecviz.h
        src/objects_menu.c
        src/objects_menu.h
        src/robotmodel.c
        src/robotmodel.h
        src/robotsim_main.c
        src/sensors.c
        src/sensors.h
        src/stat.c
        src/stat.h
        src/vizmode.h)
