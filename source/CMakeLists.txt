target_sources(Inversion PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/FileIO.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/GravimetryInversion.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/MeasurementData.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Norms.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Result.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/utils.cpp"
        )

target_sources(Interpolation PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/FileIO.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/GravimetryInversion.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/MeasurementData.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Norms.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Result.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/utils.cpp"
        )

target_sources(Inversion PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/FileIO.h"
        "${CMAKE_CURRENT_LIST_DIR}/GravimetryInversion.h"
        "${CMAKE_CURRENT_LIST_DIR}/MeasurementData.h"
        "${CMAKE_CURRENT_LIST_DIR}/Norms.h"
        "${CMAKE_CURRENT_LIST_DIR}/Result.h"
        "${CMAKE_CURRENT_LIST_DIR}/utils.h"
        )
target_sources(Interpolation PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/FileIO.h"
        "${CMAKE_CURRENT_LIST_DIR}/GravimetryInversion.h"
        "${CMAKE_CURRENT_LIST_DIR}/MeasurementData.h"
        "${CMAKE_CURRENT_LIST_DIR}/Norms.h"
        "${CMAKE_CURRENT_LIST_DIR}/Result.h"
        "${CMAKE_CURRENT_LIST_DIR}/utils.h"
        )