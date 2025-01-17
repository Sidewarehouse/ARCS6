include_directories(
        ${CMAKE_CURRENT_LIST_DIR}
        ${CMAKE_CURRENT_LIST_DIR}/c
)

set(
        ARCS_EQUIP_files
        ${CMAKE_CURRENT_LIST_DIR}/EquipBase.cc
        ${CMAKE_CURRENT_LIST_DIR}/EquipBase.hh
        ${CMAKE_CURRENT_LIST_DIR}/EquipParams.hh
        ${CMAKE_CURRENT_LIST_DIR}/EquipTemplate.cc
        ${CMAKE_CURRENT_LIST_DIR}/EquipTemplate.hh
        ${CMAKE_CURRENT_LIST_DIR}/InterfaceFunctions.hh
        ${CMAKE_CURRENT_LIST_DIR}/c/EquipFunctions.c
        ${CMAKE_CURRENT_LIST_DIR}/c/EquipFunctions.h
)