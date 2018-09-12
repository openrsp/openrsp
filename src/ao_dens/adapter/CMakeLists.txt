target_sources(
  openrsp
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/openrsp_callback_f.F90
    ${CMAKE_CURRENT_LIST_DIR}/OpenRSPGetRSPFun_f.F90
    ${CMAKE_CURRENT_LIST_DIR}/OpenRSPGetResidue_f.F90
  )
