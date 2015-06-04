#!/usr/bin/env python

import read_LSDALTON_output as readLS

path = "/home/ctcc2/Documents/CODE-DEV/xyz2top/xyz2top/tests/files"
filename = "valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ.out"
str = readLS.get_optmized_MOL_string(path + filename)

# valinomycin_geomOpt_DFT-b3lyp_6-31Gs_ADMM2-KT3X_3-21G.out
# valinomycin_geomOpt_DFT-b3lyp_6-31Gs_ADMM2-PBEX_3-21G.out
# valinomycin_geomOpt_DFT-b3lyp_6-31Gs_noDF.out
# valinomycin_geomOpt_DFT-b3lyp_6-31Gs.out
# valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ_ADMM2-KT3X_3-21G.out
# valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ_ADMM2-PBEX_3-21G.out
# valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ_noDF.out
# valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ.out
