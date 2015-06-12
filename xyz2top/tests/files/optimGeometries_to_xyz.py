#!/usr/bin/env python

import read_LSDALTON_output as readLS

filenames = ["valinomycin_geomOpt_DFT-b3lyp_6-31Gs_ADMM2-KT3X_3-21G.out",
#             "valinomycin_geomOpt_DFT-b3lyp_6-31Gs_ADMM2-PBEX_3-21G.out",
#             "valinomycin_geomOpt_DFT-b3lyp-noDF_6-31Gs.out",
#             "valinomycin_geomOpt_DFT-b3lyp_6-31Gs.out",
#             "valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ_ADMM2-KT3X_3-21G.out",
#             "valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ_ADMM2-PBEX_3-21G.out",
             "valinomycin_geomOpt_DFT-b3lyp-noDF_cc-pVTZ.out",
             "valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ.out"]

path = "/home/ctcc2/Documents/CODE-DEV/xyz2top/xyz2top/tests/files/"
filename = "valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ.out"

for filename in filenames:
    molOptimized = readLS.parse_molecule_optimized(path + filename)
 #   print molOptimized
    xyz = molOptimized.getContent_format_XYZ()
    print xyz
    #    print xyz
    with open(filename+".xyz", 'w') as outfile:
        outfile.write(xyz)

# valinomycin_geomOpt_DFT-b3lyp_6-31Gs_ADMM2-KT3X_3-21G.out
# valinomycin_geomOpt_DFT-b3lyp_6-31Gs_ADMM2-PBEX_3-21G.out
# valinomycin_geomOpt_DFT-b3lyp_6-31Gs_noDF.out
# valinomycin_geomOpt_DFT-b3lyp_6-31Gs.out
# valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ_ADMM2-KT3X_3-21G.out
# valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ_ADMM2-PBEX_3-21G.out
# valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ_noDF.out
# valinomycin_geomOpt_DFT-b3lyp_cc-pVTZ.out
