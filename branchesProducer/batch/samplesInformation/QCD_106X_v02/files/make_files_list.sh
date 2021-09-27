#!/bin/bash

path=root://t3se01.psi.ch:1094//store/user/fleble/SVJ/QCD/samples/UL2018/PFNANOAODSKIM

genRange 632  ${path}/PFNANOAODSKIM_QCD_pt_170to300_part- .root > QCD_pt_170to300_files.txt
genRange 1228 ${path}/PFNANOAODSKIM_QCD_pt_300to470_part- .root > QCD_pt_300to470_files.txt
genRange 1119 ${path}/PFNANOAODSKIM_QCD_pt_470to600_part- .root > QCD_pt_470to600_files.txt
genRange 1408 ${path}/PFNANOAODSKIM_QCD_pt_600to800_part- .root > QCD_pt_600to800_files.txt
genRange 810  ${path}/PFNANOAODSKIM_QCD_pt_800to1000_part- .root > QCD_pt_800to1000_files.txt
genRange 434  ${path}/PFNANOAODSKIM_QCD_pt_1000to1400_part- .root > QCD_pt_1000to1400_files.txt
genRange 242  ${path}/PFNANOAODSKIM_QCD_pt_1400to1800_part- .root > QCD_pt_1400to1800_files.txt
genRange 136  ${path}/PFNANOAODSKIM_QCD_pt_1800to2400_part- .root > QCD_pt_1800to2400_files.txt
genRange 77   ${path}/PFNANOAODSKIM_QCD_pt_2400to3200_part- .root > QCD_pt_2400to3200_files.txt
genRange 37   ${path}/PFNANOAODSKIM_QCD_pt_3200toInf_part- .root > QCD_pt_3200toInf_files.txt
