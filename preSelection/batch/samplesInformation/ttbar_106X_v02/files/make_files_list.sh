#!/bin/bash

path=root://t3se01.psi.ch:1094//store/user/fleble/SVJ/ttbar/samples/UL2018/PFNANOAOD

genRange 3103 ${path}/PFNANOAOD_ttbar_part- .root > ttbar_files.txt
