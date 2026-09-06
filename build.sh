#!/bin/bash
cmake --fresh -B build -G "Visual Studio 17 2022" -A x64 -DCMAKE_INSTALL_PREFIX=%pwd%/2DFEM_INSTALL
cmake --build build --config Release 
cmake --install build --config Release