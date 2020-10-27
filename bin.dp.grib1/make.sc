#!/bin/bash

make clean

rm -f rams_wgrib dgrib-*

gcc -o rams_wgrib rams_wgrib.c

make -f Makefile

