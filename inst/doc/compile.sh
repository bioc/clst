#!/bin/bash

R CMD Sweave clstDemo.Rnw && \
texi2pdf --tidy clstDemo.tex
