#!/usr/bin/bash
lualatex --interaction=nonstopmode 0extract-links.tex
rm 0extract-links.aux
rm 0extract-links.log
