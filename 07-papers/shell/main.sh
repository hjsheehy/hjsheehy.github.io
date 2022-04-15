#!/bin/bash
rm -f main.tex
dir=$(pwd)
dir=$(basename $dir)
echo "\documentclass[11pt,oneside]{article}

\input{../tex/structure}

\begin{document}

\input{../tex/${dir}}

\end{document}" >> main.tex
