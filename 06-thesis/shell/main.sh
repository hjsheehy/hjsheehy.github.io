#!/bin/bash
rm -f main.tex
dir=$(pwd)
dir=$(basename $dir)
echo "\input{../tex/structure}

\begin{document}

\input{../tex/${dir}}

\end{document}" >> main.tex
