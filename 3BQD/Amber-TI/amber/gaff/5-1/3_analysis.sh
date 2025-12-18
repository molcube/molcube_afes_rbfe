#!/bin/bash

top=${PWD}

mbar_file_path=$(find $(pwd) -name "mbar.dat" -print -quit)
list_string=$(paste -sd, "$mbar_file_path" | sed 's/^/[/' | sed 's/$/]/')

# Check if on macOS or Linux
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS sed command
    sed -i '' "s/KEYWORD/$list_string/g" ${top}/scripts/stdti_step2dats.py
else
    # Linux sed command
    sed -i "s/KEYWORD/$list_string/g" ${top}/scripts/stdti_step2dats.py
fi

trail=('t01' 't02' 't03' 't04')
for ri in ligand complex; do
    for t in ${trail[@]}; do
        cd ${ri}/${t}
        echo "${ri} ${t}"
        if [ ! -d results ]; then
             mkdir results
        fi
        for l in $(cat mbar.dat); do
            cd ${l}
            python3 ${top}/scripts/stdti_step2dats.py *prod.out

            cp *.dat ../results
            cd ..
        done
        cd ../..
    done
done


python3 scripts/MakeLaTeX_unified.py

echo 'Till this point, the necessary data has been collected.'
echo 'Below is just to use LaTeX to generate PDF file. '
echo 'We can also download the latex/ folder and generate PDF in local machine.'

cd latex
latex Main.tex
latex Main.tex

if [ -e Main.dvi ]; then
   dvips -j0 -Ppdf -G0 -tletter Main.dvi
fi
if [ -e Main.ps ]; then
  ps2pdf -dCompatibilityLevel=1.3 -dMAxSubsetPct=100 -dSubsetFonts=true -dEmbedAllFonts=true -sPAPERSIZE=letter -dEPSCrop Main.ps
fi

exit
