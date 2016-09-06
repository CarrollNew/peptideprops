# peptideprops
ExPASy-like protein/mRNA properties calculator for company internal needs.

*requires RNAfold*

wget https://www.tbi.univie.ac.at/RNA/download/arch/arch_extra/viennarna-2.2.8-1-x86_64.pkg.tar.xz (Arch Linux)

tar -xf viennarna-2.2.8-1-x86_64.pkg.tar.xz

sudo mkdir /usr/bin/ViennaRNA

sudo cp -r usr/bin/* /usr/bin/ViennaRNA

*usage example*

python main.py sample_data/sample_input.json sample_data/sample_output.json

Testing requires BioPython and ViennaRNA installed.

*possible test launch commands:*

python -m unittest discover . bio_tests.py

python -m unittest discover . plotter.py
