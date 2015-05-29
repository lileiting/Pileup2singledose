# Pileup2singledose
A procedure to call single dose SNP

Installation
------

    git clone https://github.com/lileiting/PileupToSNP.git

How to use it
------

First step, create a pileup matrix

    perl get_matrix -m male.pileup -f female.pileup -p progeny1.pileup,progeny2.pileup,progeny3.pileup -o pileup.matrix.txt

Second step, call single dose

    perl pileuptosnp.pl pileup.matrix.txt
