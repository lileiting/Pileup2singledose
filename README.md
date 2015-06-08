# Pileup2singledose
A procedure to call single dose SNP from [pileup format](http://samtools.sourceforge.net/pileup.shtml)

Installation
------

    git clone https://github.com/lileiting/Pileup2singledose.git

How to use it
------

First step, create a pileup matrix

    perl createpileupmatrix.pl -m male.pileup -f female.pileup -p progeny1.pileup,progeny2.pileup,progeny3.pileup -o pileup.matrix.txt

Second step, call single dose

    perl pileup2singledose.pl pileup.matrix.txt > genotypes.matrix.txt

Third step, filter results

    perl filtergenotypes.pl -i genotypes.matrix.txt -t 0.2 -o genotypes.matrix.t0.2.txt

