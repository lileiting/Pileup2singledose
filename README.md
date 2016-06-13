# Pileup2singledose
A procedure to call single dose SNP from [pileup format](http://samtools.sourceforge.net/pileup.shtml)

* Author: [Lei-Ting Li](https://github.com/lileiting)
* Email: lileiting@gmail.com
* LICENCE: [BSD](http://opensource.org/licenses/bsd-license.php)


## Installation

    git clone https://github.com/lileiting/Pileup2singledose.git

## How to use it

First step, create a pileup matrix

    perl create_pileup_matrix.pl -m male.pileup -f female.pileup -p progeny1.pileup,progeny2.pileup,progeny3.pileup -o pileup.matrix.txt

Second step, call single dose

    perl pileup2singledose.pl pileup.matrix.txt > genotypes.matrix.txt

Third step, filter results

    perl filter_genotypes.pl -i genotypes.matrix.txt -t 0.2 -o genotypes.matrix.t0.2.txt

Fourth step, convert genotypes codes from h, a, b to lm, ll, nn, np

    perl convert_genotypes.pl genotypes.matrix.t0.2.txt > genotypes.matrix.t0.2.cp.txt

Fifth step, convert format for [BinMarkers](https://github.com/lileiting/BinMarkers)

    perl format4binmarkers.pl genotypes.matrix.t0.2.cp.txt > genotypes.matrix.t0.2.cp.txt.matrix

## Citation
Zhang J., Sharma A., Yu Q., Wang J., Li L., Zhu L., Zhang X., Chen Y., Ming R. Comparative structural analysis of Bru1 region homeologs in Saccharum spontaneum and S. officinarum. BMC Genomics 2016, 17:446. doi:[10.1186/s12864-016-2817-9](http://dx.doi.org/10.1186/s12864-016-2817-9); pmid:[27287040](http://www.ncbi.nlm.nih.gov/pubmed/27287040)
