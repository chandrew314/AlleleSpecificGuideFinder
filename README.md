# AlleleSpecificGuideFinder

This program is designed to return allele-specific guides to be used for CRISPR/Cas9 heterozygous gene knockouts. There are slight differences in the nucleotide sequences between two different species' genes. AlleleSpecificGuideFinder capitalizes on these differences to generate CRISPR/Cas9 sgRNAs that target one species' allele.
## Requirements

AlleleSpecificGuideFinder utililizes [CRISPOR](https://github.com/maximilianh/crisporWebsite), which is described in [this paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1012-2).

You will need to be able to run CRISPOR in order to run this program as AlleleSpecificGuideFinder uses CRISPOR-generated guides to determine whether the guides are allele specific or not.

Additionally, AlleleSpecificGuideFinder will output the target cut percentage of each guide as well as the indexes of the specific nucleotides that differ between one species and the other. Guides can be categorized as in-PAM or near-PAM based on these differences: in-PAM guides are based on differences within the PAM site (indexes 22 and 21), and near-PAM guides are based on differences near the PAM site (indexes 16-19).

AlleleSpecificGuideFinder utilizes these packages (in addition to the one CRISPOR uses):
1. Biopython
2. Biotite
3. Functools

## How to use

The script `findingsgRNAs.py` must be run in the 'crisporWebsite' directory.
```
Usage: python3.6 findingsgRNAs.py [options] species1FastaInFile species1Genome species1GuideOutFile species2FastaInFile species2Genome species2GuideOutFile guideOutputFile

positional arguments:

species1FastaInFile  = location of the FASTA file for one species

species1Genome       = genome identifier for this species

species1GuideOutFile = tab-separated file that will hold all CRISPOR-generated guides for this species

species2FastaInFile  = location of the FASTA file for another species

species2Genome       = genome identifier for this species

species2GuideOutFile = tab-separated file that will hold all CRISPOR-generated guides for this species

guideOutputFile      = tab-separated file that will contain the output of this program

species1/2Genome uses the same genome identifiers as CRISPOR.

optional arguments:
-h, --help      show this help message and exit

-PAM [PAM]      PAM motif to be used in CRISPOR; default is NGG

-m              this will only output in-PAM and near-PAM guides

-control        this will only output guides that can target both species

[-m | -control] are mutually exclusive since guides that can target both species are not in-PAM/near-PAM.
```
## Quirks

When one of the two exon sequences are much longer/shorter than the other (ie. one species has more exons than the other), the Bio.Align package used will usually align the sequences with a length of dashes. Any guides that are in this region will have an 'analogous' guide that is made entirely of dashes. Therefore, the entire guide is unique; the -m option will take these guides off.
