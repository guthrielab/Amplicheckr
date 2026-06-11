## To use: 
### Command line:
```
python3 main.py -v hantavirus .\testdata_hantavirus\tprimerhanta.csv .\testdata_hantavirus\tgenome.fna
```
### HTML:
```
python3 main.py -v hantavirus .\testdata_hantavirus\tprimerhanta.csv .\testdata_hantavirus\tgenome.fna -H
```
#### Swap " \\ " with " / " if you are on MacOS or linux systems.

## Test data info 

[Primer sets #1 and #2](https://www.hug.ch/sites/interhug/files/structures/centre_maladies_virales_emergentes/Documents/andes-virus-rt-pcr-protocol-hug_080526.pdf)
#### This is a case where since a "set" has two probes for only 1 forward and reverse primer, they are split into two sets with identical fwd and rev but different probes.

[Primer set #3](https://pmc.ncbi.nlm.nih.gov/articles/PMC2704762/)

Genomes are the 27 most recent entries on NCBI for andes hantavirus. 

Note that entries [OQ09224(3-5)](https://www.ncbi.nlm.nih.gov/nuccore/OQ092243) are from a different strain known as VARS/22-01 and are intentionally included to show an example of what the warning messages look like.
