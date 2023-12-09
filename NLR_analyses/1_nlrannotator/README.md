## Annotating the NLR genes in the maize NAM lines

An automated NLR gene annotation pipeline was used for detecting NLR motifs as follows.

`java -jar NLR-Annotator-v2.1b.jar -t 4 -i B73.PLATINUM.pseudomolecules-v1.fasta.mod -x mot.txt -y store.txt -o B73.PLATINUM.pseudomolecules-v1.out -b B73.PLATINUM.pseudomolecules-v1.bed`

The full list of commands is available in `scripts/run_nlrannotator.sh`. The resulting motifs were intersected with the NAM gene annotations and all intersecting genes were considered NLR genes. The results are provided for all NAM lines, for example in `data/B73.nlr.bed`. An independently annotated set of NLR genes in [Maize_NLRome_GeneTable.txt](https://github.com/daniilprigozhin/NLRCladeFinder/blob/main/Maize_NLRome/Maize_NLRome_GeneTable.txt) was obtained from the repository [https://github.com/daniilprigozhin/NLRCladeFinder/](https://github.com/daniilprigozhin/NLRCladeFinder/).
