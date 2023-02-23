# crisprHAL: A generalizable Cas9/sgRNA prediction model using machine transfer learning with small high-quality datasets

The CRISPR/Cas9 nuclease from Streptococcus pyogenes SpCas9 can be used with single guide RNAs (sgRNAs) as a sequence-specific antimicrobial agent and as a genome-engineering tool. However, current bacterial sgRNA activity models poorly predict SpCas9/sgRNA activity and are not generalizable, possibly because the underlying datasets used to train the models do not accurately measure SpCas9/sgRNA cleavage activity and cannot distinguish cleavage activity from toxicity. We solved this problem by using a two-plasmid positive selection system to generate high-quality biologically-relevant data that more accurately reports on SpCas9/sgRNA cleavage activity and that separates activity from toxicity. We developed a new machine transfer learning architecture (crisprHAL) that can be trained on existing datasets and that shows marked improvements in sgRNA activity prediction accuracy when transfer learning is used with small amounts of high-quality data. The crisprHAL model recapitulates known SpCas9/sgRNA-target DNA interactions and provides a pathway to a generalizable sgRNA bacterial activity prediction tool.

```
python crisprHAL.py [Enzyme] [Input File csv] [Optional Compare]
```

Enzyme: "TevSpCas9" or "SpCas9"

Input: CSV file input name (format shown below)

Optional Compare: "Compare" if your CSV file contains activity scores for comparison (format shown below)

```
sgRNA Input Required:
• 28 nucleotides long
• 20 nucleotide target site, 3 nucleotide NGG PAM, 5 nucleotides downstream
• CTCGATTGAGGGGCTGGGAA TGG GTGAT
```

Example Input CSV File (No "Compare" Option):
```

```

Example Input CSV File with "Compare" Option:
```
sgRNA,score
GGCTTTAGACAATGCAGTATTGGCGGTC,26.4036360965
CTCGATTGAGGGGCTGGGAATGGGTGAT,8.21062788839
ATCTTTATCGGTCTTAGCAAAGGCTTTG,30.1092205446
CGGGCCAGACTGGCCGAGACGGGTCGTT,11.0586722001
CAGCATCATGCCGTGATCGCCGGGAAAG,0.67958058668
GGGGAAAGGACAAGCAGTGCGGGCAAAA,29.2338752707
GAGAGAATTTTGACCTGATCCGGCTCGC,2.28311187737
CGTGATGCAACTGTGTAATGCGGCTGAC,27.8665335936
GAACATCACCGCCTCACGTCCGGTTTTG,26.3480104838
TCGATTGAGGGGCTGGGAATGGGTGATC,41.2590972746
CCGTGTAAGGGAGATTACACAGGCTAAG,4.25926295656
GCCGCACTGATGCCTGCATTGGGCAGCA,30.1979667509
AACCCGCTGCAAGACCACGAAGGAGAGA,27.5492041756
GTCTTCCTGTCCGTCGATATCGGTTGCA,36.6200060519
CTGCTGATGATAAGCAAGAGAGGCGAGT,19.0683445911
AAACATCGGAAGCCTGAATGCGGATACG,13.027694508
```
