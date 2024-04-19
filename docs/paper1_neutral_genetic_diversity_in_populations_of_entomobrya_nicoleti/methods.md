# METHODS

- [METHODS](#methods)
	- [Reference genome assembly and annotation](#reference-genome-assembly-and-annotation)
		- [Genome assembly](#genome-assembly)
		- [Annotation](#annotation)
			- [Repetitive content](#repetitive-content)
			- [Gene annotation](#gene-annotation)
		- [References](#references)
	- [Mapping, filtering and variant calling of resequencing data](#mapping-filtering-and-variant-calling-of-resequencing-data)
		- [Mapping](#mapping)
		- [Variant calling and filtering](#variant-calling-and-filtering)
		- [References](#references-1)
	- [Acknowledgements](#acknowledgements)

## Reference genome assembly and annotation

### Genome assembly

Before the assembly process any potential remaining adapters were removed from the PacBio HiFi sequencing data using HiFiAdapterFilt (Sim et al 2022). A draft genome assembly was carried out using hifiasm (Cheng et al 2021; Cheng et al 2022) with a similarity threshold of 0.1 and an initial purge level of 3. Further haplotypic duplication and overlaps were identified using purge_dups (Guan et al 2020) with minimap2 (Li 2018) and the resulting draft genome was assessed with BUSCO (Manni et al 2021) using arthropoda_odb10 before scaffolding the assembled contigs with Hi-C data using the Juicer pipeline (Durand et al 2016¹) (to create a Hi-C contact map) followed by the 3D-DNA pipeline (Dudchenko et al 2017) (to correct misassembles, anchor, order and orient fragments of DNA). The resulting chromosome level assembly was then reviewed and manually curated using Juicebox Assembly Tools (Durand et al 2016²).

### Annotation

#### Repetitive content

Repetitive content in the genome was identified using RepeatMasker (Smit, Hubley & Green 2013-16) combining the results from two seperate runs, with the results from RepeatModeler (Flynn et 2020) and an arthropod dataset from RepBase (Bao, Kojima & Kohany 2015) respectively, and soft-masked using BEDTools (Quinlan & Hall 2010).

#### Gene annotation

The protein coding gens were annotated using Braker (Bruna, Lomsadze & Borodovsky 2020; Bruna et al 2021; Buchfink, Xie & Huson 2015; Gabriel et al 2021; Gotoh 2008; Hoff et al 2016; Hoff et al 2019; Iwata & Gotoh 2012; Lomsadze et al 2005; Stanke et al 2006; Stanke et al 2008). RNA sequence data was generated for this study and aligned to the genome using STAR (Dobin et al 2013), however Braker found insufficient support for annotation, instead protein coding sequences were extracted using AGAT (Dainat) from Orchesella cincta data downloaded from NCBI (Genome, GCA_001718145.1) and used for the annotation.

### References

- **HiFiAdapterFilt**  
Sim, S.B., Corpuz, R.L., Simmonds, T.J. *et al*. HiFiAdapterFilt, a memory efficient read processing pipeline, prevents occurrence of adapter sequence in PacBio HiFi reads and their negative impacts on genome assembly. *BMC Genomics* **23**, 157 (2022). <https://doi.org/10.1186/s12864-022-08375-1>
- **hifiasm**  
Cheng, H., Concepcion, G.T., Feng, X. *et al.* Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nat Methods* **18**, 170–175 (2021). <https://doi.org/10.1038/s41592-020-01056-5>  
Cheng, H., Jarvis, E.D., Fedrigo, O. *et al.* Haplotype-resolved assembly of diploid genomes without parental data. *Nat Biotechnol* **40**, 1332–1335 (2022). <https://doi.org/10.1038/s41587-022-01261-x>
- **BUSCO**  
Manni, M., Berkeley, M. R., Seppey, M., Simão, F. A., Zdobnov, E. M., BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes, *Molecular Biology and Evolution*, Volume 38, Issue 10, October 2021, Pages 4647–4654, <https://doi.org/10.1093/molbev/msab199>
- **purge_dups**  
Guan, D., McCarthy, S. A., Wood, J., Howe, K., Wang, Y., Durbin, R. Identifying and removing haplotypic duplication in primary genome assemblies, *Bioinformatics*, Volume 36, Issue 9, May 2020, Pages 2896–2898. <https://doi.org/10.1093/bioinformatics/btaa025>
- **minimap2**  
Li, H. Minimap2: pairwise alignment for nucleotide sequences, *Bioinformatics*, Volume 34, Issue 18, September 2018, Pages 3094–3100, <https://doi.org/10.1093/bioinformatics/bty191>
- **Juicer**  
Durand¹, N. C., Shamim, M. S., Machol, I, Rao, S. S. P., Huntley, M. H., Lander, E. S., & Aiden, E. L. "Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments." *Cell Systems* 3(1), 2016. <https://doi.org/10.1016/j.cels.2016.07.002>
- **3D-DNA**  
Dudchenko, O., Batra, S.S., Omer, A.D., Nyquist, S.K., Hoeger, M., Durand, N.C., Shamim, M.S., Machol, I., Lander, E.S., Aiden, A.P., et al. (2017). *De novo assembly of the Aedes aegypti genome using Hi-C yields chromosome-length scaffolds*. Science. Apr 7; 356(6333):92-95. doi: <https://doi.org/10.1126/science.aal3327>. Epub 2017 Mar 23.
- **JuiceBox**  
Durand² NC, Robinson JT, Shamim MS, Machol I, Mesirov JP, Lander ES, Aiden EL. Juicebox Provides a Visualization System for Hi-C Contact Maps with Unlimited Zoom. *Cell Syst.* 2016 Jul;3(1):99-101. <https://doi.org/10.1016/j.cels.2015.07.012>
- **RepeatModeler**  
Flynn JM, Hubley R, Goubert C, Rosen J, Clark AG, Feschotte C, Smit AF (2020) RepeatModeler2 for automated genomic discovery of transposable element families. Proceedings of the National Academy of Sciences of the United States of America 117: 9451–9457. <https://doi.org/10.1073/pnas.1921046117>
- **SeqKit**  
Shen, W., Sipos, B., and Zhao L. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. [doi:10.1002/imt2.191](https://doi.org/10.1002/imt2.191)  
Shen, W., Le, S., Li, Y., and Hu, F. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLOS ONE. [doi:10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962)
- **RepBase**  
Bao W, Kojima KK, Kohany O (2015) Repbase Update, a database of repetitive elements in eukaryotic genomes. Mobile DNA-UK 6: 1–11. <https://doi.org/10.1186/s13100-015-0041-9>
- **RepeatMasker**  
Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0.
2013-2015 <http://www.repeatmasker.org>
- **BEDTools**  
Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features, *Bioinformatics*, Volume 26, Issue 6, March 2010, Pages 841–842, <https://doi.org/10.1093/bioinformatics/btq033>
- **STAR**  
Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: [10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635). Epub 2012 Oct 25. PMID: 23104886; PMCID: [PMC3530905](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc3530905/).
- **AGAT** (change to relevant version)  
Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. (Version v1.2.0). Zenodo. <https://www.doi.org/10.5281/zenodo.3552717>
- **BRAKER** (Specific to the Entomobrya nicoleti run. No RNA, only DB based on Orchesella cincta)  
Bruna, T., Lomsadze, A., & Borodovsky, M. (2020). GeneMark-EP+: eukaryotic gene prediction with self-training in the space of genes and proteins. NAR Genomics and Bioinformatics, 2(2), lqaa026.  
Bruna, T., Hoff, K.J., Lomsadze, A., Stanke, M., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR Genomics and Bioinformatics 3(1), lqaa108.  
Buchfink, B., Xie, C., & Huson, D. H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature Methods, 12(1), 59.  
Gabriel, L., Hoff, K. J., Bruna, T., Borodovsky, M., & Stanke, M. (2021). TSEBRA: transcript selector for BRAKER. BMC Bioinformatics, 22:566.  
Gotoh, O. (2008). A space-efficient and accurate method for mapping and aligning cDNA sequences onto genomic sequence. Nucleic acids research, 36(8), 2630-2638.  
Hoff, K. J., Lange, S., Lomsadze, A., Borodovsky, M., & Stanke, M. (2016). BRAKER1: unsupervised RNA-Seq-based genome annotation with GeneMark-ET and AUGUSTUS. Bioinformatics, 32(5), 767-769.  
Hoff, K. J., Lomsadze, A., Borodovsky, M., & Stanke, M. (2019). Whole-genome annotation with BRAKER. In Gene Prediction (pp. 65-95). Humana, New York, NY.  
Iwata, H., & Gotoh, O. (2012). Benchmarking spliced alignment programs including Spaln2, an extended version of Spaln that incorporates additional species-specific features. Nucleic acids research, 40(20), e161-e161.  
Lomsadze, A., Ter-Hovhannisyan, V., Chernoff, Y. O., & Borodovsky, M. (2005). Gene identification in novel eukaryotic genomes by self-training algorithm. Nucleic acids research, 33(20), 6494-6506.  
Stanke, M., Schöffmann, O., Morgenstern, B., & Waack, S. (2006). Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources. BMC Bioinformatics, 7(1), 62.  
Stanke, M., Diekhans, M., Baertsch, R., & Haussler, D. (2008). Using native and syntenically mapped cDNA alignments to improve de novo gene finding. Bioinformatics, 24(5), 637-644.
- **NCBI Orchesella cincta**  
Genome [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; 2012 – Accession No. GCA_001718145.1, Orchesella cincta, [LJIJ01](https://www.ncbi.nlm.nih.gov/nuccore/LJIJ00000000.1), 2016 [cited YYYY Mmm DD]. Available from: <https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001718145.1/>

## Mapping, filtering and variant calling of resequencing data

### Mapping

Prior to mapping, the re-sequencing data was processed with AdapterRemoval (Schubert et al. 2016) to search for and remove any residuals of known adapter sequences. Furthermore, reads were trimmed at the 5' and 3' termini for Ns and bases with a quality score lower than 25, reads shorter than 20 base-pairs were also discarded, and for paired-end reads overlapping mates were merged into a single read and the base quality was recalculated if the overlap was at least 11 base pairs in length (default) with a maximum mismatch ratio of 1/3 (default). The re-sequencing data was mapped to the reference genome using the BWA-MEM algorithm from BWA (Li 2013). Duplicates were marked and alignment files filtered using Samtools (Danecek 2021). During filtering all unmapped reads, non-primary reads, failed reads, duplicates and supplementary reads were removed (excluded bit-flag value: 3844). Assessment of alignment data quality post-filtering was done using Qualimap (García-Alcalde et al 2012).

### Variant calling and filtering

Variants were called for each of the populations using freebayes (Garrison & Garth 2012) with the parameters: -n 3 -p 100 --min-alternate-fraction 0 --min-alternate-count 2 --pooled-discrete. VCF files were then filtered using BCFTools (Danecek et al 2021) according to the following: 1) removing variants 5 bp on each side of indel polymorphisms, 2) removing all indels, 3) removing sites with coverage lower than 301 and higher than 599, and 4) removing sites that were not bi-allelic.

For all subsequent analyses we used genomic regions not including protein coding genes (exons and introns) and repetitive regions. This was done to analyze regions as neutral as possible to avoid patterns generated by selection.

### References

- **BWA**  
Li, H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. [arXiv:1303.3997v2](http://arxiv.org/abs/1303.3997) [q-bio.GN]
- **AdapterRemoval**  
Schubert, M., Lindgreen, S. & Orlando, L. AdapterRemoval v2: rapid adapter trimming, identification, and read merging. *BMC Res Notes* **9**, 88 (2016). <https://doi.org/10.1186/s13104-016-1900-2>
- **Samtools + BCFtools**  
Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., Li, H., Twelve years of SAMtools and BCFtools, *GigaScience*, Volume 10, Issue 2, February 2021, giab008, <https://doi.org/10.1093/gigascience/giab008>
- **Qualimap**  
García-Alcalde, F., Okonechnikov, K., Carbonell, J., Cruz, L. M., Götz, S., Tarazona, S., Dopazo, J., Meyer, T. F., Conesa, A., Qualimap: evaluating next-generation sequencing alignment data, *Bioinformatics*, Volume 28, Issue 20, October 2012, Pages 2678–2679, <https://doi.org/10.1093/bioinformatics/bts503>
- **Freebayes**  
Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. arXiv preprint [arXiv:1207.3907](https://arxiv.org/abs/1207.3907) [q-bio.GN] 2012
- **Popoolation2**  
Kofler, R., Pandey, R. V., Schlötterer, C., PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq), *Bioinformatics*, Volume 27, Issue 24, December 2011, Pages 3435–3436, <https://doi.org/10.1093/bioinformatics/btr589>

## Acknowledgements

All of the computing for this project was performed on the [GenomeDK](https://genome.au.dk/) cluster and workflows were created and organized using [gwf](https://gwf.app/). We would like to thank GenomeDK and Aarhus University for providing computational resources and support that contributed to these research results.
