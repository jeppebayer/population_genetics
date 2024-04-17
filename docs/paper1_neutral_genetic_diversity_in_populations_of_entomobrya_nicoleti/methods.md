# METHODS

## Reference genome assembly and annotation

Order
HiFiAdapterFilt


The genome assembly was carried out using Hifiasm (REF) (parameters????). We identified and removed haplotypic duplications using purge_dups (REF) before scaffolding assembled contigs with Hi-C data using JUICEBOX (REF).

The protein coding genes were annotated using Braker3 (REF). RNA data was generated for this study, while protein coding amino acid sequences from Orchesella cincta ??? were downloaded from ncbi.

Repetitive content was annotated by ...

### References

**HiFiAdapterFilt**  
Sim, S.B., Corpuz, R.L., Simmonds, T.J. *et al.* HiFiAdapterFilt, a memory efficient read processing pipeline, prevents occurrence of adapter sequence in PacBio HiFi reads and their negative impacts on genome assembly. *BMC Genomics* **23**, 157 (2022). <https://doi.org/10.1186/s12864-022-08375-1>

**hifiasm**  
Cheng, H., Concepcion, G.T., Feng, X. *et al.* Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nat Methods* **18**, 170–175 (2021). <https://doi.org/10.1038/s41592-020-01056-5>  
Cheng, H., Jarvis, E.D., Fedrigo, O. *et al.* Haplotype-resolved assembly of diploid genomes without parental data. *Nat Biotechnol* **40**, 1332–1335 (2022). <https://doi.org/10.1038/s41587-022-01261-x>

**BUSCO**  
Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov, BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes, *Molecular Biology and Evolution*, Volume 38, Issue 10, October 2021, Pages 4647–4654, <https://doi.org/10.1093/molbev/msab199>

**purge_dups**
Dengfeng Guan, Shane A McCarthy, Jonathan Wood, Kerstin Howe, Yadong Wang, Richard Durbin, Identifying and removing haplotypic duplication in primary genome assemblies, *Bioinformatics*, Volume 36, Issue 9, May 2020, Pages 2896–2898. <https://doi.org/10.1093/bioinformatics/btaa025>

**Juicer**  
Neva C. Durand, Muhammad S. Shamim, Ido Machol, Suhas S. P. Rao, Miriam H. Huntley, Eric S. Lander, and Erez Lieberman Aiden. "Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments." *Cell Systems* 3(1), 2016. <https://doi.org/10.1016/j.cels.2016.07.002>

**3D-DNA**  
Dudchenko, O., Batra, S.S., Omer, A.D., Nyquist, S.K., Hoeger, M., Durand, N.C., Shamim, M.S., Machol, I., Lander, E.S., Aiden, A.P., et al. (2017). *De novo assembly of the Aedes aegypti genome using Hi-C yields chromosome-length scaffolds*. Science. Apr 7; 356(6333):92-95. doi: <https://doi.org/10.1126/science.aal3327>. Epub 2017 Mar 23.

**JuiceBox**  
Durand NC, Robinson JT, Shamim MS, Machol I, Mesirov JP, Lander ES, Aiden EL. Juicebox Provides a Visualization System for Hi-C Contact Maps with Unlimited Zoom. *Cell Syst.* 2016 Jul;3(1):99-101. <https://doi.org/10.1016/j.cels.2015.07.012>

**RepeatModeler**  
 Flynn JM, Hubley R, Goubert C, Rosen J, Clark AG, Feschotte C, Smit AF (2020) RepeatModeler2 for automated genomic discovery of transposable element families. Proceedings of the National Academy of Sciences of the United States of America 117: 9451–9457. <https://doi.org/10.1073/pnas.1921046117>

**SeqKit**  
Wei Shen*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. [doi:10.1002/imt2.191](https://doi.org/10.1002/imt2.191)  
Wei Shen, Shuai Le, Yan Li*, and Fuquan Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLOS ONE. [doi:10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962)

**RepBase**  
Bao W, Kojima KK, Kohany O (2015) Repbase Update, a database of repetitive elements in eukaryotic genomes. Mobile DNA-UK 6: 1–11. <https://doi.org/10.1186/s13100-015-0041-9>

**RepeatMasker**  
Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0.
2013-2015 <http://www.repeatmasker.org>

**STAR**  
Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: [10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635). Epub 2012 Oct 25. PMID: 23104886; PMCID: [PMC3530905](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc3530905/).

**AGAT** (change to relevant version)  
Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  
(Version v0.8.0). Zenodo. <https://www.doi.org/10.5281/zenodo.3552717>

**BRAKER** (Specific to this run)  
Hoff, K. J., Lange, S., Lomsadze, A., Borodovsky, M., & Stanke, M. (2016). BRAKER1: unsupervised RNA-Seq-based genome annotation with GeneMark-ET and AUGUSTUS. Bioinformatics, 32(5), 767-769.  
Bruna, T., Hoff, K.J., Lomsadze, A., Stanke, M., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR Genomics and Bioinformatics 3(1), lqaa108.  
Hoff, K. J., Lomsadze, A., Borodovsky, M., & Stanke, M. (2019). Whole-genome annotation with BRAKER. In Gene Prediction (pp. 65-95). Humana, New York, NY.  
Bruna, T., Lomsadze, A., & Borodovsky, M. (2020). GeneMark-EP+: eukaryotic gene prediction with self-training in the space of genes and proteins. NAR Genomics and Bioinformatics, 2(2), lqaa026.  
Lomsadze, A., Ter-Hovhannisyan, V., Chernoff, Y. O., & Borodovsky, M. (2005). Gene identification in novel eukaryotic genomes by self-training algorithm. Nucleic acids research, 33(20), 6494-6506.  
Buchfink, B., Xie, C., & Huson, D. H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature Methods, 12(1), 59.  
Gotoh, O. (2008). A space-efficient and accurate method for mapping and aligning cDNA sequences onto genomic sequence. Nucleic acids research, 36(8), 2630-2638.  
Iwata, H., & Gotoh, O. (2012). Benchmarking spliced alignment programs including Spaln2, an extended version of Spaln that incorporates additional species-specific features. Nucleic acids research, 40(20), e161-e161.  
Stanke, M., Diekhans, M., Baertsch, R., & Haussler, D. (2008). Using native and syntenically mapped cDNA alignments to improve de novo gene finding. Bioinformatics, 24(5), 637-644.  
Stanke, M., Schöffmann, O., Morgenstern, B., & Waack, S. (2006). Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources. BMC Bioinformatics, 7(1), 62.  
Gabriel, L., Hoff, K. J., Bruna, T., Borodovsky, M., & Stanke, M. (2021). TSEBRA: transcript selector for BRAKER. BMC Bioinformatics, 22:566.

## Mapping, filtering and variant calling of resequencing data

Prior to mapping, the re-sequencing data processed with AdapterRemoval (Schubert et al. 2016) to search for and remove any residuals of known adapter sequences. Furthermore, reads were trimmed at the 5' and 3' termini for Ns and bases with a quality lower than 25, reads shorter than 20 base-pairs were also discarded, and for paired-end reads overlapping mates were merged into a single read and the base quality was recalculated if the overlap was at least 11 base pairs in length (default) with a maximum mismatch ratio of 1/3 (default). The re-sequencing data was mapped to reference genomes using the BWA-MEM algorithm from BWA (Li H. 2013).  

FILTERING

Variants were called using Freebayes (VCF files). PARAMETERS ...

VCF files were filtered by 1) removing variants 5bp on each side of indel polymorphisms, 2) removing sites with coverage lower than 300 and higher than 600, and 3) removing sites that were not bi-allelic.

For all subsequent analyses we used genomic regions not including protein coding genes (exons and introns) and repetitive regions. This was done to analyze regions as neutral as possible to avoid patterns generated by selection.

### References

**BWA**  
Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. [arXiv:1303.3997v2](http://arxiv.org/abs/1303.3997) [q-bio.GN]
