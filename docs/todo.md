# TO-DO

```bash
bcftools query -R /faststorage/project/EcoGenetics/people/Jeppe_Bayer/population_genetics/steps/collembola/Entomobrya_nicoleti/genetic_load/DoS/Grassland/EntNic_aaRJ-C225/deleterious.bed -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\t%INFO/AF\t%INFO/ANN\n' /faststorage/project/EcoGenetics/people/Jeppe_Bayer/population_genetics/steps/collembola/Entomobrya_nicoleti/genetic_load/DoS/ancestral_allele/Grassland/EntNic_aaRJ-C225/EntNic_aaRJ-C225_SnpGap_DP300_600_biallelic_AO.freebayes_n3_p100_minaltfrc0_minaltcnt2.ann.aa.vcf.gz | awk 'BEGIN{FS=OFS="\t"}{split($7, impactarray, "|"); print $1, $2, $3, $4, $5, $6, impactarray[2]}' | awk 'BEGIN{FS=OFS="\t"}{if($6 > 0.9){print}}'
```

- [ ] Draft assembly IsoVir
- [x] Draft assembly OrcCin
- [x] Draft assembly PogFla
- [ ] Hi-C scaffold IsoVir
- [ ] Hi-C scaffold OrcCin
- [ ] Hi-C scaffold PogFla
