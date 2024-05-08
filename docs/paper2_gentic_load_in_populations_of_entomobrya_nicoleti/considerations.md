# CONSIDERATIONS

## Genetic load measures

N<sub>del</sub>: Number of deleterious mutations segregating in population (MAF > 0)

L<sub>drift</sub>: Number of deleterious mutations that are nearly fixed in population (MAF > 0.9)

Extracting data from VCF in CDS regions with DoS value in bottom 10%

| CHROM | POS | REF | ALT | AA | AF | CLASS | |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| HiC_scaffold_7  | 16461637 | C | A | C | 0.03 | stop_gained | <- Ignored as stop gained is its own category, part of LoF |
| HiC_scaffold_7 | 16461571 | A | T | . | 0.09 | missense_variant | <- Ignored as ancestral allele could not be determined |
| HiC_scaffold_7 | 16461452 | TT | TA | T | 0.08 | missense_variant | <- Assume ancestral state and refence are identical, count towards N<sub>del</sub> as frequency > 0 |
| HiC_scaffold_7 | 16461091 | ATG | ACG | A | 0.03 | missense_variant | <- Ignore? |
| HiC_scaffold_1 | 19028065 | T | G | T | 0.98 | missense_variant | <- Count towards N<sub>del</sub> and L<sub>drift</sub> as frequency > 0.9 |
| HiC_scaffold_7 | 16461577 | G | C | G | 0.03 | missense_variant | <- Count towards N<sub>del</sub> as frequency > 0 |

## SnpEff

## Effect sequence ontology

missense_variant:  
Variant causes a codon that produces a different amino acid. MODERATE

stop_lost:  
Variant causes stop codon to be mutated into a non-stop codon. HIGH

start_lost:  
Variant causes start codon to be mutated into a non-start codon. HIGH

stop_gained:  
Variant causes a STOP codon. HIGH

synonymous_variant:  
Variant causes a codon that produces the same amino acid. LOW

start_retained:  
Variant causes start codon to be mutated into another start codon. LOW

stop_retained_variant:  
Variant causes stop codon to be mutated into another stop codon. LOW

### Impact

HIGH:  
The variant is assumed to have high (disruptive) impact in the protein, probably causing protein truncation, loss of function or triggering nonsense mediated decay.

MODERATE:  
A non-disruptive variant that might change protein effectiveness.

LOW:  
Assumed to be mostly harmless or unlikely to change protein behavior.