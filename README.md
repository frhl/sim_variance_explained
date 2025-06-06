# Interactive Genetic Architecture Explorer

This Shiny app goes along with our paper **"Deviations from genetic additivity driven by rare variants at biobank scale"** by Lassen et al. It lets you play around with different genetic architectures and see what happens when things don't follow simple additive patterns.

## What it does

- **Two ways to set up your data**: 
  - **HWE Mode**: The usual way, assuming Hardy-Weinberg equilibrium
  - **Custom Counts**: Set your own genotype counts (no HWE required!)
- **HWE Testing**: Shows you in real-time if your population is violating HWE assumptions
- **Pretty plots**: Visualizes how genotype relates to phenotype, plus shows you when additive models break down
- **Variance breakdown**: See how much variance comes from additive vs non-additive effects

## The main genetic patterns you can explore

1. **Additive**: Het effect is half the homozygote effect - boring but common
2. **Partially Recessive**: Hets have some effect, but homozygotes have disproportionately more
3. **Mendelian Recessive**: Hets do nothing, only homozygotes matter - this breaks additive models badly

## Fun things to try

**Rare disease scenario (like in our paper)**
```
Set custom counts: 200,000 normals, 100 hets, 5 homozygotes
Effects: [0, 0, 2]
```
This is what real rare variant data looks like in biobanks!

**Breaking HWE on purpose**
```
Try: 180,000 normals, 1,000 hets, 1,000 homozygotes  
```
Watch the app tell you "HWE VIOLATED" - this is why we needed new methods!

**Regular trait that's kinda recessive**
```
Use HWE mode with MAF=0.3 and effects [0, 0.5, 2]
```

## How to run it

```r
install.packages(c("shiny", "ggplot2", "dplyr", "patchwork", "bslib"))
# Then just run the code!
```

## Why this matters for our research

The whole point of our paper was that existing methods assume HWE, but rare variants often violate this. This app shows you:

- Why you can't just use allele frequencies when HWE doesn't hold
- How badly additive models fail for rare recessive stuff  
- What our orthogonal encoding actually captures (the red vs blue bars)
- Why detecting rare variant effects is so hard in the first place

Basically, if you want to understand why we had to develop new methods, play around with this for a bit!

## Citation

> Lassen, F.H., Venkatesh, S.S., Baya, N.A., Lindgren, C.M., and Palmer, D.S. "Deviations from genetic additivity driven by rare variants at biobank scale."
