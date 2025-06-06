# Interactive Genetic Architecture Visualization

This Shiny application complements the research paper **"Deviations from genetic additivity driven by rare variants at biobank scale"** by Lassen et al., allowing users to explore different genetic architectures and visualize deviations from additivity.

## Features

- **Interactive Controls**: Adjust minor allele frequency (0.01-0.5), phenotypic effects for each genotype, and population size
- **Dose-response Plot**: Shows relationship between genotype (0, 1, 2 copies) and phenotypic effect with additive model overlay
- **Variance Decomposition**: Displays proportion of variance explained by additive vs. non-additive components

## Key Genetic Architectures

1. **Additive**: Heterozygous effect = half homozygous effect (additive model fits perfectly)
2. **Partially Recessive**: Heterozygous effect between 0 and homozygous effect (mixed variance)  
3. **Mendelian Recessive**: Heterozygous effect = 0 (additive model fails, especially at low MAF)

## Example Settings

**Rare Recessive Disease (MAF = 0.01)**
```
Wildtype: 0, Heterozygous: 0, Homozygous: 2
```

**Common Partially Recessive Trait (MAF = 0.3)**
```
Wildtype: 0, Heterozygous: 0.5, Homozygous: 2
```

## Installation & Usage

```r
install.packages(c("shiny", "ggplot2", "dplyr", "patchwork"))
shinyApp(ui = ui, server = server)
```

## Research Connection

This tool illustrates key findings from the paper:
- Why rare recessive variants are difficult to detect with additive models
- How genetic architecture affects variance partitioning  
- The critical relationship between allele frequency and statistical power

## Citation

> Lassen, F.H., Venkatesh, S.S., Baya, N.A., Lindgren, C.M., and Palmer, D.S. "Deviations from genetic additivity driven by rare variants at biobank scale."
