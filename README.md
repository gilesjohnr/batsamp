[![DOI](https://zenodo.org/badge/289552613.svg)](https://zenodo.org/badge/latestdoi/289552613)

## Code and data to reproduce analyses in 'Optimizing Non-Invasive Sampling of an Infectious Bat Virus'

This repository provides R code to reproduce models and simulations in the following manuscript:

Giles, John R., Alison J. Peel, Konstans Wells, Raina K. Plowright, Hamish McCallum, and Olivier Restif. 2018. “Optimizing Non-Invasive Sampling of an Infectious Bat Virus.” bioRxiv. [https://doi.org/10.1101/401968](https://www.biorxiv.org/content/10.1101/401968v1.abstract).

Data of Hendra virus prevalence used here are redistributed under the Creative Commons Attribution 4.0 license. Original data sources are made publicly available by the Queensland Department of Agriculture and Fisheries and New South Wales Department of Primary Industries here:

  * [Individual-level prevalence data](https://www.data.qld.gov.au/dataset/hendra-virus-infection-in-australian-black-flying-foxes-pteropus-alecto)
  * [Roost-level prevalence data](https://www.data.qld.gov.au/dataset/hev-infection-flying-foxes-eastern-australia)


### Abstract
Outbreaks of infectious viruses resulting from spillover events from bats have brought much attention to the ecological origins of bat-borne zoonoses, resulting in an increase in ecological and epidemiological studies on bat populations. Field sampling methods often collect pooled samples of bat excreta from large plastic sheets placed under roosts. However, because multiple individuals may contribute to a pooled sample, positive bias is introduced, making studies of viral dynamics difficult. We therefore explored the general issue of bias in spatial sample pooling using Hendra virus in Australian bats as a case study.

We analyzed bias in different under-roost sampling designs by fitting Generalized Additive Models to field data from individually captured bats and pooled urine samples. We then used theoretical simulation models of bat density and under-roost sampling to understand the mechanistic drivers of bias and assessed the respective accuracy of four sampling designs to estimate viral prevalence.

The design used in most field studies estimated viral prevalence with a magnitude 3.2 times higher than individual-level data, with positive bias 5--7 times higher than other designs, which resulted from spatial auto-correlation among sampling sheets and clustering of bats in the roost. Given plausible population sizes of bats, our simulation models indicate using a stratified random design to collect 30--40 pooled urine samples from 80--100 sheets, each with an area of 0.75--1m$^2$, would allow estimation of the true prevalence with minimum sampling bias and false negatives.

These results show that widely used under-roost sampling techniques are highly sensitive to viral presence, but lack specificity, providing limited information regarding viral dynamics. Our simulation results indicate that substantially improved estimation of true prevalence requires specific, yet minor, changes to existing sampling approaches such as reducing sheet size, increasing the number of sheets, and spreading sheets out within the roost area. Our findings provide insight into how spatial sample pooling is vulnerable to bias for a wide range of systems in disease ecology, where optimal sampling design is influenced by three things: individual prevalence, population density and patterns of aggregation.

### Troubleshooting

For questions contact package maintainers John Giles (gilesjohnr@gmail.com).

### License

Code here are available under the GNU General Public License v3.0.
