
## 4. Burden tests for line proliferation efficiency


### 4.1 Perform burden tests for proliferation

Script [`01_Burden_tests.R`](01_Burden_tests.R). 

Throughout the differentiation process the changes in cell line proportions in each pool were computed in `/lustre/scratch123/hgi/mdt1/projects/otar2065/OTAR2065_WGS_iPSC_premac_micro/code/OTAR2065_WGS_iPSC_premac_micro/`. These line proliferation measurements were regressed on the cell line burdens of deleterious variants (and missense non-deleterious and synonymous variants as negative controls), plus line- and donor-level covariates. Linear mixed models were fitted and the effect size and significance of the genetic burdens were examined.  

#### 4.1.0 Explore and process proliferation data

Proliferation estimates are given for the following three comparisons:

- from iPSC → young macrophage precursor (preMAC)
- from young → old preMac 
- from (young - old) preMac → microglia (in IFN, LPS, and untreated)

Proliferation was calculated as the log-[proportion of line L in stage S2 in pool P / proportion of line L in stage S1 in pool P], scaled by the mean of the pool P (`scaled_log_fraction_[comparison]`). For each comparison, per line per pool (per batch in some pools), we have the log-fractions. In preMac → microglia, in addition, we have the log-fractions per microglia treatment and preMac age. 

Precomputed genotype principal components (accounting for population stratification of the donors) in `/OTAR2065_sc_eQTL/data/genotype/plink_genotypes/all_pools.no_outliers.genotype.MAF05.eigenvec` were added to the cell line data.


#### 4.1.1 Compute genetic burdens per line  
The numbers of deleterious (del), protein-truncating (ptv), and synonymous (syn) variants per gene per line were available in `/lustre/scratch123/hgi/projects/otar2065/resources/Puigdevall_Neuroseq_efficiency_2023/mutBurdenTabs/`.  

##### 4.1.1.1 Global genetic burdens per line
The burden of del, ptv, and syn variants per line across all genes was computed. 

##### 4.1.1.2 Genetic burdens per line in SKAT-O genes
Previous SKAT-O analyses were performed by Marta Perez to find genes whose individual deleterious variants, or their deleterious burden, associate with changes in proliferation from one differentiation stage to another. The burdens of del, ptv, and syn variants per line in the genes that were significant in these tests for each comparison (for deleterious variants), were computed. 

##### 4.1.1.3 Genetic burdens per line in signif genes of burden tests
Burden tests of deleterious variants were performed by Marta Perez on the significant SKAT-O genes to obtain proliferation-associated genes whose variants are likely to go in the same direction. The burdens of del, ptv, and syn variants per line in the genes that were significant in these tests for each comparison, were computed. In addition, the significant genes were separated by effect directionality (positive and negative beta) and their burdens were aggregated separately per line. 

##### 4.1.1.4 Add line burdens to proliferation data
All burdens were added to the lines for which proliferation estimates exist. 



#### 4.1.2 Fit linear mixed models
Burden tests were performed regressing the `scaled_log_fraction_[comparison]` on the computed burdens in **4.1.1 Compute genetic burdens per line** across lines, using `lmerTest::lmer()`. The following were the linear mixed models used to assess:

The effect of the **global** burden of del / ptv / syn mutations per line on:

- proliferation from iPSC → preMac:
    
    `scaled_log_fraction_premac_vs_iPSC` ~ `Global_Burden_[Del/Ptv/Syn]` + `sex` + `prop_unadjusted_min_value` + `PC1` + `PC2` + `(1 | pool)`
    
- proliferation from young → old preMac:
    
    `scaled_log_fraction_old_vs_young_premac` ~ `Global_Burden_[Del/Ptv/Syn]` + `sex` + `prop_unadjusted_min_value` + `PC1` + `PC2` + `(1 | pool)`
    
- proliferation from preMac → microglia (across different preMac ages and microglia treatments):
    
    `scaled_log_fraction_microglia_vs_premac` ~ `Global_Burden_[Del/Ptv/Syn]` + `sex` + `prop_unadjusted_min_value` + `PC1` + `PC2` + `(1 | pool)` + `(1|treatment)`

The effect of the burden of del / ptv / syn mutations per line in **SKAT-O** significant genes associated with iPSC → preMac / young → old preMac / preMac → microglia, based on their **deleterious** variants, on:

- proliferation from iPSC → preMac:
    
    `scaled_log_fraction_premac_vs_iPSC`  ~ `SKAT0_del_premac_vs_iPSC_Burden_[Del/Ptv/Syn]` + `sex` + `prop_unadjusted_min_value` + `PC1` + `PC2` + `(1 | pool)`
    
- proliferation from young → old preMac:
    
    `scaled_log_fraction_old_vs_young_premac` ~ `SKAT0_del_old_vs_young_premac_Burden_[Del/Ptv/Syn]` + `sex` + `prop_unadjusted_min_value` + `PC1` + `PC2` + `(1 | pool)`
    
- proliferation from preMac → microglia:
    
    `scaled_log_fraction_microglia_vs_premac` ~ `SKAT0_del_microglia_vs_premac_Burden_[Del/Ptv/Syn]` + `sex` + `prop_unadjusted_min_value` + `PC1` + `PC2` + `(1 | pool)` + `(1|treatment)`
    

The effect of the burden of del / ptv / syn mutations per line in significant genes associated with iPSC → preMac / young → old preMac / preMac → microglia, based on their **deleterious burden**, on:

- proliferation from iPSC → preMac:

  `scaled_log_fraction_premac_vs_iPSC`  ~ `Burden_del_premac_vs_iPSC_Burden_[Del/Ptv/Syn]` + `sex` + `prop_unadjusted_min_value` + `PC1` + `PC2` + `(1 | pool)`

- proliferation from young → old preMac:

  `scaled_log_fraction_old_vs_young_premac` ~ `Burden_del_old_vs_young_premac_Burden_[Del/Ptv/Syn]` + `sex` + `prop_unadjusted_min_value` + `PC1` + `PC2` + `(1 | pool)`

- proliferation from preMac → microglia:

  `scaled_log_fraction_microglia_vs_premac` ~ `Burden_del_microglia_vs_premac_Burden_[Del/Ptv/Syn]` + `sex` + `prop_unadjusted_min_value` + `PC1` + `PC2` + `(1 | pool)` + `(1|treatment)`


The effect of the burden of del / ptv / syn mutations per line in significant genes **positively** (`pos`) OR **negatively** (`neg`) associated with iPSC → preMac / young → old preMac / preMac → microglia, based on their **deleterious burden**, on:

- proliferation from iPSC → preMac:

  `scaled_log_fraction_premac_vs_iPSC`  ~ `Burden_del_premac_vs_iPSC_[pos/neg]_Burden_[Del/Ptv/Syn]` + `sex` + `prop_unadjusted_min_value` + `PC1` + `PC2` + `(1 | pool)`
  
- proliferation from young → old preMac:

  `scaled_log_fraction_old_vs_young_premac` ~ `Burden_del_old_vs_young_premac_[pos/neg]_Burden_[Del/Ptv/Syn]` + `sex` + `prop_unadjusted_min_value` + `PC1` + `PC2` + `(1 | pool)`

- proliferation from preMac → microglia:

  `scaled_log_fraction_microglia_vs_premac` ~ `Burden_del**_microglia_vs_premac_[pos/neg]_Burden_[Del/Ptv/Syn]` + `sex` + `prop_unadjusted_min_value` + `PC1` + `PC2` + `(1 | pool)` + `(1|treatment)`



### 4.2 Functional enrichment analysis for proliferation-associated genes

Script [`02_Functional_enrichment`](02_Functional_enrichment.R). 

Biological functionality support was added to our proliferation-associated genes through functional enrichment tests.

Sets of functional genes considered:

- Essential human genes from [DepMap](https://depmap.org/portal/) implicated in cell viability (n = 1,247).
- Human genes with mouse ortholog(s) implicated in macrophage survival from [Covarrubias et al. (2020)](https://doi.org/10.1016/j.celrep.2020.108541) through:
    - MAGeCK analysis (n = 417 significant genes in mouse; 393 human orthologs).
    - Mann–Whitney test (n = 21,377 assessed mouse genes; 15,815 human genes with mouse ortholog(s) implicated in macrophage survival, of which 3,175 were significant).
    - Union of MAGeCK and Mann-Whitney significant genes in human (n = 3,212).
    
#### 4.2.1 Overrepresentation analysis (ORA)

The following were the enrichment tests performed:

**Enrichment among SKAT-O genes**

- ***SKAT-O for deleterious variants***
    - SKAT-O genes associated with any prolif. comparison: the enrichment of the 4 sets among genes that were significant in the SKAT-O test for **deleterious** variants for any of the comparisons tested (iPSC → young preMAC, young → old preMac, OR preMac → microglia).
    - SKAT-O genes associated with prolif. from iPSC to premac: same as above but for genes significantly associated with iPSC → young preMAC.
    - SKAT-O genes associated with prolif. from young to old premac: same as above but for genes significantly associated with young → old preMac.
    - SKAT-O genes associated with prolif. from premac to microglia: same as above but for genes significantly associated with preMac → microglia.
    
- ***SKAT-O for non-deleterious (ptv) variants***
    - SKAT-O genes associated with any prolif. comparison: the enrichment of the 4 sets among genes that were significant in the SKAT-O test for **non-deleterious** variants for any of the comparisons tested (iPSC → young preMAC, young → old preMac, OR preMac → microglia).
    - SKAT-O genes associated with prolif. from iPSC to premac: same as above but for genes significantly associated with iPSC → young preMAC.
    - SKAT-O genes associated with prolif. from young to old premac: same as above but for genes significantly associated with young → old preMac.
    - SKAT-O genes associated with prolif. from premac to microglia: same as above but for genes significantly associated with preMac → microglia.

**Enrichment among Burden test genes**

- ***Burden tests for deleterious variants***
    - Burden test genes associated with any prolif. comparison: the same 4 gene sets were tested for their enrichment among genes that were significantly associated with any proliferation based on their **deleterious** burdens.
    - Burden test genes associated with prolif. from iPSC to premac: same as above but for genes significantly associated with iPSC → young preMAC.
    - Burden test genes associated with prolif. from young to old premac: same as above but for genes significantly associated with young → old preMac.
    - Burden test genes associated with prolif. from premac to microglia: same as above but for genes significantly associated with preMac → microglia.

- ***Burden tests for non-deleterious (ptv) variants***
    - Burden test genes associated with any prolif. comparison: the same 4 gene sets were tested for their enrichment among genes that were significantly associated with any proliferation based on their **non-deleterious** burdens.
    - Burden test genes associated with prolif. from iPSC to premac: same as above but for genes significantly associated with iPSC → young preMAC.
    - Burden test genes associated with prolif. from premac to microglia: same as above but for genes significantly associated with preMac → microglia.

- ***Burden tests for synonymous variants***
    - Burden test genes associated with any prolif. comparison: the same 4 gene sets were tested for their enrichment among genes that were significantly associated with any proliferation based on their **synonymous** burdens.
    - Burden test genes associated with prolif. from premac to microglia: same as above but for genes significantly associated with preMac → microglia.
    
*Note*: enrichment among genes with significant ptv and syn burdens are missing for some proliferation comparisons because these burden tests didn't result in significant genes. 


**Enrichment among SKAT-O signif negative genes**

- ***SKAT-O for deleterious variants***
    - SKAT-O genes decreasing any prolif. comparison: the enrichment of the previous 4 sets (+ Mann-Whitney genes divided by sign) among genes that were significant in the SKAT-O test for **deleterious** variants for any of the comparisons tested (iPSC → young preMAC, young → old preMac, OR preMac → microglia), and with **negative** beta in the burden test.
    - SKAT-O genes decreasing prolif. from iPSC to premac: same as above but for genes significantly decreasing iPSC → young preMAC.
    - SKAT-O genes decreasing prolif. from young to old premac: same as above but for genes significantly decreasing young → old preMac.
    - SKAT-O genes decreasing prolif. from premac to microglia: same as above but for genes significantly decreasing preMac → microglia.
    

**Enrichment among SKAT-O signif positive genes**

- ***SKAT-O for deleterious variants***
    - SKAT-O genes increasing any prolif. comparison: the enrichment of the previous 6 sets among genes that were significant in the SKAT-O test for **deleterious** variants for any of the comparisons tested (iPSC → young preMAC, young → old preMac, OR preMac → microglia), and with **positive** beta in the burden test.
    - SKAT-O genes increasing prolif. from iPSC to premac: same as above but for genes significantly increasing iPSC → young preMAC.
    - SKAT-O genes increasing prolif. from young to old premac: same as above but for genes significantly increasing young → old preMac.
    - SKAT-O genes increasing prolif. from premac to microglia: same as above but for genes significantly increasing preMac → microglia.
  

**Enrichment among Mann–Whitney genes** 

- The enrichment of **SKAT-O** and **Burden** test significant genes for any comparison and each one individually based on **deleterious** variants, were assessed for their enrichment among the significant macrophage genes in Covarrubias et al. according to the Mann-Whitney test.



#### 4.2.2 Gene Set Enrichment Analysis (GSEA)

**Enrichment among ranked SKAT-O genes**

- ***SKAT-O test for deleterious variants***
    - SKAT-O genes ranked for association with prolif. from iPSC to premac: genes assessed in SKAT-O for iPSC → young preMAC based on **deleterious** variants were ranked by decreasing -log10(adj. *p*) and the previous gene sets were tested for enrichment.
    - SKAT-O genes ranked for association with prolif. from young to old premac: same as above but for genes assessed in SKAT-O for young → old preMac.
    - SKAT-O genes ranked for association with prolif. from premac to microglia: same as above but for genes assessed in SKAT-O for preMac → microglia.


**Enrichment among ranked Burden test genes**

- ***Burden test for deleterious variants***
    - Burden test genes ranked for association with prolif. from iPSC to premac: genes assessed in burden tests for iPSC → young preMAC based on **deleterious** variants were ranked by decreasing -log10(adj. *p*) x signed coeff. and the previous gene sets were tested for enrichment.
    - Burden test genes ranked for association with prolif. from young to old premac: same as above but for genes assessed in burden tests for young → old preMac.
    - Burden test genes ranked for association with prolif. from premac to microglia: same as above but for genes assessed in burden tests for preMac → microglia.
    
    
**Enrichment among ranked Covarrubias MAGeCK genes**

- The genes analyzed in Covarrrubias et al. with MAGeCK were ranked by decreasing -log10(adj. *p*) and the sets of significant **SKAT-O** and **Burden** test genes based on **deleterious** variants were assessed for enrichment.

**Enrichment among ranked Covarrubias Mann-Whitney genes** 

- The genes tested in Covarrrubias et al. with Mann-Whitney were ranked by decreasing -log10(*p*) x coeff. and the sets of significant **SKAT-O** and **Burden** test genes based on **deleterious** variants were assessed for enrichment.


