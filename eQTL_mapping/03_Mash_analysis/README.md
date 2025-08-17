
## 3. Multivariate adaptive shrinkage analysis (*mash*)

The genetic variant effects on the expression of genes in macrophages treated with IFN-gamma, sLPS, and Ctrl, in 6 and 24 hrs each, were obtained from [Panousis et al., 2023](https://www.biorxiv.org/content/10.1101/2023.05.29.542425v1). *Mash* was run on these macrophage and our microglia eQTLs to improve effect estimates by learning effect sharing patterns among conditions in both cell types. 

### 3.0 Prepare input data for *mash*
Script [`00_prepare_input_data.R`](00_prepare_input_data.R). 

**Step 1: define strong eQTL set across all microglia and macrophage conditions**

- *Step 1.1: find significant lead eQTLs per cell type*

  The lead eQTLs that were significant in at least one condition (q <0.05) were extracted separately for microglia (subsetting to non-proliferative cluster) and macrophages. Within each cell type and condition, all lead eQTLs were unique and all genes had only one lead variant. In microglia there were 12,500 unique eQTLs lead and significant in at least one treatment, whereas in macrophages 19,902 unique eQTLs were lead and significant in at least one of the 6 conditions (treatment x hour).

  Because the variant IDs of eQTLs were given with `[Ref]-[Alt]` alleles in microglia and with variant `rsID` in macrophages, we searched the rsIDs of the 12,500 strong microglia eQTLs in [`find_variant_rsIDs.sh`](find_variant_rsIDs.sh) (by variant chr and position) to later search them in macrophage nominal data (*Step 1.2*). After discarding the ones without rsID available (not findable among macrophage conditions) or with shared rsID (= ambiguous cases where we cannot resolve which allele the rsID in macrophage refers to in microglia), 9,160 microglia strong eQTLs remained. 
  
  Analogously, we mapped the `Ref` and `Alt` allele(s) for each variant in the 19,902 macrophage eQTLs (by using variant chr and position) in [`find_variant_alleles.sh`](find_variant_alleles.sh) script. After extracting the allele information we expanded all eQTLs from `[geneID]`-`[chr]`-`[pos]`-`[rsID]` to `[geneID]`-`[chr]`-`[pos]`-`[Ref]`-`[Alt]`. Variants with no alleles were discarded as they cannot be found among microglia conditions anyways. All `Alt` alleles of multiallelic variants were incorporated in the eQTL IDs so that we could search them all in the microglia nominal data (see *Step 1.2*). 
  
  After allele expansion, 35,338 unique allele-based eQTLs were obtained. Of the multiple `Alt` alleles of multiallelic variants, we kept the single one present in the microglia nominal data. Macrophage multiallelic variants with >1 allele present in microglia data were discarded because for those we don't know which one was used in Panousis et al. In total, 9,937 macrophage strong eQTLs were kept. 
  
  The union of the 9,160 microglia and the 9,937 macrophage lead signif. eQTLs gave the initial set of strong eQTLs: 18,491 strong eQTLs that are lead and significant in at least one condition in microglia or macrophages. 


- *Step 1.2: extract effects and se for strong eQTLs across all conditions and cell types*

  We first joined all eQTLs with effect and std. error across the three treatments in microglia, giving 11,251,961 unique eQTLs. Of these, we removed the ones sharing `[geneID]-[chr]-[pos]-[Ref]` as these are eQTLs from multiallelic variants for which we captured >1 `Alt` allele in our microglia data, therefore ambiguous when mapped in the macromap data. A total of 11,021,081 eQTLs remained after this step.

  Next, we searched the rsID of such eQTLs by variant chr and position in [`find_variant_rsIDs.sh`](find_variant_rsIDs.sh). Again, eQTLs without/sharing variant rsIDs were removed, leaving 9,688,468 eQTLs. 

  Secondly, all 63,348,409 macrophage eQTLs with effects across all 6 conditions were joined. The 9,688,468 microglia eQTLs for all 3 treatments were searched among these macrophage eQTLs based on their rsIDs, giving 8,719,455 eQTLs with effects and standard errors across all 9 microglia/macrophage conditions. All of these had only one `Alt` allele per `rsID`. 

  Utilizing this set, we searched the 18,491 eQTLs in the strong set, finding **16,427 eQTLs**.


**Step 2: extract random sample from nominal eQTL data**

*Mash* requires to first fit the MVN model using all tests or a large random subset of them (not the strong signals) to learn **data sparsity**, i.e., that many tests are zero effects/ null, and “shrinks” (“corrects”) posterior estimates towards 0. This is particularly relevant in cases where most effects are null or nearly null.

For that, we selected 5 random eQTLs for each of the genes included in the strong set, resulting in a set of **32,165 eQTLs**. 




### 3.1 Run *mash*
Script [`01_run_mash.R`](01_run_mash.R). 

#### 3.1.1 Set up covariance matrices
**3.1.1.1 Canonical covariance matrices**

Canonical covariance matrices for simple effect sharing patterns among conditions were computed. The correlation between conditions was estimated with the null tests (32K random eQTLs generated in **3.0 Prepare input data for *mash***) to account for non-independent measurements across conditions. 


**3.1.1.2 Data-driven covariance matrices**

*3.1.1.2.1 PCA for initial data-driven matrices*

Initial PC-based covariance matrices were computed by Single Value Decomposition with the set of 16K strong eQTLs.
PCA was run to further explore eQTL differences across conditions. 

*3.1.1.2.2 More refined data-driven matrices with ED*

The PC-based matrices from *step 3.1.1.2.1* were used to initialize the Extreme Deconvolution (ED) algorithm to obtain more refined data-driven matrices used to fit the model in the next step.


#### 3.1.2 Fit the mixture MVN model
The mixture MVN model was fitted to the random eQTL set with the canonical and (strong) data-driven covariance matrices, and incorporating the correlation matrix obtained in **step 3.1.1.1**. The mixture proportions for each cov matrix - scaling factor combination were returned.


#### 3.1.3 Extract posterior summaries
The posterior summaries for the strong eQTLs were calculated based on the mixture model *mash* learned with all the random tests in **3.1.2 Fit the mixture MVN model** (same mixture proportions). 

The proportion of significant effects shared (by sign and magniture) by each pair of conditions was calculated.

