# cell.dynamics_aging.asynchrony

This foder contains source code for data analysis in manuscript "Cellular age drives asynchronous transcriptome aging" by Ming Yang, Benjamin R. Harrison, Daniel E.L. Promislow.

## cell_demography_model

This folder contains code scripts for **Figure 1B, Supplemental Fig. S1, and S2**, based on the cell demography model we developed to model cell population dynamics in **Supplemental Methods** and cell lifespan data in mice in **Supplemental Table S1**.

## TMS_transcriptome_variation
- `00_setup_mouseID.R`, select cell types from TMS with avaible mouse cell lifespan estimates in Sender and Milo (2021)
- `01_countland_subsample.transcripts.R`, down-sampling cells to have an equal number of trancripts per tissue:cell type combination.
- `02.1_find_invariant_gene_per.mouse.R`, identity genes whose expression level didn't change significant with organismal age. A detailed pipeline is shown in **Supplemental Fig. S4**.
- `02.2_invar.gene_lmer_beta_age_per.mouse.R`, use age-invariant genes to compute gene expression variability per age group per mouse individual, **Supplemental Fig. S5**
- `03.1_change.in.expr_change.in.cell.age.R`, quantify the age effect on gene expression variability measured by coefficient variation (CV), plot the relationship between beta_age and change in mean cell age as shown in **Figure 2B**.
- `03.2_change.in.expr_change.in.cell.age_endothelia.cells.only.R`, exmamine tissue microenvironment effect in endothelial cells, **Supplemental Fig. S10**.
- `04_run_dist_spearman.cor-1_per.mouse.R`, effect of age on transcriptome variability measued by cell-cell correltaions, **Supplemental Fig. S6**.
- `04_run_GCL_decibel_gcl_per.mouse.R`, effect of age on transcriptome variability measued by GCL, **Supplemental Fig. S7**.
- `04.1_alternative.metric_log2ratio_change.in.cell.age_per.mouse.R`, integrate effect of age on transcriptome variability via different metrics, **Supplemental Table S2**.
- `05_mean.cell.age_obs.CV_per.mouse.R`, `05_mean.cell.age_obs.gcl_per.mouse.R`, and `05_mean.cell.age_obs.spearman.cor_per.mouse.R`, plot the relationship between cellular age and transcriptome variability as shown in **Figure 2A and Supplemental Fig. S8**.
- `05_mean.cell.age_obs.CV_per.mouse_female.R`, plot the relationship between cellular age and transcriptome variability in female mice as shown in **Supplemental Fig. S9**.


## GO_term_expression
GO term expression analysis for **Figure 3 and Supplemental Fig. S11**.

## dNdS_analysis

- `mouse_rat`, , dN/dS analysis with mouse-rat orthologous genes,  **Figure 4, Supplemental Fig. S12 and S13**.
- `mouse_human`, dN/dS analysis with mouse-human orthologous genes,  **Supplemental Fig. S14**.

## src

- `src_barplot_lifespan_byTissue.R`, plot the distribution of cell lifespan data used in our study as shown in **Figure 1A**.
- `src_compare_facs_droplet_depth.R`, TMS dataset overview in **Supplemental Fig. S3**
- `src_TMS_genes_id.mapping.R`, TMS data gene id mapping.
- `src_cell.demography_cell.age.difference.R`, generate cell age difference between organismal ages.
- `src_mouse_rat_biomart_to_retrieve_dnds.R`, extract mouse-rat orthologous genes.
- `src_mouse_human_biomart_to_retrieve_dnds.R`, extract mouse-human orthologous genes.
- `src_GOterms_map2_mouse.genes_BP.CC.MF.R`, extract GO terms and their assocaited gene members in mice.
