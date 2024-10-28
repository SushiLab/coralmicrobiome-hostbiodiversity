# coralmicrobiome-hostbiodiversity
This is the repository associated with the manuscript **&#x00BB;Host-level biodiversity shapes the dynamics and networks within the coral reef microbiome&#x00AB;**. It contains all code for analyses and figures presented in the manuscript.

## Preparation
To reproduce the analysis, you need a local installation of git (clone the repository with: `git clone git@github.com:SushiLab/reef-microbiomics-paper.git`) and R. Download the data from ==[Insert Zenodo URL]== into a directory named ``data`` at the same level as ``code``. As computational resources might be limited, we provide intermediary analyses files on Zenodo, too.

---

## Scripts
```
00    --> generate colour palette
01-15 --> curate data, including decontamination, merging, and rarefying
16-18 --> compute diversity metrics
19-21 --> community-based analyses: compare communities
22-23 --> community-based analyses: framework microbial exchange
24-27 --> community-based analyses: interaction networks
28-31 --> member-based analyses: interaction networks
32-35 --> member-based analyses: differential abundance analysis
36-40--> functional prediction analyses
```

#### ``00_Palettes.R``
1. Define colour palette (``colour_palette.RDS``)

#### ``01_Merge-ASV-Tables_Host-Associated.R``
1. Remove reads (i) unclassified at the domain-level, (ii) belonging to the domain Eukaryota, (iii) belonging to the family Mitochondria, (iv) belonging to the order Chloroplast from the host-associated data
2. Merge the two ASV tables using the ``mergeSequenceTables()`` and ``collapseNoMismatch()`` functions by dada2
3. Rename ASVs and sample names
4. Generate the ASV data table (``asv_dat_original_host-associated.tsv``) and ASV abundance table (``asv_abtab_original_host-associated.tsv``)

#### ``02_Explore-Original-Data_Host-Associated.R``
1. Compute and plot sequencing depth
2. Subsample to 1000 reads
3. Compare microbial communities (square root transformed Bray-Curtis distances)
4. Analyse sample-to-blank distances
5. Compare mucus vs control communities

#### ``03_Explore-Decontam-Thresholds_Host-Associated.R``
1. Run ``decontam()`` with thresholds 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5 and generate decontaminated ASV tables and stats files
2. Explore decontam scores
3. Compute and plot percentage of original reads/ASVs kept after ``decontam()`` using various thresholds (=> default threshold of 0.1 works well)

#### ``04_Explore-Decontaminated-Data_Host-Associated.R``
1. Subsample to 1000 reads
2. Analyse sample-to-blank distances for both the original and the decontaminated (at 0.1) ASV abundance tables
3. Compare microbial communities (square root transformed Bray-Curtis distances) after decontamination
4. Run PERMANOVA ($vegdist(\sqrt{asv\_abtab}) \sim species$) and PERMDISP (=> ``species`` significant)
5. Compare microbial communities (square root transformed Bray-Curtis distances) by taxonomic group

#### ``05_Generate-Cleaned-Tables_Host-Associated.R``
1. Remove all samples that were lost during the experiment as well as all positive and negative controls
2. Write cleaned metadata (``metadat_cleaned_host-associated.tsv``), ASV data table (``asv_dat_cleaned_host-associated.tsv``), and ASV abundance table (``asv_abtab_cleaned_host-associated.tsv``; threshold = 0.1)
3. Compute and plot sequencing depth

#### ``06_Merge-ASV-Tables_Free-living.R``
1. Remove reads (i) unclassified at the domain-level, (ii) belonging to the domain Eukaryota, (iii) belonging to the family Mitochondria, (iv) belonging to the order Chloroplast from the free-living and exuded data
2. Merge the two ASV tables using the ``mergeSequenceTables()`` and ``collapseNoMismatch()`` functions by dada2
3. Rename ASVs and sample names
4. Generate the ASV data table (``asv_dat_original_free-living.tsv``) and ASV abundance table (``asv_abtab_original_free-living.tsv``)

#### ``07_Explore-Original-Data_Free-living.R``
1. Compute and plot sequencing depth
2. Subsample to 1000 reads
3. Compare microbial communities (square root transformed Bray-Curtis distances)

#### ``08_Explore-Decontam-Thresholds_Free-living.R``
1. Run ``decontam()`` with thresholds 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5 and generate decontaminated ASV tables and stats files
2. Explore decontam scores
3. Compute and plot percentage of original reads/ASVs kept after ``decontam()`` using various thresholds (=> more stringent threshold of 0.2 works well)

#### ``09_Explore-Decontaminated-Data_Free-living.R``
1. Subsample to 1000 reads
2. Analyse sample-to-blank distances for both the original and the decontaminated (at 0.2) ASV abundance tables
3. Compare microbial communities (square root transformed Bray-Curtis distances) after decontamination
4. Run PERMANOVA ($vegdist(\sqrt{asv\_abtab}) \sim species$) and PERMDISP (=> ``species`` significant)
5. Compare microbial communities (square root transformed Bray-Curtis distances) by taxonomic group

#### ``10_Generate-Cleaned-Tables_Free-living.R``
1. Remove all samples that were lost during the experiment as well as all positive and negative controls
2. Write cleaned metadata (``metadat_cleaned_free-living.tsv``), ASV data table (``asv_dat_cleaned_free-living.tsv``), and ASV abundance table (``asv_abtab_cleaned_free-living.tsv``; threshold = 0.2)
3. Compute and plot sequencing depth

#### ``11_Merge-ASV-Tables.R``
1. Merge the decontaminated (host-associated and free-living) ASV tables using the ``mergeSequenceTables()`` and ``collapseNoMismatch()`` functions by dada2
2. Remove singleton ASVs (present in a single sample at an abundance of 1)
3. Rename ASVs and sample names
4. Generate the ASV data table (``asv_dat_original.tsv``) and ASV abundance table (``asv_abtab_original.tsv``)

#### ``12_Merge-Metadata_Seq-Depth.R``
1. Merge the two metadata tables (``metadat_cleaned.tsv``)
2. Compute and plot sequencing depth

#### ``13_Rarefy.R``
1. Rarefy 50 times to 2500 reads per sample (``asv_abtab_rarefied_{1,50}.tsv``)
2. Write metadata excluding dropped samples (``metadat_rarefied.tsv``)

#### ``14_Compile-Sequence-Files.R``
1. Load all 50 rarefied ASV abundance tables
2. Extract the sequence information
3. Write fasta files (``asv_abtab_rarefied_{1,50}.fasta``)

#### ``15_Annotate-Taxonomy.R``
1. Assign taxonomic info to differentially abundant ASVs using ``silva_nr99_v138.1_train_set.fa``
2. Write taxonomy file (``asv_dat_taxinfo.tsv``)

#### ``16_Compute-BC-Dissimilarity.R``
1. Load all 50 rarefied ASV abundance tables and the associated metadata
2. Compute the square-root transformed Bray-Curtis distances for each of the 50 ASV tables and save as a list
3. Set the diagonal to NA
4. Reduce the list to a single distance table by taking the sum of the corresponding elements and calculating the mean by dividing by the length of the dataframe
5. Write the square-root transformed Bray-Curtis distance table (``asv_bctab.tsv``)

#### ``17_Compute-Richness.R``
1. Load all 50 rarefied ASV abundance tables and the associated metadata
2. Compute the hill number diversity indices of each sample for each of the 50 ASV tables and save as a list
3. Reduce the list to a single distance table by taking the sum of the corresponding elements and calculating the mean by dividing by the length of the dataframe
4. Write the richness table (``asv_richtab.tsv``)

#### ``18_Compute-Shared.R``
1. Load all 50 rarefied ASV abundance tables and the associated metadata
2. Compute the $designdist(asv\_abtab,method="J",terms="binary")$ for each of the 50 ASV tables and save as a list
3. Set the diagonal to NA
4. Reduce the list to a single distance table by taking the sum of the corresponding elements and calculating the mean by dividing by the length of the dataframe
5. Write the number of shared ASVs table (``asv_shared.tsv``)

#### ``19_BDiv_Compare-Communities.R``
**Compare the host-associated and exuded microbial communities by species/sample type**
1. Load the square-root transformed Bray-Curtis distance table
2. Perform a PCoA on the distances
3. Run a PERMANOVA to test for differences explained by species ($asv\_bctab \sim species$) and PERMDISP (=> ``species`` significant)

#### ``20_ADiv_Host-Exudates.R``
**Correlate alpha-diversity changes of exuded vs host-associated microbial communities**
1. Load the richness table
2. Plot exudate vs host Hill numbers by species and get correlation coefficient

#### ``21_BDiv_Complexity_Host-Exudates.R``
**Compare the host-associated and exuded microbial communities of hosts living in complex vs reduced reef communities**
1. Select samples from phase 2
2. For each source and species, run a PERMANOVA to test for differences explained by the ambient community complexity ($asv\_bctab \sim treatment\_binary$)

#### ``22_BDiv_Specificity-Exudation.R``
**Explore host-specificity and exudation**
1. Transform the distance table into a 3-column dataframe (keeping all data from all three phases)
2. Prepare the host-to-tank distance dataframe
3. Prepare the exudates-to-incubation control dataframe
4. Plot exudate to control (x-axis) vs host to tank (y-axis)
5. Draw an ellipse at 0.5 (half of the points lie within the ellipse)
6. Use data to define cutoffs for groups (mean)
7. Test whether the groups are significantly different ($manova(cbind(bc\_exudatetocontrol,bc\_hosttotank) \sim species)$) (=> ``species`` significant)
8. Add standard deviation to each data point (as one point is the mean of three distances)

#### ``23_BDiv_Complexity_Specificity-Exudation.R``
**Explore host-specificity and exudation as a function of the complexity of the reef community**
1. Select samples from phase 2
2. Transform the distance table into a 3-column dataframe
3. Prepare the host-to-tank distance dataframe
4. Prepare the exudates-to-incubation control dataframe
5. Compare both distances of hosts living in a complex vs reduced reef community by species with ``p.adjust.method="holm"``
6. Test whether the mean of the groups is significantly different as a response to treatment ($compare\_means(bc\_dissim \sim treatment\_binary)$)

#### ``24_BDiv_Similarity-Between-Hosts.R``
1. Select species samples from phase 2
2. Extract distances from ``asv_bctab.tsv`` across treatment and source
3. Group by host and average distance across samples
4. Convert dissimilarity to similarity
5. Write the edge files
6. Compute stats on the full data comparing treatments and write the stats files

#### ``25_BDiv_Richness-By-Host.R``
1. Select species samples from phase 2
2. Extract richness from ``asv_richtab.tsv`` across treatment and source
3. Group by host and average richness across samples
4. Write the node files
5. Compute stats on the full data comparing treatments and write the stats files

#### ``26_BDiv_Similarity-Richness-Network.R``
1. Load edge and node files
2. Compute similarity ratio between complex and reduced by calculating $0.5-\frac{similarity_{reduced}}{similarity_{complex}+similarity_{reduced}}$ (=> negative values mean the similarity between hosts is higher in a reduced habitat, positive values in a complex habitat; the higher the values the more the similarities between the two treatments differ)
3. Compute richness ratio between complex and reduced by calculating $0.5-\frac{richness_{reduced}}{richness_{complex}+richness_{reduced}}$ (=> negative values mean the richness of a host is higher in a reduced habitat, positive values in a complex habitat; the higher the values the more the richness between the two treatments differ)
4. Initiate the network plot using the packages ``igraph`` and ``ggraph``
5. Plot the network with colours denoting in which treatment the values are higher and line width/circle size indicating the magnitude of the difference

#### ``27_BDiv_Similarity-Richness-Network_Supplementary.R``
1. Load edge and node files
2. Initiate the network plot using the packages ``igraph`` and ``ggraph``
3. Plot the network with line width and colours denoting the magnitude of the compositional interactions and circle size indicating the richness

#### ``28_ASVs_Shared-Between-Hosts.R``
1. Select species samples from phase 2
2. Extract distances from ``asv_shared.tsv`` across treatment and source
3. Group by host and average distance across samples
4. Write the edge files
5. Compute stats on the full data comparing treatments and write the stats files

#### ``29_ASVs_Unique-By-Host.R``
1. Load all 50 rarefied ASV abundance tables and the associated metadata
2. Select species samples from phase 2
3. Group by host and compute number of unique ASVs across treatment and source using the function ``upset_data()``
4. Reduce the list to a single distance table by taking the sum of the corresponding elements and calculating the mean by dividing by the length of the dataframe
5. Write the node files

#### ``30_ASVs_Shared-Unique-Network.R``
1. Load edge and node files
2. Compute shared ratio between complex and reduced by calculating $0.5-\frac{shared_{reduced}}{shared_{complex}+shared_{reduced}}$ (=> negative values mean the number of shared ASVs between hosts is higher in a reduced habitat, positive values in a complex habitat; the higher the values the more the number of shared ASVs between the two treatments differ)
3. Compute unique ratio between complex and reduced by calculating $0.5-\frac{unique_{reduced}}{unique_{complex}+unique_{reduced}}$ (=> negative values mean the number of unique ASVs of a host is higher in a reduced habitat, positive values in a complex habitat; the higher the values the more the number of unique ASVs between the two treatments differ)
4. Initiate the network plot using the packages ``igraph`` and ``ggraph``
5. Plot the network with colours denoting in which treatment the values are higher and line width/circle size indicating the magnitude of the difference

#### ``31_ASVs_Shared-Unique-Network_Supplementary.R``
1. Load edge and node files
2. Initiate the network plot using the packages ``igraph`` and ``ggraph``
3. Plot the network with line width and colours denoting the level of shared ASVs and circle size indicating the unique ASVs

#### ``32_ASVs_DA.R``
**Analyse the differentially abundant ASVs between source and treatment in two separate analyses**
1. Select samples from phase 2
2. Load the modified ancombc R files
3. For each source separately, subset and transform the data by species
4. Construct TreeSummarizedExperiment object
5. Use ``ancombc2()`` with fixed effect $treatment\_binary$ and ``p_adj_method="holm"``
6. Write ANCOMBC2 output files (``diffab_species.tsv``)

#### ``33_ASVs_DA_Volcano-Plot.R``
**Plot the DA output**
1. Plot log-fold changes (based on the natural log) with significant differences in red (after correcting for multiple comparisons using the Holm method)

#### ``34_ASVs_DA_iTOL_LFC.R``
**Prepare the compilated log-fold change file for iTOL**
1. Merge the DA output
2. Pull the significant log-fold changes between complex and reduced communities (``asv_diffab_lfc.tsv``)

#### ``35_ASVs_DA_iTOL_Tree.R``
**Examine the phylogeny of the differentially abundant ASVs (across host)**
1. Align sequences through multiple sequence aligner using Muscle
2. Test maximum-likelihood models
3. Select the two superior ones and compute trees
4. Bootstrap 100 times

#### ``36_FunPred_KEGG-Download.R``
1. Download the Pathway list, the KO list, and the KO to Pathway link files
2. Reformat and write files (``pathway_list.tsv``, ``ko_list.tsv``, ``ko_pathway_link.tsv``)

#### ``picrust_rarefied_data.sh``
```
#! /bin/bash
#SBATCH --job-name picrust      # how the job will be named
#SBATCH --output picrust.log    # filename for the output messages
#SBATCH --error picrust.log     # filename for the error messages
#SBATCH -n 20                   #​ number of CPUs
#SBATCH --mem-per-cpu=5G        #​ memory per CPU
#SBATCH --time=14-23:00:00      # maximum time to complete (after which it is killed)

eval "$(conda shell.bash hook)"

conda activate picrust2_env

mkdir /[INSERT PATH]/data/picrust_output

for i in /[INSERT PATH]/data/picrust_input/*.fasta ; \
do j=$( echo $i | sed 's/.fasta//' | sed 's/.*picrust_input[/]//' ) ; \
picrust2_pipeline.py -s "$i" -i /[INSERT PATH]/data/picrust_input/"$j".tsv \
-o /[INSERT PATH]/data/picrust_output/"$j" -p 20 ; done
```

#### ``37_FunPred_Compute-KO-per-Sample-Abundances.R``
1. Load all 50 generated predicted unstratified metagenome files (``pred_metagenome_unstrat.tsv.gz``)
2. Get a list of the functions (K numbers) that are present across all 50 rarefaction files
3. Extract these functions from each abundance file and save as a list
4. Reduce the list to a single abundance table by taking the sum of the corresponding elements and calculating the mean by dividing by the length of the dataframe
5. Write the abundance table (``ko_abtab_sample.tsv``)

#### ``picrust_pathway_to_ko.sh``
```
#! /bin/bash
#SBATCH --job-name picrust      # how the job will be named
#SBATCH --output picrust.log    # filename for the output messages
#SBATCH --error picrust.log     # filename for the error messages
#SBATCH -n 20                   #​ number of CPUs
#SBATCH --mem-per-cpu=5G        #​ memory per CPU
#SBATCH --time=14-23:00:00      # maximum time to complete (after which it is killed)

eval "$(conda shell.bash hook)"

conda activate picrust2_env

pathway_pipeline.py -i /[INSERT PATH]/data/ko_abtab_sample.tsv \
-o /[INSERT PATH]/data/picrust_output \
--no_regroup --map /[INSERT PATH]/data/KEGG_pathway_to_ko.tsv -p 20
```

#### ``38_FunPred_Pathway_Species_.R``
1. Load the KEGG pathway list (``pathway_list_cat.csv``) with hand annotated categories
2. Load the KEGG pathway converted abundance table (``path_abun_unstrat.tsv.gz``) and extract pathways associated with metabolism, but remove global and overview maps
3. Compute the square-root transformed Bray-Curtis distances
4. Perform a PCoA on the distances
5. Run a PERMANOVA to test for differences explained by species ($asv\_bctab \sim species$) and PERMDISP (=> ``species`` significant) and test pairwise using `pairwise.adonis2()`

#### ``39_FunPred_Compute-KO-per-ASV-Abundances.R``
1. Load all 50 generated predicted KO abundance files (``KO_predicted.tsv.gz``)
2. Get a list of the KOs that are present across all 50 rarefaction files
3. Extract these KOs from each abundance file and save as a list
4. Reduce the list to a single abundance table by taking the sum of the corresponding elements and calculating the mean by dividing by the length of the dataframe
5. Write the abundance table (``ko_abtab_asv.tsv``)

#### ``40_FunPred_Upset_Sinularia.R``
1. Load the KEGG pathway list (``pathway_list_cat.csv``) with hand annotated categories
2. Load the ko-to-pathway link file (``ko_pathway_link.tsv``) and extract pathways associated with metabolism, but remove global and overview maps
3. Load the KO abundance table (``ko_abtab_asv.tsv``)
4. Select the ASVs of interest (two closely related members within the family Flavobacteriaceae that are differentially abundant as a response to treatment in *Sinularia* sp.)
5. Summarise the KOs to pathways
6. Remove pathways that are absent
7. Convert pathway abundance table to presence-absence and annotate
8. Generate upset plot

---

## Metadata
==Name file and provide link==

- ``id`` (ID): unique sample identifier (sample name)
- ``sample`` (categorical): medium sampled
  - levels host-associated: biofilm; extraction blank; extraction mock; mucus; PBS control; PCR blank; PCR mock; tissue; Tris-NaCl control
  - levels free-living: extraction blank; extraction mock; exudates; incubation control; microscopy-slide control; PCR blank; PCR mock; tank water
- ``type`` (categorical): whether it is a true sample or a control
  - levels: control; sample
- ``fragment_id`` (categorical): fragment ID; all fragments (clones) of a species were generated from a single mother colony; each fragment underwent a predetermined set of three treatments; while each fragment existed only once, this ID is not unique in the metadata as fragments were sampled three times (once per phase)
  - naming: characters 1&ndash;3 = species name; 4&ndash;9 = treatments in chronological order; 9 = replicate number
- ``species`` (taxa): taxonomic name of sampled species to deepest taxonomic level known
  - levels: Caulerpa sp.; Dictyoceratida; Haliclona cnidata; Montipora digitata; Peyssonnelia sp.; Porites verrucosa; Sinularia sp.; Xenia sp.
- ``taxonomic_group`` (categorical): taxonomic group to which sampled species belongs
  - levels: macroalgae; soft coral; sponge; stony coral
- ``phase`` (categorical): phase during which the sample was collected
  - levels: P0; P1; P2
- ``treatment`` (categorical): treatment the fragment was in at the time of sampling
  - levels: monoculture; polyculture controlled; polyculture seminatural
- ``treatment_binary`` (categorical): binary treatment classification
  - levels: biodiverse; degraded
- ``replicate`` (replicate): replicate number; each set of treatments was undergone by three fragments
- ``tank`` (categorical): tank the fragment was kept in at the time of sampling; corresponding to the official tank name given by the Ocean2100 facility
  - levels: A3; A5; B1; B3; B6; C2; C3; C4; C6; Ex
- ``incubation_spot`` (categorical): spot in incubation chamber on which the incubation took place (only free-living)
- ``incubation_time_start`` (time): time fragment was placed into incubation jar (only free-living)
- ``incubation_time_end`` (time): time incubation jar was removed from incubation spot (only free-living)
- ``sampling_date`` (date): date of sample collection
- ``sampling_time`` (time): time of sample collection
- ``filter_type`` (categorical): type of Sterivex filter used to collect sample (only free-living)
  - levels: SVGPL10RC; SVGV010RS
- ``weight_brutto`` (numeric trait): buoyant or wet weight of the fragment with hook, screw nut, silk, microscopy slide, label etc; buoyant weight measured for stony corals; wet weight measured for soft corals, macroalgae, and sponges (only free-living)
- ``weight_netto`` (numeric trait): buoyant or wet weight calculated from weight_brutto minus weight_tara (hook, screw nut, silk, microscopy slide, label etc) (only free-living)
- ``extraction_date`` (date): date of DNA extraction
- ``extraction_round`` (numeric): round during which sample was extracted
- ``plate_number`` (numeric): number of 96-well plate (1&ndash;4)
- ``plate_position`` (alphanumeric): position of sample on 96-well plate (rows: A&ndash;H; columns: 01&ndash;12)
- ``plate_position_sequential`` (numeric): unique numeric position of sample on 96-well plate (plate 1: 1&ndash;96; plate 2: 97&ndash;192; etc)
- ``biomass_extraction_mg`` (numeric trait): biomass used for homogenisation and biofilm enrichment, respectively (sponge and macroalgae); coral-mucus biomass was not weighed (instead, 100 $\mu$l of mucus was used) (only host-associated)
- ``gDNA_conc_ng_ul_QuantiFluor`` (numeric trait): genomic DNA concentration, measured with the Promega QuantiFluor dsDNA System (E2670); the accuracy at these low levels is not great and is simply indicative of an expected concentration range (only host-associated)
- ``gDNA_conc_ng_ul`` (numeric trait): genomic DNA concentration, measured with the Qubit (1X HS kit)
- ``primers`` (categorical): common name of the used forward and reverse primers (main reference for forward primer: Parada, A. E., Needham, D. M., & Fuhrman, J. A. (2016). Every base matters: assessing small subunit rRNA primers for marine microbiomes with mock communities, time series and global field samples. *Environmental microbiology*, *18*(5), 1403-1414. [doi.org/10.1111/1462-2920.13023](https://doi.org/10.1111/1462-2920.13023); for reverse primer: Apprill, A., McNally, S., Parsons, R., & Weber, L. (2015). Minor revision to V4 region SSU rRNA 806R gene primer greatly increases detection of SAR11 bacterioplankton. *Aquatic Microbial Ecology*, *75*(2), 129-137. [doi.org/10.3354/ame01753](https://doi.org/10.3354/ame01753))
  - levels: 515F-Y/806R_B
- ``primer_fwd_seq`` (categorical): nucleotide sequence of the 16S forward primer used (5' to 3' and without Nextera adapter sequence)
  - levels: GTGYCAGCMGCCGCGGTAA
- ``primer_rev_seq`` (categorical): nucleotide sequence of the 16S reverse primer used (5' to 3' and without Nextera adapter sequence)
  - levels: GGACTACNVGGGTWTCTAAT
- ``barcodes`` (id): unique barcode identifier (for both i5 and i7 sequences); these are unique dual index adapters from Illumina
- ``barcode_rev_i7_seq`` (id): nucleotide sequence for the i7 index as the reverse complement (used for sample sheet)
- ``barcode_fwd_i5_seq`` (id): nucleotide sequence of the i5 index (used for sample sheet)
- ``aDNA_conc_ng_ul`` (numeric trait): concentration of the PCR product after clean-up, measured with the Nanodrop Eight
- ``aDNA_A260_A280`` (numeric trait): purity of the PCR product after clean-up, measured with the Nanodrop Eight; these measurements are particularly inaccurate for samples with low DNA concentration
- ``aDNA_A260_A230`` (numeric trait): purity of the PCR product after clean-up, measured with the Nanodrop Eight; these measurements are particularly inaccurate for samples with low DNA concentration
- ``aDNA_smear_conc_nmol_l`` (numeric trait): smear concentration of the region of interest as measured with the Fragment Analyser (NGS SS kit); smear range from 400&ndash;500 bases  
- ``pooling_conc`` (numeric trait): concentration in nmol/l at which sample was added to the library pool
- ``sequencing_date`` (date): date of sequencing
- ``sequencing_run`` (numeric): run during which sample was sequenced
- ``notes`` (comment): any additional notes collected throughout the experiment
