# cd4_isoform_counts  <- readRDS("./results/input_data/cd4_isoform_counts.rds")
# isoform_fdata <- readRDS("./results/input_data/isoform_fdata.rds")
# aw_meta_cv2 <- readRDS("./results/input_data/aw_meta_cv2.rds")
# 
# cd4_isoform_counts <- cd4_isoform_counts[
#   ,colnames(cd4_isoform_counts) %in% rownames(aw_meta_cv2)]
# 
# cd4_meta <- aw_meta_cv2[rownames(aw_meta_cv2) %in% colnames(cd4_isoform_counts),]
# 
# cd4_meta <- cd4_meta[match(colnames(cd4_isoform_counts),rownames(cd4_meta)),]
# 
# identical(colnames(cd4_isoform_counts),rownames(cd4_meta))
# identical(rownames(cd4_isoform_counts),rownames(isoform_fdata))
# 
# se <- SummarizedExperiment(
#   assays = SimpleList(
#     counts=as.matrix(cd4_isoform_counts)
#   ),
#   colData = DataFrame(cd4_meta),
#   rowData = isoform_fdata
# )


### Overflow

#### Test IsoformSwitchAnalyzer


# Fix issue with isoform mismatch
# gtf <- import(gtf_path)
# 
# gtf_filt <- gtf[gtf$transcript_id %in% counts_ids]
# 
# export(gtf_filt, "gtf_matched_to_quant.gtf")

# gtf_fixed <- rtracklayer::import("../data/gtf_matched_to_quant.gtf")
# 
# cd4_se_filt <- cd4_se[rownames(cd4_se) %in% gtf_fixed$transcript_id, ]
# cd16_se_filt <- cd16_se[rownames(cd16_se) %in% gtf_fixed$transcript_id, ]



# cd4_iso_res <- run_isoform_switching(
#     se = cd4_se_filt, 
#     assay="counts", 
#     condition_col = "pftfev1fvc_actual",
#     covariates=c("gli_age", "gli_sex", "gli_height", "fis_bin"),
#     condition_levels = NULL,
#     gtf_path = "../data/gtf_matched_to_quant.gtf",
#     genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
#     verbose = TRUE,
#     removeNonConvensionalChr = TRUE,
#     ignoreAfterPeriod    = TRUE
# )
# 
# cd16_iso_res <- run_isoform_switching(
#     se = cd16_se_filt, 
#     assay="counts", 
#     condition_col = "pftfev1fvc_actual",
#     covariates=c("gli_age", "gli_sex", "gli_height", "fis"),
#     condition_levels = NULL,
#     gtf_path = "../data/gtf_matched_to_quant.gtf",
#     genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
#     verbose = TRUE,
#     removeNonConvensionalChr = TRUE,
#     ignoreAfterPeriod    = TRUE
# )


#### Test for DEXSeq

# se <- tidyexposomics:::.update_assay_colData(fev1_fvc_expom,"CD4 T cell Isoforms")
# 
# SummarizedExperiment::assay(se,"counts") <- SummarizedExperiment::assay(se,"counts") |> 
#   as.data.frame() |> 
#   mutate_all(as.integer) |> 
#   as.matrix()
# 
# dex_out <- run_dexseq_parallel(
#     se = se,
#     condition_col = "pftfev1fvc_actual",
#     covariates = c("gli_age", "gli_sex", "gli_height"),
#     gtf_path = "../../asd_transcriptional_atlas/data/Homo_sapiens.GRCh38.111.gtf/Homo_sapiens.GRCh38.111.gtf",
#     cores = 12,
#     backend = "Multicore",
#     verbose = TRUE
# )

#### Clean the GTF


# # load the gtf
# gtf <- rtracklayer::import("../../asd_transcriptional_atlas/data/Homo_sapiens.GRCh38.111.gtf/Homo_sapiens.GRCh38.111.gtf")
# 
# tx_annot <- rowData(isoform_se) |> as.data.frame()
# 
# gtf_mapped <- as.data.frame(gtf) %>%
#   inner_join(tx_annot %>%
#                dplyr::select(tx),
#              by=c("transcript_id"="tx"))
# 
# gtf_mapped$isoform_id <- gtf_mapped$transcript_id
# 
# gc()
# 
# #write_tsv(gtf_mapped,file="../results/gtf_mapped.gtf")
# 
# # gtf_mapped <- gtf_mapped %>%
# #   filter(type=="exon") %>%
# #   mutate(ref_gene_id=gene_id) %>%
# #   dplyr::select(-type)
# # #mutate(type=as.character(type))
# 
# gtf_granges <- GRanges(gtf_mapped)
# 
# gc()
# rtracklayer::export(gtf_granges,"./data/Homo_sapiens.GRCh38.111.gtf/gtf_mapped.gtf")
