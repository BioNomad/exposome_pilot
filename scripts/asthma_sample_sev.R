# asthma_sample_sev <- aw_meta_cv2 |> 
#   dplyr::select(
#     asuiscore,
#     rsuiscore,
#     actscore,
#     pftfev1fvc_actual
#   ) |> 
#   rownames_to_column("sample") |> 
#   mutate(
#     asuiscore_prank = 1 - percent_rank(asuiscore),
#     rsuiscore_prank = 1 - percent_rank(rsuiscore),
#     actscore_prank = 1 - percent_rank(actscore),
#     pftfev1fvc_prank = 1 - percent_rank(pftfev1fvc_actual),
#     severity_sum = asuiscore_prank + rsuiscore_prank + actscore_prank + pftfev1fvc_prank
#   ) |> 
#   # filter out sample outliers in omics data
#   filter(!sample %in% c("s937", "s3331","s3692")) |> 
#   filter(!is.na(severity_sum)) |> 
#   arrange(desc(severity_sum))
# 
# write.csv(asthma_sample_sev,
#           file="./results/asthma_sample_sev.csv",
#           quote = F)
