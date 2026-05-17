# fetch_movebank.R
# Must be run after 'move2' is installed
library(move2)

cat("Querying Movebank API for public study metadata...\n")

# Try to download public study info without credentials
study_info <- tryCatch({
  movebank_download_study_info()
}, error = function(e) {
  cat("Error downloading study info:", e$message, "\n")
  return(NULL)
})

if (!is.null(study_info)) {
  # Filter to studies that are public or where the user is an owner
  public_studies <- study_info[study_info$public == TRUE, ]
  
  # Search for dual-tag studies by looking at the sensor_type_ids column.
  # GPS = 653, Light-level = 673, Argos = 82
  # We want a study that has both Light (673) and GPS (653) or Fastloc/Argos
  has_light <- sapply(public_studies$sensor_type_ids, function(x) {
    if (is.null(x)) return(FALSE)
    673 %in% unlist(strsplit(as.character(x), ","))
  })
  
  has_gps <- sapply(public_studies$sensor_type_ids, function(x) {
    if (is.null(x)) return(FALSE)
    653 %in% unlist(strsplit(as.character(x), ","))
  })
  
  dual_tag_studies <- public_studies[has_light & has_gps, ]
  
  if (nrow(dual_tag_studies) > 0) {
    cat("\nFound", nrow(dual_tag_studies), "public dual-tag studies (Light + GPS)!\n")
    
    # Sort by number of locations (descending)
    dual_tag_studies <- dual_tag_studies[order(dual_tag_studies$number_of_events, decreasing = TRUE), ]
    
    print(head(dual_tag_studies[, c("id", "name", "number_of_events", "sensor_type_ids")], 10))
    
    # Save the best candidate study ID to a text file so we know what to download
    best_id <- dual_tag_studies$id[1]
    best_name <- dual_tag_studies$name[1]
    
    cat("\nTop candidate study:", best_name, "(ID:", best_id, ")\n")
    writeLines(as.character(best_id), "scratch/movebank_target_id.txt")
  } else {
    cat("\nCould not find any public studies with BOTH Light (673) and GPS (653) sensors.\n")
    # Let's search just by name then
    keyword_match <- grep("light|gls|dual", public_studies$name, ignore.case = TRUE)
    cat("But found", length(keyword_match), "studies with 'light' or 'gls' in the title.\n")
    print(head(public_studies[keyword_match, c("id", "name")]))
  }
}
