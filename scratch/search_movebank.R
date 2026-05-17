library(move2)

# Movebank allows searching for studies
cat("Searching Movebank for public studies with 'dual', 'double', or 'light' and 'gps' in the name...\n")

# This requires no account for public data.
studies <- movebank_download_study_info()

# Filter for public studies
public_studies <- studies[studies$public == TRUE | studies$i_am_owner == TRUE, ]

# Search for relevant keywords
keywords <- c("light", "gls", "dual", "double", "gps", "argos")
# Actually, just search for known species
marine_studies <- public_studies[grep("albatross|seal|turtle", public_studies$name, ignore.case = TRUE), ]

cat(sprintf("Found %d public marine studies.\n", nrow(marine_studies)))
head(marine_studies$name, 20)

# We can pick a specific one, or just save the list
saveRDS(marine_studies, "scratch/movebank_public_marine_studies.rds")
cat("Saved study list to scratch/movebank_public_marine_studies.rds\n")
