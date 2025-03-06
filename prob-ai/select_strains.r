
## https://googlesheets4.tidyverse.org/
# install.packages("googlesheets4")

library("dplyr")
library("tidyverse")
library("googlesheets4")

####################################################################
writeLines(" - configuring access to Google Drive")
apik = "AIzaSyBbVOwEz8_91ZeQKzdOfEFAx8ohDV83yjE"
gs4_auth_configure(api_key = apik)

####################################################################
writeLines(" - reading the probiotic spreadsheet from Google Drive")

gs4_deauth()
url = "https://docs.google.com/spreadsheets/d/18B0et5eFGacukmWygwskXvBha-YnOw_tG4YAlIaYHkc/edit?gid=312024716#gid=312024716"
print(paste("The file to be read is:", url))

probs = read_sheet(url)
print(paste("N. of records read from file:", nrow(probs)))

####################################################################
writeLines(" - data preprocessing ")
probs <- probs |> 
  rename(probiotic = `Confirmed probiotic (Y/N) (in vivo trials)`, ncbi_genome = `NCBI genome (Y/N)`, target = `target (plant vs animal)`) |>
  mutate(probiotic = tolower(probiotic), ncbi_genome = tolower(ncbi_genome))

####################################################################
writeLines(" - data statistics (EDA)")
probs |>
  group_by(probiotic, ncbi_genome) |>
  summarise(N = n()) |>
  spread(key = probiotic, value = N) |>
  print()

####################################################################
writeLines(" - selecting probiotic strains to download")
probs <- probs |>
  select(Organism.Name, Strain, BioSample, BioProject, Assembly, target, probiotic, ncbi_genome) |>
  filter(probiotic == "yes", ncbi_genome == "yes")

print(paste("N. of selected probiotic strains:", nrow(probs)))

print("DONE!!")