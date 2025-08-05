# Matt Prill
# MALDI-TOF MS Msc
# Merging the Data
# Here I merge the 2 data sources and harmonise them.



# Libraries ----
library(tidyverse)
library(reticulate)
# ----
# ----
# ----
# Merging the Metadata ----
# CC22
cc22_meta <- read_csv("data/CC22/CC22_phenotype.csv") %>% 
  mutate(source = "cc22") %>%  # To distinguish between datasets once merged
  mutate(
    Age = as.numeric(Age),
    across(-c(Age, toxicity, biofilm, CCI), as.factor)  # Set other cols as factor
  )

# Cork
cork_meta <- read_csv("data/Cork/metadata.csv") %>%
  mutate(source = "cork") %>%
  rename(   # Standardise columns and conform to CC22
    Infection.Site.Bacteraemia = `Infection Site Bacteraemia`,                        # Standardise columns and conform to CC22
    Infection.Site.SSTI = `Infection Site SSTI`,
    Infection.Site.Septic.Arthiritis = `Infection Site Septic Arthiritis`,            # Standardise columns and conform to CC22
    Infection.Site.Osteomyelitis = `Infection Site Osteomyelitis`,
    Infection.Site.Infective.Endocarditis = `Infection Site Infective Endocarditis`,  # Standardise columns and conform to CC22
    Infection.Site.Prosthesis = `Infection Site Prosthesis`,
    Prosthesis.type = `Prosthesis type`,                                              # Standardise columns and conform to CC22
    Infection.Site.Abscess = `Infection Site Abscess`,
    Type.of.abscess = `Type of abscess`,                                              # Standardise columns and conform to CC22
    Collection = `MSSA/MRSA`,
    Outcome30days = `Mortality_30_Days`,                                              # Standardise columns and conform to CC22
    Outcome90days = `Mortality_90_Days`
  ) %>%
  mutate(
    Age = as.numeric(Age),
    across(-Age, as.factor)  # Set other Cols at factor
  )

# Combine both
merged_meta <- bind_rows(cc22_meta, cork_meta)  # Now need to merge columns better e.g. alive/death and site of infection
#write.csv(merged_meta, "merged_meta.csv")



merged_meta <- merged_meta %>%
  mutate(
    Outcome30days = fct_collapse(Outcome30days,
                                 Alive = c("Alive", "No"),  # Sets mortality = no as 'Alive'
                                 Dead = c("Death", "Yes")   # Sets mortality = yes as 'Dead'
    ),
    Outcome90days = Outcome90days %>%
      fct_explicit_na(na_level = "Death") %>%  # Sets NAs as Death
      fct_collapse(
        Alive = c("Alive", "No"),  # Sets mortality = no as 'Alive'
        Dead = c("Death", "Yes")   # Sets mortality = yes as 'Dead'
      )
  )



# The order of the arguments enforces a hierarchy. For example,
# where 'Infection.Site.Infective.Endocarditis' and 'Infection.Site.Osteomyelitis'
# both == 'yes', the focusofinfection is kept as Infection.Site.Infective.Endocarditis

merged_meta <- merged_meta %>%
  mutate(
    FocusOfInfection = case_when(
      Infection.Site.Infective.Endocarditis == "Yes" ~ "InfectiveEndocarditis",
      Infection.Site.Osteomyelitis == "Yes" ~ "BoneJoint",
      Infection.Site.Septic.Arthiritis == "Yes" ~ "BoneJoint",
      Infection.Site.SSTI == "Yes" ~ "SoftTissueInfection",
      Infection.Site.Abscess == "Yes" ~ "DeepTissueAbscess",
      is.na(FocusOfInfection) & Prosthesis.type == "Line" ~ "Line",
      is.na(FocusOfInfection) & Prosthesis.type %in% c("Joint", "Orthopaedic Metalwork") ~ "BoneJoint",
      is.na(FocusOfInfection) & Prosthesis.type == "Cardiac Valve/ Device" ~ "CardiacDevice",
      is.na(FocusOfInfection) & Prosthesis.type == "Urinary Catheter" ~ "Urogenital",
      is.na(FocusOfInfection) & Prosthesis.type == "Other" ~ "ProstheticMaterial",
      is.na(FocusOfInfection) ~ "Unknown",
      TRUE ~ FocusOfInfection  # Keeps existing arguments the same
    ),
    FocusOfInfection = as.factor(FocusOfInfection)
  )


# Omit non SAB isolates & Remove redundant columns and
merged_meta <- merged_meta %>%
  filter(Infection.Site.Bacteraemia != "No" | is.na(Infection.Site.Bacteraemia)) %>% 
  select(sangerID:source) 
# REMEMBER TO REMOVE ROWS WITHOUT CORRESPONDING MASS SPECTRA




# write.csv(merged_meta, "data/merged_data/merged_meta.csv", row.names = FALSE)
# ----
# ----
# ----
# Final Meta ----
merged_meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")





# ----
# ----
# ----
# Merging Mass Spectra ----
# Remember to remove rows with no corresponding mass spectra and viceversa

# cc22 is an excel file
cc22_ms <- read.csv("data/CC22/CC22_relative_intensity.csv")
str(cc22_ms)  # first col = m/z, next are isolates, (already relative intensity?)

# cork is an R.data file
load("data/Cork/Cork_peaksDF.Rda")  # Check what object is loaded
str(peaksDFs)      # 181 dataframes, first column = m/z, second = intensity (absolute)



# Conforming cc22 to cork format
cc22_list <- lapply(names(cc22_ms)[-1], function(colname) {  # Separates all of the isolates columns into different dataframes
  data.frame(
    mass = cc22_ms$position,        # Takes the first column (m/z) and joins to each df
    intensity = cc22_ms[[colname]]  # Assigns each isolates intensities to intensity
  )
})
names(cc22_list) <- names(cc22_ms)[-1]  # Ensures that each object in the data is labelled with the isoalte name

# Combining the data into a single list
combined_spectra <- c(cc22_list, peaksDFs)
str(combined_spectra)




combined_spectra[[ "ASARM183" ]]  # Example for extracting particular isolate mass spectra
combined_spectra[[ "A030" ]]  # Example for extracting particular isolate mass spectra
max(combined_spectra[["ASARM183"]]$intensity, na.rm = TRUE)  # Checking max m/z 
combined_spectra[[ "ASARM183" ]]  # Example for extracting particular isolate mass spectra




# ----
# ----
# ----
# Filter spectra and metadata  for Matching Entries ----

# Retain the 'repeats' and standardise names to the meta
rename_map <- c(
  "A078 repeat" = "A078",
  "A095 #2" = "A095",
  "A115 Real" = "A115"
)

# Replace the non repeats (assumes theyre faulty)
for (old_name in names(rename_map)) {
  new_name <- rename_map[[old_name]]
  
  if (old_name %in% names(combined_spectra)) {
    combined_spectra[[new_name]] <- combined_spectra[[old_name]]
    combined_spectra[[old_name]] <- NULL
  }
}


# Step 1: Get names
spectra_names <- names(combined_spectra)
metadata_names <- merged_meta$sampleID

# Step 2: Find the intersection
matched_names <- intersect(spectra_names, metadata_names)

# Step 3: Filter both datasets
combined_spectra_filtered <- combined_spectra[matched_names]
combined_metadata_filtered <- merged_meta[merged_meta$sampleID %in% matched_names, ]


# Example PLot
plot(
  combined_spectra_filtered[[ "ASASM430" ]]$mass,
  combined_spectra_filtered[[ "ASASM430" ]]$intensity,
  type = "l",
  xlab = "m/z",
  ylab = "Intensity",
  xlim = c(2000,10000)
)


# ----
# ----
# ----
# Final Result ----
#save(combined_spectra_filtered, file = "data/merged_data/combined_spectra_filtered.RData")
#write.csv(combined_metadata_filtered, "data/merged_data/combined_metadata_filtered.csv")

# CHECK 
load("data/merged_data/combined_spectra_filtered.RData")  # called combined_spectra_filtered
str(combined_spectra_filtered)
length(combined_spectra_filtered)

combined_metadata_filtered <- read.csv("data/merged_data/combined_metadata_filtered.csv")
head(combined_metadata_filtered)
nrow(combined_metadata_filtered)




# Same again but with the raw CC22 ----

py_config()


py_run_string("
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
import pandas as pd

filedir = 'data/CC22/CC22_collection_rep4/'

position = []
heightRaw = []
sampleID = []

for filename in listdir(filedir):
    if filename.endswith('.txt'):
        with open(filedir + filename) as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if 'binary' in line:
                x = [float(val) for val in lines[i + 2].split()[2:]]
                y = [float(val) for val in lines[i + 5].split()[2:]]  # <— convert to float
                position.append(np.array(x))
                heightRaw.append(np.array(y))
                sampleID.append(filename[:-4].replace(' ', ''))
                break

position = np.stack(position, axis=0)
heightRaw = np.stack(heightRaw, axis=0)

# Create DataFrame with raw (non-relative) intensity, as float
df = pd.DataFrame(position[0], columns=['position'])
for i, sid in enumerate(sampleID):
    df[sid] = heightRaw[i]

df.to_csv('CC22_rawSpectra.csv', index=False)
")



cc22_raw <- read.csv("CC22_rawSpectra.csv")  # raw data
str(cc22_raw) ## same as before (cc22_ms) expect not relative or corrected

# cork is an R.data file
load("data/Cork/Cork_peaksDF.Rda")  # Check what object is loaded
str(peaksDFs)      # 181 dataframes, first column = m/z, second = intensity (absolute)


merged_meta <- read.csv("data/merged_data/combined_metadata_filtered.csv")


# Merging and filtering as before
# Convert cc22_raw to Cork format
cc22_raw_list <- lapply(names(cc22_raw)[-1], function(colname) {
  data.frame(
    mass = cc22_raw$position,
    intensity = cc22_raw[[colname]]
  )
})

names(cc22_raw_list) <- names(cc22_raw)[-1]

# Combine raw CC22 and Cork spectra
combined_spectra_raw <- c(cc22_raw_list, peaksDFs)

# Standardize repeated isolate names
rename_map <- c(
  "A078 repeat" = "A078",
  "A095 #2" = "A095",
  "A115 Real" = "A115"
)

for (old_name in names(rename_map)) {
  new_name <- rename_map[[old_name]]
  if (old_name %in% names(combined_spectra_raw)) {
    combined_spectra_raw[[new_name]] <- combined_spectra_raw[[old_name]]
    combined_spectra_raw[[old_name]] <- NULL
  }
}



# Get matching names between spectra and metadata
spectra_names_raw <- names(combined_spectra_raw)
metadata_names <- merged_meta$sampleID

matched_names_raw <- intersect(spectra_names_raw, metadata_names)  # match to already filtered metadata

# Filter both spectra and metadata
combined_spectra_raw_filtered <- combined_spectra_raw[matched_names_raw]
combined_metadata_raw_filtered <- merged_meta[merged_meta$sampleID %in% matched_names_raw, ]

# Checks
all(mapply(nrow, combined_spectra_raw_filtered) == mapply(nrow, combined_spectra_filtered))  # exactly same no. of data points

length(combined_spectra_raw_filtered)  # perfect
nrow(combined_metadata_raw_filtered)  # all matches

str(combined_spectra_raw_filtered[["ASARM176"]])  # more proof
max(combined_spectra_raw_filtered[["ASARM176"]]$intensity, na.rm = TRUE)  # more proof

# Final Result

#save(combined_spectra_raw_filtered, file = "data/merged_data/combined_spectra_raw_filtered.RData")
load("data/merged_data/combined_spectra_raw_filtered.RData")  # called combined_spectra_filtered


