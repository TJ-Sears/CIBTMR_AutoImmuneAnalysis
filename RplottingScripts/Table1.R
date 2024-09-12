# Load necessary libraries
library(knitr)
library(kableExtra)
library(dplyr)

# Create the data frame
characteristics_data <- data.frame(
  Characteristics = c("No. of patients", "No. of centers", "Patient age (year) - median (min-max)", 
                      "Sex - no. (%)", "Male", "Race/ethnicity - no. (%)", "Caucasian, non-Hispanic", 
                      "KPS - no. (%)", "0-90", "HCT-CI- no. (%)", "0-2", "3+", "Missing", 
                      "Pre-transplant therapies/ (%)", "HMA alone", "Chemo alone", "HMA plus Chemo", 
                      "Neither", "Missing", "MDS IPSS-R score pre transplant - no. (%)", "Very low", 
                      "Low", "Intermediate", "High", "Very high", "Missing", "Time from diagnosis to HCT (month) - median (range)", 
                      "Donor type - no. (%)", "HLA-identical sibling", "Other related", 
                      "Well-matched unrelated (8/8)", "Partially-matched unrelated (7/8)", 
                      "Mis-matched unrelated (≤6/8) or unknown", "stem cell source - no. (%)", 
                      "Bone marrow", "Peripheral blood", "regimen intensity - no. (%)", "Myeloablative", 
                      "Reduced intensity", "Non-myeloablative", "Missing", "Year of HCT - no. (%)", 
                      "2014", "2015", "2016", "2017", "2018", "Median follow up of survivors (months) - median (range)"),
  Patients = c("494", "93", "66 (22-78)", "", "315 (64)", "", "494 (100)", "", "251(51)", "", "171(35)", 
               "315(64)", "8(2)", "", "340 (69)", "15 (3)", "34 (7)", "93 (19)", "12 (2)", "", "57 (12)", 
               "123 (25)", "160 (32)", "74 (15)", "22 (4)", "58 (12)", "18(2-263)", "", "65 (13)", "32 (6)", 
               "353 (71)", "39 (8)", "5 (1)", "", "59(12)", "435(88)", "", "127(26)", "308(62)", 
               "44(9)", "15(3)", "", "91(18)", "156(32)", "118(24)", "122(25)", "7(1)", "34.5(3.2-62.7)")
)

# Create a column to identify subheaders
characteristics_data <- characteristics_data %>%
  mutate(is_subheader = ifelse(Patients == "", TRUE, FALSE))

# Extract row indices for styling
subheader_rows <- which(characteristics_data$is_subheader)
regular_rows <- which(!characteristics_data$is_subheader)

# Remove the is_subheader column for the final table output
characteristics_data <- characteristics_data %>%
  select(-is_subheader)

# Create the table
kable(characteristics_data, "html", escape = FALSE, caption = "Patient Characteristics") %>%
  kable_styling(full_width = FALSE, position = "left", font_size = 10) %>%
  column_spec(1, width = "30em") %>%
  row_spec(subheader_rows, bold = TRUE, background = "#f2f2f2", extra_css = "padding-top: 2px; padding-bottom: 2px;") %>%
  row_spec(regular_rows, extra_css = "padding-left: 20px; padding-top: 2px; padding-bottom: 2px;") %>%
  row_spec(0, extra_css = "padding-top: 2px; padding-bottom: 2px;")
