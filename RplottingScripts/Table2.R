#TABLE 2

# Load necessary libraries
library(knitr)
library(kableExtra)

# Create the data frame
hla_data <- data.frame(
  Allele = c("HLA-A29:02", "HLA-B27:05", "HLA-B35:01", "HLA-B50:01", "HLA-B51:01", 
             "HLA-B13:02", "HLA-B07:02", "HLA-C12:03", "HLA-C05:01", "DRB1_0801", 
             "DRB1_1302", "DRB1_0202", "DRB1_0803", "DRB1_0302", "DRB1_0304", 
             "DRB1_1401", "DQB10201", "DQB10202"),
  Primary_Risk_Associations = c("Birdshot uvitis", "Ankylosing spondylitis", "Autoimmune hepatitis, subacute thyroiditis", 
                                "Myasthenia gravis", "Behçet's syndrome", "Psoriasis", 
                                "Ankylosing spondylitis", "psoriasis", "multiple sclerosis", "Primary biliary cirrhosis", 
                                "Autoimmune hepatitis type 1, early childhood myasthenia gravis", "Osteoarthritis", 
                                "Primary biliary cirrhosis, sarcoidosis", "Autoimmune hepatitis type 1", "Graves' disease", 
                                "Psoriasis vulgaris", "Celiac Disease", "Celiac Disease"),
  PMID = c("PMID: 28314830, PMID: 25434765, PMID: 33262772, PMID: 4127279", 
           "PMID: 4123836, PMID: 25861975, PMID: 28188227, PMID: 24838411", 
           "PMID: 38089552, PMID: 17383147, PMID: 15307945, PMID: 1822362", 
           "PMID: 27802446, PMID: 17480220", 
           "PMID: 34690278, PMID: 9568801", 
           "PMID: 8932582, PMID: 21475522, PMID: 22121262", 
           "PMID: 27254288, PMID: 19655145", 
           "PMID: 25087609, PMID: 31341424, PMID: 16235096, PMID: 35693758", 
           "PMID: 19412418", 
           "PMID: 16941709", 
           "PMID: 15003812, PMID: 11182227", 
           "PMID: 12048293", 
           "PMID: 7699227, PMID: 20685690", 
           "PMID: 17105585", 
           "PMID: 10468973, PMID: 9768636, PMID: 16254435", 
           "PMID: 9990359", 
           "PMID: 17919990, PMID: 17190762", 
           "PMID: 17919990, PMID: 17190762")
)

# Create the table
kable(hla_data, "html", caption = "Table 2: Autoimmune HLA Alleles") %>%
  kable_styling(full_width = FALSE, position = "left") %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(3, width = "30em")
