#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
####HEADER CODE ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#purpose
#initial commit date: 8/18/21
#original author: Tobias Koch
#email r.tobiaskoch@gmail.com or toby.koch@yale.edu

#PACKAGE INSTALL ####
packages = c("seqinr", "tidyverse", "lubridate", "xlsx", "rlist")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

rm(packages, package.check)

set.seed(1234)

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####FILE INPUT ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#fasta_file manually put into the folder
fasta_file = list.files("data_input/fasta_file/")[1]
fasta_input = read.fasta(paste("data_input/fasta_file/",
                               fasta_file,
                               sep = "")
                         )

#download metadata for fasta
seq_file = list.files("data_input/metadata/")[1]
seq_input = read.csv(paste("data_input/metadata/",
                            seq_file,
                            sep = "")
                     ) #make sure to download as csv

#lab names from glab metadata sheet
#only need to redownload if a new lab is added to the list
lab_file = list.files("data_input/",
                                  pattern = "*Lab names.csv")[1]
lab_names = read.csv(paste("data_input/",
                           lab_file,
                           sep = "")
                      )

rm(seq_file, fasta_file, lab_file)

#abbrevations of states
states = read.csv("state_abbreviation.csv") %>%
  select(-Abbrev) %>% #remove unnecessary column
  mutate(Country = "USA",
         Continent = "North America") %>%
  add_row(State = c("Puerto Rico", "US Virgin Islands", "Dominican Republic"), 
          Code = c("PR", "VI","DR"),
          Country = c("USA", "USA", "Dominican Republic"),
          Continent = rep("North America", 3)
          )

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####USER INPUT ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#column for argument this is your GISAID username
gisaid_id = "kb2228"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# DATA CLEANING AND FORMATTING ####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
####FORMATTING COLUMNS ####
submission = seq_input %>%
  left_join(states, by = c("Division..state." = "State")) %>%
  left_join(lab_names, by = "Source") %>%
  mutate(Collection.date = as.Date(Collection.date)) %>%
  transmute(Submitter       = gisaid_id,
          fn              = Sample.ID,
          covv_virus_name = str_c("hCov-19/", Country.y,
                                  "/",
                                  Code,
                                  "-",
                                  Sample.ID,
                                  "/",
                                  year(Collection.date)
                                  ),
          covv_type = "betacoronavirus",
          covv_passage = "Original",
          covv_collection_date = as.character(
                                      as.Date(Collection.date)
                                      ),
          covv_location = str_c(Continent,
                            "/",
                            Country.y,
                            "/",
                            Division..state.,
                            "/",
                            Location..county.
                            ),
          covv_add_location = "",
          covv_host = "Human",   #the next few columns are static that don't change for submissions
          covv_add_host_info = "",
          covv_sampling_strategy = "",
          covv_gender = "unknown",
          covv_patient_age = "unknown",
          covv_patient_status = "unknown",
          covv_specimen = "",
          covv_outbreak = "",
          covv_last_vaccinated = "",
          covv_treatment = "",
          covv_seq_technology = "Illumina NovaSeq",
          covv_assembly_method = "iVar",
          covv_coverage = "20x",
          covv_orig_lab = Originating.Lab,
          covv_orig_lab_addr = Address,
          covv_provider_sample_id = "unknown",
          covv_subm_lab = "Grubaugh Lab - Yale School of Public Health",
          covv_subm_lab_addr = "60 College St, New Haven, CT 06510, USA",
          covv_subm_sample_id = "unknown",
          covv_authors = Authors,
          covv_comment = ""
        )


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# FASTA FILE HEADER SWITCH####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
fa_df = names(fasta_input) %>% #extracts names from fasta file
  tibble() %>% #converts to a tibble
  rename(old = ".") %>% #renames for easier reading
  extract(old, "fasta", regex = "(Yale\\-[0-9]{4})") %>% #extracts Yale-ID from fasta file
  left_join(submission[,2:3], by = c("fasta" = "fn")) #joins by Yale-ID from submission file 
                                                      #to get the Virus Name for the GISAID Submission


#converts to a list to replace the names in the fasta file
fa_names = as.list(fa_df$covv_virus_name)

#creates copy for QC check
fasta_input2 = fasta_input

#rewrites the fasta file name
names(fasta_input2) = fa_names

#removes any NA's in the fasta sample name
#this would indicate a sample that was in run that shouldn't have been and isn't in our data
#for samples that are in the metadata file but not in the fasta file please see removed
fasta_input3 = list.subset(fasta_input2, 
                           !is.na(names(fasta_input2)
                                 )
                           )
#creates dataframe filtered names for filtering the metadata submission
fa_df_filter = data.frame(covv_virus_name = names(fasta_input3))
                           


#writes to output for new file
write.fasta(sequences = fasta_input3, 
            names = names(fasta_input3),
            file.out = "data_output/fasta_header_switch.fasta")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#FILTER METADATA####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
submission2 = submission %>%
  semi_join(fa_df_filter, by = "covv_virus_name") #return rows in submission that have an id match in fasta

write.csv(submission2, "data_output/GISAID_submission.csv")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#QC CHECK####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#REMOVED SAMPLES DATA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#counts number of removed samples.
#These are samples that are in the metadata but not in the fasta file
removed =
  anti_join(submission,fa_df, by = c("fn" = "fasta"))

write.csv(removed, "data_output/removed_check.csv")


#REMOVED SAMPLES LAB SUMMARY
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#creates summary table of the labs that have metadata being submitted
#initially and after samples were removed
lab_summary = data.frame(table(submission$covv_orig_lab) )

lab_summary_final = data.frame(table(submission2$covv_orig_lab) ) 

lab_summary_comb = full_join(lab_summary,lab_summary_final, by = "Var1") %>%
  mutate(removed = Freq.x - Freq.y)

write.csv(lab_summary_comb,"data_output/lab_summary.csv")

#LAB NAMES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#missing lab names from glab sheet on the metadata
missing_lab = submission %>%
  filter(is.na(covv_orig_lab))

#NEW FASTA ID MATCH
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#makes new id for QC
new_id = str_extract(c(names(fasta_input2)),
                       "Yale\\-[0-9]{4}")

#makes old id for QC
old_id = str_extract(c(names(fasta_input)),
                     "Yale\\-[0-9]{4}")

#creates logical vector to check if the yale-ids match
id_check = new_id == old_id


#gets count of number of samples that match
match = sum(id_check, na.rm = T)

#gets count of number of samples that match
no_match = length(fasta_input2) - sum(id_check, na.rm = T)

#gets count of number of samples are in fasta but missing from the metadata
na = sum(is.na(id_check))

#checks if number of samples in gisaid match the number of samples from the fasta file
lmatch = nrow(submission2) == length(fasta_input3)


#printed report of QC
qc_report = print(paste(
      "\n",
      nrow(seq_input), "samples were in the initial metadata file.\n",
      match, "samples yale-ids matched between the metadata and the fasta file.\n",
      nrow(removed), "samples were removed from the sequencing run due to low coverage. \n",
      nrow(missing_lab), "samples had no match to existing lab names.\n",
      na, "sample id's from fasta had no id in the metadata.\n",
      "\n Statement: The length of fasta file and metadata file being equal is", lmatch)
)

write.csv(qc_report, "data_output/qc_report.txt")  



  
  
