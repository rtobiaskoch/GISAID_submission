####HEADER CODE ####
#purpose
#initial commit date: 8/18/21
#original author: Tobias Koch
#email r.tobiaskoch@gmail.com or toby.koch@yale.edu

#PACKAGE INSTALL ####
packages = c("readr", "tidyverse", "lubridate","EpiEstim", "xlsx","stats","zoo",
             "DescTools", "googledrive","googlesheets4")



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

rm(packages)

set.seed(1234)

commandArgs()


####FILE INPUT ####

seq_filename = list.files(pattern = "*.metadata.xlsx")
seq_input = read.xlsx(seq_filename, )
mdata_filename = list.files(pattern = "data_GLab_SC2")
mdata_input = read.xlsx(mdata_filename, sheet = "Sample metadata")

#column for argument this is your GISAID username
Submitter = Arg


####STATIC COLS ####
Type
Passage details/history
Host
Gender
Patient age
Patient status
Sequencing technology
Assembly method
Coverage
Sample ID given by originating laboratory
Submitting lab
Submitting lab address

####DYNAMIC COLS ####
