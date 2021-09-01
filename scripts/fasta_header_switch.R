
#PACKAGE INSTALL ####
packages = c("seqinr", "tidyverse")


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

#fasta_file manually put into the folder
fasta_file = list.files(pattern = ".fa*")
fasta_input = read.fasta(fasta_file)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# FASTA FILE HEADER SWITCH####
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
names(fast) = names(fasta_input) %>%
  tibble() %>%
  rename(old = ".") %>%
  extract(old, "new", regex = "(Yale\\-[0-9]{4})") %>%
  
  

  
  