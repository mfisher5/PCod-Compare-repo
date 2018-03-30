####################### Convert stacks-called genotypes to 0101/0202 #######################
#
# For use with OutFLANK, so that the script convert_genepop_to_SNPmat.py can be run. 
#  The input file is the matrix of genotypes from genepop, *not including locus names*
#  The output should be copied and pasted back into the genepop file
#
# MF 3/24/2018
#
#############################################################################################


# Read in Data ------------------------------------------------------------
mydata <- read.table("batch_8_eastt_genotypes_forR.txt", header = FALSE, sep = " ", colClasses = "character")


# For Loop to Replace Genotypes -------------------------------------------
nloci <- ncol(mydata) - 1
nloci
for(i in seq(2,nloci+1)){
  genotypes <- unique(mydata[,i])
  if('0101' %in% genotypes | '0102' %in% genotypes | '0103' %in% genotypes | '0104' %in% genotypes){
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0303", "0202")
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0103", "0102")
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0404", "0202")
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0104", "0102")
  }
  else if('0202' %in% genotypes && '0303' %in% genotypes | '0203' %in% genotypes){
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0202", "0101")
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0303", "0202")
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0203", "0102")
  }
  else if ('0202' %in% genotypes && '0404' %in% genotypes | '0204' %in% genotypes){
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0202", "0101")
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0404", "0202")
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0204", "0102")
  } 
  else if ('0303' %in% genotypes && '0404' %in% genotypes | '0304' %in% genotypes){
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0303", "0101")
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0404", "0202")
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0304", "0102")
  }
  else{
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0303", "0101")
    mydata[,i] <- replace(mydata[,i], mydata[,i]=="0404", "0202")
  }
}


# Write out Data ----------------------------------------------------------
write.table(mydata, "east_genotypes_matrix_alleles12_edit.txt", quote=FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE)
