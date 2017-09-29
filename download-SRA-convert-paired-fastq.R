#Download experiment files (.sra)
#need SraRunInfo.csv before execution
#from https://www.ncbi.nlm.nih.gov/sra: Send to (top right corner), Select “File”,  Select format “RunInfo”,  Click on “Create File"
sri<-read.csv("SraRunInfo.csv", stringsAsFactors=FALSE) #read the .csv
files<-basename(sri$download_path) #in variable files, store the names of all the .sra files (by removing the whole path anme)
for(i in 1:length(files)) download.file(sri$download_path[i], files[i]) #download all the files and save them with their respective sra names

for(f in files) { #for all these files
	args_str <- paste("--split-3",f,sep=" ") #prepare arguments for console command
	system2("fastq-dump", args=args_str) #run console command to convert from .sra filetype to paired .fastq
}

