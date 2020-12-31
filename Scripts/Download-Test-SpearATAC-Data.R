#Download SpearATAC Test Data
aws <- "https://jeffgranja.s3.amazonaws.com/SpearATAC_2021"

dataFiles <- c(
	"K562_R1.fragments.tsv.gz",
	"K562_R1.fragments.tsv.gz.tbi",
	"K562_R1.sgRNA.rds",
	"K562_R1.singlecell.csv",
	"K562_R2.fragments.tsv.gz",
	"K562_R2.fragments.tsv.gz.tbi",
	"K562_R2.sgRNA.rds",
	"K562_R2.singlecell.csv"
)

for(i in seq_along(dataFiles)){
	download.file(
		url = file.path(aws, dataFiles[i]),
		destfile = file.path("data", dataFiles[i])
	)
}