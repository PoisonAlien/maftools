TCGAmutations::tcga_load(study = "brca")
TCGAmutations::tcga_load(study = "brca", source = "Firehose")

x = subsetMaf(maf = tcga_brca, tsb = c("TCGA-AN-A046", "TCGA-AC-A23H"), mafObj = F)[,.(Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode)]
x[,Chromosome := as.character(Chromosome)]
x[,Start_Position := as.numeric(as.character(Start_Position))]
x[,End_Position := as.numeric(as.character(End_Position))]

y = subsetMaf(maf = tcga_brca_mc3, tsb = c("TCGA-AN-A046", "TCGA-AC-A23H"), mafObj = F)[,.(Chromosome, Start_Position, End_Position,Tumor_Sample_Barcode)]
y[,Chromosome := as.character(Chromosome)]
y[,Start_Position := as.numeric(as.character(Start_Position))]
y[,End_Position := as.numeric(as.character(End_Position))]

x[,id := paste(Chromosome, Start_Position, End_Position, sep = ":")]
y[,id := paste(Chromosome, Start_Position, End_Position, sep = ":")]

data.table::setkey(x = x, Chromosome, Start_Position, End_Position)
data.table::setkey(x = y, Chromosome, Start_Position, End_Position)


xINy = x[x$id %in% y$id] #Common elements among x and y
nrow(xINy)
xUNIQ = x[!x$id %in% y$id] #Elements unique to x
nrow(xUNIQ)

xEQyInd = data.table::foverlaps(
  x = x,
  y = y,
  type = "equal",
  nomatch = NULL,
  mult = "all",
  which = TRUE
)

xEQy = x[unique(xEQyInd$xid),]

xextra = xEQy[!id %in% xINy$id]
