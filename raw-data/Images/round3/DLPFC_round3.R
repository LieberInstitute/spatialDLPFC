path1 = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/spatialDLPFC/Images/round3"  #dont include forward slash at end
listOFfiles = Sys.glob(paste0(path1,"/100014*round3*.tif"))
write.table(listOFfiles,file = paste0(path1,"/DLPFC_round3.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)# stores the text file in the main data directory
