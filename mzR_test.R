# massProcesser demo data
file_dir <- "D:/fudan/Projects/2024/MetaboProcess/Progress/build_package/test_data/massProcesser/"
file_path <- list.files(file_dir, pattern = ".mzXML")
file_path <- paste0(file_dir, file_path)

filename <- file_path[1]
fileh <- mzR::openMSfile(filename, backend = NULL)
mzR::instrumentInfo(fileh)
mzR::runInfo(fileh)
mzR::header(fileh)
allSpect <- mzR::peaks(fileh, 1:2)
mzR::close()
