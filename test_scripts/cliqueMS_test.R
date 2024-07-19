mzfile <- system.file("standards.mzXML", package = "cliqueMS")
data(positive.adinfo, package = "cliqueMS")
positive.adinfo <- positive.adinfo[positive.adinfo$adduct %in%
                                     c("[M+H]+", "[M+H-H2O]+", "[M+Na]+", "[M+H-NH3]+", "[M+K]+", "[M+NH4]+"),]
data(negative.adinfo, package = "cliqueMS")
mzraw <- MSnbase::readMSData(files = mzfile, mode = "onDisk")
cpw <- xcms::CentWaveParam(ppm = 15, peakwidth = c(5, 20), snthresh = 10)
mzData <- xcms::findChromPeaks(mzraw, cpw)
en.anClique <- cliqueMS::createanClique(mzData)
set.seed(2)
ex.cliqueGroups <- cliqueMS::getCliques(mzData, filter = TRUE)
ex.Isotopes <- cliqueMS::getIsotopes(ex.cliqueGroups, ppm = 10)
ex.Adducts <- cliqueMS::getAnnotation(ex.Isotopes, ppm = 10,
                                      adinfo = positive.adinfo, polarity = "positive", normalizeScore = TRUE)
test2 <- ex.Adducts@peaklist[, c("mz", "rt", "maxo", "sample", "cliqueGroup", "isotope", "mass1", "an1")]
ex.Adducts@cliques
