reg_install_packages <- c('BiocManager', 'tidyverse', 'devtools', 'data.table')
install.packages(setdiff(reg_install_packages, rownames(installed.packages())) , repos = "http://cran.us.r-project.org")

if (!require("EnsDb.Mmusculus.v79", character.only = TRUE)){
    BiocManager::install("EnsDb.Mmusculus.v79")
}

if (!require("riboWaltz",character.only = TRUE)){
    devtools::install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)
}

all_packages <- c(reg_install_packages, 'riboWaltz', "EnsDb.Mmusculus.v79")

if (length(setdiff(all_packages, rownames(installed.packages()))) == 0){
    fileConn<-file("scripts/install_verification.txt")
    writeLines('All installed', fileConn)
    close(fileConn)
} else {
    fileConn<-file("scripts/failed_installations.txt")
    for (pkg in setdiff(all_packages, rownames(installed.packages()))){
        writeLines(pkg, fileConn)
    }
    close(fileConn) 
}


