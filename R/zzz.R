.onLoad <- function(libname = find.package("DPpackage"), pkgname = "DPpackage"){
  if(Sys.getenv("R_ARCH") == "\x64"){
    apnd <- paste0(Sys.getenv("R_LIBS_USER"), '/DPpackage/libs/x64')
  } else {
    apnd <- paste0(Sys.getenv("R_LIBS_USER"), '/DPpackage/libs/i386')
  }
  
  PATH.1 <- paste(Sys.getenv("PATH"), gsub("/", "\\\\", apnd), sep=";")
  Sys.setenv("PATH"=PATH.1)
}

.onUnload <- function(libname, pkgname){
  
}
