### helperfuncs.R                   
### A collection of auxiliary functions for various tasks
###
### Copyright: Jorge Gonzalez, 2012.
### Last modification: 25-05-2012.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### Author contact information:
###
###      Jorge Gonzalez B.
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22 
###      Santiago
###      Chile
###      Voice: +56-2-3545467  URL  : http://www.mat.puc.cl/~jgonzale
###      Fax  : +56-2-3547729  Email: jgonzale@mat.puc.cl
###


is.whole <- function(x){
  return((x %% 1) == 0)
}

#' Take a matrix and sum blocks of rows
#'
#' The original data set contains very long column headers. This function
#' does a keyword search over the headers to find those column headers that
#' match a particular keyword, e.g., mean, median, etc.
#' @param mat Input matrix
#' @param blocksize Size of the row blocks
#' @param w (Optional) Vector for weighted sum
#' @return Matrix
#' @export
rowBlockSum <- function(mat,blocksize,w=NULL){
  jumps <- dim(mat)[1]/blocksize
  
  if(!is.whole(jumps))
    stop("Number of rows have to be divisible by Block Size")
  
  if(!is.null(w)){ # If have non-null weights
    if(length(w) != jumps){ # And have incorrect size
      stop("Incorrect size of weights vector")
    }
  } else{
    w = rep(1, jumps)
  }
  
  ret <- Reduce('+', lapply(0:(jumps-1), 
                              function(i){ 
                                w[i+1]*mat[(1:(blocksize)+i*(blocksize)), ] 
                              } 
                            ))
  
  return(ret)
}

#' Transform a table from String to Data Frame.
#'
#' The usual way to compare results is using tables from books. This method
#' accept as input a pasted String, usually copied from a LaTeX table, and
#' create a Data Frame. Assumes the input is separated evenly by a character,
#' default is an empty space.
#' 
#' @param string Input string
#' @param n_cols Number of columns in the table
#' @param sep Separator character
#' @param header Have a header?
#' @return data.frame
#' @export
pasted_table_to_df <- function(string, n_cols, sep=" ", header=TRUE){
  raw_data <- strsplit(string, " ")[[1]]
  raw_data <- gsub("???", "-", raw_data)
  if(header){
    df <- data.frame(matrix(as.numeric(raw_data[-c(1:n_cols)]), 
                            ncol=n_cols, byrow=T))
    colnames(df) <- raw_data[1:n_cols]
  }
  else{
    df <- data.frame(matrix(as.numeric(raw_data), ncol=n_cols, byrow=T))
  }
  
  return(df)
}