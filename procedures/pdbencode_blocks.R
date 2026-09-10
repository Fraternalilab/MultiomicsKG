#!/usr/bin/env Rscript

library("PDBencode")
library("bio3d")

# options("warn"=1)

encode_range = function(chain, range) {
  # message(paste("Chain", chain))
  ca.inds = bio3d::atom.select(str_bio3d, "calpha", chain=chain, resno=range, verbose=FALSE)
  ## in case the chain does not contain enough protein sequence
  if (length(ca.inds$atom) < 4) {
    message("  Chain not encoded: <4 CA atoms")
    return()
  }
  # subset PDB structure with indices
  str_bio3d_ca = trim.pdb(str_bio3d, ca.inds)
  
  #_______________________________________________________________________________
  ## encode structure
  sa_char.v = encode(str_bio3d_ca$xyz)
  ## from vector of chars to (vector of) string
  sa_string.v = paste(sa_char.v, collapse = '')
  # message(paste("Encoded sequence range:", sa_string.v))
  return(sa_string.v)
}

encode_chain = function(chain) {
  message(paste("Chain", chain))
  ca.inds = bio3d::atom.select(str_bio3d, "calpha", chain=chain, verbose=TRUE)
  ## in case the chain does not contain enough protein sequence
  if (length(ca.inds$atom) < 4) {
    message("  Chain not encoded: <4 CA atoms")
    return()
  }
  # subset PDB structure with indices
  str_bio3d_ca = trim.pdb(str_bio3d, ca.inds)
  
  #_______________________________________________________________________________
  ## encode structure
  sa_char.v = encode(str_bio3d_ca$xyz)
  ## from vector of chars to (vector of) string
  sa_string.v = paste(sa_char.v, collapse = '')
  message(paste("Encoded sequence:", sa_string.v))
  return(sa_string.v)
}

load_structure = function(file_path) {
  # filepath = paste(dir, file_name, sep = '/')
  # message(paste("Loading structure: ", file_path))
  str_bio3d <<- bio3d::read.pdb(file_path)
}
