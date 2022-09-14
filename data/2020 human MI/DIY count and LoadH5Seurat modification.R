library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(scater)
library(reticulate)
library(rhdf5)
library(Matrix)
setwd("D:/LocalGitHub/scRNA-seq/2020 human MI")


### Generate spDIY from Matrix
obj_HDF5<- Connect("raw/snRNA-seq-submission.h5seurat", mode = "r+")
x <- obj_HDF5[["assays/RNA/data"]]
spDIY <- Matrix::sparseMatrix(
  i=x[["indices"]][]+1,
  p=x[["indptr"]][],
  x=x[["data"]][]  )

dim(spDIY)

t_seurat_spDIY = SeuratDisk::Transpose(spDIY)
#took two days...used in AssembleAssay function modification in another file
#and that returns an object "obj", which is used in the loading function below


## Call LoadH5Seurat
file= "raw/snRNA-seq-submission.h5seurat"
hfile <- h5Seurat$new(filename = file, mode = 'r')

myLoader <- function(
    x = hfile
    assays = NULL
    reductions = NULL
    graphs = NULL
    neighbors = NULL
    images = NULL
    meta.data = TRUE
    commands = TRUE
    # misc = TRUE
    misc = F
    tools = TRUE
    verbose = TRUE
    ...
) {
  index <- x$index()
  obj.all <- all(vapply(
    X = c(assays, reductions, graphs),
    FUN = is.null,
    FUN.VALUE = logical(length = 1L)
  ))
  # Load Assays
  assays <- GetAssays(assays = assays, index = index)
  if (!DefaultAssay(object = index) %in% names(x = assays)) {
    active.assay <- names(x = assays)[1]
    warning(
      "Default assay not requested, using ",
      active.assay,
      " instead",
      call. = FALSE,
      immediate. = TRUE
    )
  } else {
    active.assay <- DefaultAssay(object = index)
  }
  assay.objects <- vector(mode = 'list', length = length(x = assays))
  names(x = assay.objects) <- names(x = assays)
  #🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞🐞
  for (assay in names(x = assays)) {
    # assay.objects[[assay]] <- AssembleAssay(
    #   assay = assay,
    #   file = x,
    #   slots = assays[[assay]],
    #   verbose = verbose
    # )
    assay.objects[[assay]] <- obj
  }
  default.assay <- list(assay.objects[[active.assay]])
  names(x = default.assay) <- active.assay
  object <- new(
    Class = 'Seurat',
    assays = default.assay,
    active.assay = active.assay,
    meta.data = data.frame(row.names = Cells(x = x)),
    version = package_version(x = x$version())
  )
  for (assay in names(x = assay.objects)) {
    if (assay != active.assay) {
      object[[assay]] <- assay.objects[[assay]]
    }
  }
  # Load DimReducs
  reductions <- GetDimReducs(
    reductions = reductions,
    index = index,
    assays = assays
  )
  for (reduc in reductions) {
    if (verbose) {
      message("Adding reduction ", reduc)
    }
    reduction <- AssembleDimReduc(
      reduction = reduc,
      file = x,
      verbose = verbose
    )
    if (isTRUE(x = getOption(x = 'SeuratDisk.dimreducs.allglobal', default = FALSE))) {
      slot(object = reduction, name = 'global') <- TRUE
    }
    object[[reduc]] <- reduction
  }
  # Load Graphs
  graphs <- GetGraphs(graphs = graphs, index = index, assays = assays)
  for (graph in graphs) {
    if (verbose) {
      message("Adding graph ", graph)
    }
    object[[graph]] <- AssembleGraph(graph = graph, file = x, verbose = verbose)
  }
  # Load Neighbors
  neighbors <- GetNeighbors(neighbors = neighbors, index = index)
  for (neighbor in neighbors) {
    if (verbose) {
      message("Adding neighbors ", neighbor)
    }
    object[[neighbor]] <- AssembleNeighbor(
      neighbor = neighbor,
      file = x,
      verbose = verbose
    )
  }
  # Load SpatialImages
  if (packageVersion(pkg = 'Seurat') >= numeric_version(x = spatial.version)) {
    images <- GetImages(images = images, index = index, assays = assays)
    for (image in images) {
      if (verbose) {
        message("Adding image ", image)
      }
      object[[image]] <- AssembleImage(
        image = image,
        file = x,
        verbose = verbose
      )
    }
  }
  # Load SeuratCommands
  if (commands) {
    if (verbose) {
      message("Adding command information")
    }
    cmds <- GetCommands(index = index, assays = assays)
    cmdlogs <- vector(mode = 'list', length = length(x = cmds))
    names(x = cmdlogs) <- cmds
    for (cmd in cmds) {
      cmdlogs[[cmd]] <- AssembleSeuratCommand(
        cmd = cmd,
        file = x,
        verbose = verbose
      )
    }
    slot(object = object, name = 'commands') <- cmdlogs
  }
  # Load meta.data
  if (meta.data) {
    if (verbose) {
      message("Adding cell-level metadata")
    }
    md <- as.data.frame(x = x[['meta.data']], row.names = Cells(x = x))
    if (ncol(x = md)) {
      object <- AddMetaData(object = object, metadata = md)
    }
  }
  # Set cell identities and object project
  Idents(object = object) <- Idents(object = x)
  Project(object = object) <- Project(object = x)
  # Load misc
  if (misc) {
    if (verbose) {
      message("Adding miscellaneous information")
    }
    slot(object = object, name = 'misc') <- as.list(x = x[['misc']])
  }
  # Load tools
  if (tools) {
    if (verbose) {
      message("Adding tool-specific results")
    }
    slot(object = object, name = 'tools') <- as.list(x = x[['tools']])
  }
  # Load no.assay information
  if (obj.all && !is.null(x = index$no.assay)) {
    if (verbose) {
      message("Adding data that was not associated with an assay")
    }
    for (graph in index$no.assay$graphs) {
      object[[graph]] <- AssembleGraph(
        graph = graph,
        file = x,
        verbose = verbose
      )
    }
    for (cmd in index$no.assay$commands) {
      object[[cmd]] <- AssembleSeuratCommand(
        cmd = cmd,
        file = x,
        verbose = verbose
      )
    }
  }
  return(object)
}

