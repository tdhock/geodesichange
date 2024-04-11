col.name.list <- list(
  loss=c(
    "penalty", "segments", "bases", "bedGraph.lines",
    "mean.pen.cost", "total.loss", 
    "mean.intervals", "max.intervals"),
  segments=c("chrom","chromStart", "chromEnd", "status", "param"),
  coverage=c("chrom", "chromStart", "chromEnd", "count")
)

geodesicFPOP_vec <- function(angle.vec, pen.num){
  if(!(
    is.numeric(pen.num) &&
    length(pen.num)==1 &&
    0 <= pen.num)){
    stop("pen.num must be non-negative numeric scalar")
  }
  if(!is.numeric(angle.vec)){
    stop("angle.vec must be integer")
  }
  z.rle.vec <- rle(angle.vec)
  chromEnd <- cumsum(z.rle.vec$lengths)
  coverage.df <- data.frame(
    chrom="chrUnknown",
    chromStart=c(0L, chromEnd[-length(chromEnd)]),
    chromEnd,
    count=z.rle.vec$values)
  geodesicFPOP_df(coverage.df, pen.num)
}

writeBedGraph <- function(count.df, coverage.bedGraph){
  if(!is.data.frame(count.df)){
    stop("count.df must be data.frame")
  }
  exp.names <- c("chrom", "chromStart", "chromEnd", "count")
  if(!identical(names(count.df), exp.names)){
    stop("count.df must have names ", paste(exp.names, collapse=", "))
  }
  if(!is.integer(count.df$chromStart)){
    stop("count.df$chromStart must be integer")
  }
  if(!is.integer(count.df$chromEnd)){
    stop("count.df$chromEnd must be integer")
  }
  if(!is.numeric(count.df$count)){
    stop("count.df$count must be numeric")
  }
  if(any(count.df$chromStart < 0)){
    stop("count.df$chromStart must always be non-negative")
  }
  if(!all(count.df$chromStart < count.df$chromEnd)){
    stop("chromStart must be less than chromEnd for all rows of count.df")
  }
  write.table(
    count.df, coverage.bedGraph,
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

geodesicFPOP_df <- function(count.df, pen.num, base.dir=tempdir()){
  if(!(
    is.numeric(pen.num) &&
    length(pen.num)==1 &&
    0 <= pen.num)){
    stop("pen.num must be non-negative numeric scalar")
  }
  data.dir <- file.path(
    base.dir,
    with(count.df, sprintf(
      "%s-%d-%d", chrom[1], min(chromStart), max(chromEnd))))
  unlink(data.dir, recursive=TRUE)
  dir.create(data.dir, showWarnings=FALSE, recursive=TRUE)
  coverage.bedGraph <- file.path(data.dir, "coverage.bedGraph")
  writeBedGraph(count.df, coverage.bedGraph)
  L <- geodesicFPOP_dir(data.dir, paste(pen.num))
  L$data <- data.table(count.df)
  L
}

geodesicFPOP_dir <- function(problem.dir, penalty.param, db.file=NULL){
  megabytes <- NULL
  if(!(
    is.character(problem.dir) &&
    length(problem.dir)==1 &&
    dir.exists(problem.dir))){
    stop(
      "problem.dir=", problem.dir,
      " must be the name of a directory",
      " containing a file named coverage.bedGraph")
  }
  if(!(
    (is.numeric(penalty.param) || is.character(penalty.param)) &&
    length(penalty.param)==1 &&
    (!is.na(penalty.param))
  )){
    stop("penalty.param must be numeric or character, length 1, not missing")
  }
  penalty.str <- paste(penalty.param)
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  pre <- paste0(prob.cov.bedGraph, "_penalty=", penalty.str)
  penalty_segments.bed <- paste0(pre, "_segments.bed")
  penalty_loss.tsv <- paste0(pre, "_loss.tsv")
  penalty_timing.tsv <- paste0(pre, "_timing.tsv")
  already.computed <- tryCatch({
    timing <- fread(
      file=penalty_timing.tsv,
      col.names=c("penalty", "megabytes", "seconds"))
    first.seg.line <- fread.first(penalty_segments.bed, col.name.list$segments)
    last.seg.line <- fread.last(penalty_segments.bed, col.name.list$segments)
    first.cov.line <- fread.first(prob.cov.bedGraph, col.name.list$coverage)
    last.cov.line <- fread.last(prob.cov.bedGraph, col.name.list$coverage)
    penalty.loss <- fread(file=penalty_loss.tsv, col.names=col.name.list$loss)
    nrow.ok <- nrow(timing)==1 && nrow(penalty.loss)==1 &&
      nrow(first.seg.line)==1 && nrow(last.seg.line)==1 &&
      nrow(first.cov.line)==1 && nrow(last.cov.line)==1
    loss.segments.consistent <-
      first.seg.line$chromEnd-last.seg.line$chromStart == penalty.loss$bases
    ## segments files are written by decoding/backtracking after
    ## dynamic progamming, so it is normal that the first line of the
    ## segment file is actually the last segment in terms of position
    ## on the chromosome.
    start.ok <- first.cov.line$chromStart == last.seg.line$chromStart
    end.ok <- last.cov.line$chromEnd == first.seg.line$chromEnd
    nrow.ok && loss.segments.consistent && start.ok && end.ok
  }, error=function(e){
    FALSE
  })
  if(!already.computed){
    seconds <- system.time({
      result <- geodesicFPOP_file(prob.cov.bedGraph, penalty.str, db.file)
    })[["elapsed"]]
    timing <- data.table(
      penalty=as.numeric(penalty.str),
      megabytes=result$megabytes,
      seconds)
    write.table(
      timing,
      penalty_timing.tsv,
      row.names=FALSE, col.names=FALSE,
      quote=FALSE, sep="\t")
    penalty.loss <- fread(file=penalty_loss.tsv, col.names=col.name.list$loss)
  }
  penalty.segs <- setkey(fread(
    file=penalty_segments.bed,
    col.names=col.name.list$segments),
    chromStart)
  L <- list(
    segments=penalty.segs,
    loss=data.table(
      penalty.loss,
      timing[, list(megabytes, seconds)]))
  L
}

geodesicFPOP_file <- function(bedGraph.file, pen.str, db.file=NULL){
  if(!(
    is.character(bedGraph.file) &&
    length(bedGraph.file)==1 &&
    file.exists(bedGraph.file)
  )){
    stop(
      "bedGraph.file=", bedGraph.file,
      " must be the name of a data file to segment")
  }
  if(!is.character(pen.str)){
    stop(paste(
      "pen.str must be a character string",
      "that can be converted to a non-negative numeric scalar"
    ))
  }
  penalty <- as.numeric(pen.str)
  if(!(
    is.numeric(penalty) &&
    length(penalty)==1 &&
    0 <= penalty && penalty <= Inf
  )){
    stop("as.numeric(pen.str)=", penalty, " but it must be a non-negative numeric scalar")
  }
  norm.file <- normalizePath(bedGraph.file, mustWork=TRUE)
  if(is.null(db.file)){
    db.file <- sprintf("%s_penalty=%s.db", norm.file, pen.str)
  }
  if(!(
    is.character(db.file) &&
    length(db.file)==1
  )){
    stop(
      "db.file=", db.file,
      " must be a temporary file name where cost function db can be written")
  }
  unlink(db.file)
  print(norm.file)
  geodesicFPOP_interface(
    norm.file,
    pen.str,
    db.file)
  result <- list()
  result$megabytes <- if(file.exists(db.file)){
    file.size(db.file)/1024/1024
  }else{
    0
  }
  unlink(db.file)
  prefix <- paste0(bedGraph.file, "_penalty=", pen.str)
  loss.tsv <- paste0(prefix, "_loss.tsv")
  if(file.size(loss.tsv)==0){
    stop(
      "unable to write to loss output file ",
      loss.tsv,
      " (disk is probably full)"
    )
  }
  result
}
