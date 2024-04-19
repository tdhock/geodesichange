col.name.list <- list(
  loss=c(
    "penalty", "segments", "n.data", "n.lines",
    "mean.pen.cost", "total.loss", 
    "mean.intervals", "max.intervals"),
  segments=c("start", "end", "param"),
  data=c("start", "end", "radians")
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
  end <- cumsum(z.rle.vec$lengths)
  data.df <- data.frame(
    start=c(0L, end[-length(end)]),
    end,
    radians=z.rle.vec$values)
  geodesicFPOP_df(data.df, pen.num)
}

writeData <- function(radians.df, data.csv){
  if(!is.data.frame(radians.df)){
    stop("radians.df must be data.frame")
  }
  exp.names <- col.name.list$data
  if(!identical(names(radians.df), exp.names)){
    stop("radians.df must have names ", paste(exp.names, collapse=", "))
  }
  if(!is.integer(radians.df$start)){
    stop("radians.df$start must be integer")
  }
  if(!is.integer(radians.df$end)){
    stop("radians.df$end must be integer")
  }
  if(!is.numeric(radians.df$radians)){
    stop("radians.df$radians must be numeric")
  }
  if(any(radians.df$start < 0)){
    stop("radians.df$start must always be non-negative")
  }
  if(!all(radians.df$start < radians.df$end)){
    stop("start must be less than end for all rows of radians.df")
  }
  write.table(
    radians.df, data.csv,
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

geodesicFPOP_df <- function(radians.df, pen.num, base.dir=tempdir()){
  if(!(
    is.numeric(pen.num) &&
    length(pen.num)==1 &&
    0 <= pen.num)){
    stop("pen.num must be non-negative numeric scalar")
  }
  data.dir <- file.path(
    base.dir,
    with(radians.df, sprintf(
      "%d-%d", min(start), max(end))))
  unlink(data.dir, recursive=TRUE)
  dir.create(data.dir, showWarnings=FALSE, recursive=TRUE)
  data.csv <- file.path(data.dir, "data.csv")
  writeData(radians.df, data.csv)
  L <- geodesicFPOP_dir(data.dir, paste(pen.num))
  L$data <- data.table(radians.df)
  L
}

geodesicFPOP_dir <- function(problem.dir, penalty.param, db.file=NULL){
  megabytes <- start <- NULL
  if(!(
    is.character(problem.dir) &&
    length(problem.dir)==1 &&
    dir.exists(problem.dir))){
    stop(
      "problem.dir=", problem.dir,
      " must be the name of a directory",
      " containing a file named data.csv")
  }
  if(!(
    (is.numeric(penalty.param) || is.character(penalty.param)) &&
    length(penalty.param)==1 &&
    (!is.na(penalty.param))
  )){
    stop("penalty.param must be numeric or character, length 1, not missing")
  }
  penalty.str <- paste(penalty.param)
  prob.cov.csv <- file.path(problem.dir, "data.csv")
  pre <- paste0(prob.cov.csv, "_penalty=", penalty.str)
  penalty_segments.csv <- paste0(pre, "_segments.csv")
  penalty_loss.tsv <- paste0(pre, "_loss.tsv")
  penalty_timing.tsv <- paste0(pre, "_timing.tsv")
  already.computed <- tryCatch({
    timing <- fread(
      file=penalty_timing.tsv,
      col.names=c("penalty", "megabytes", "seconds"))
    first.seg.line <- fread.first(penalty_segments.csv, col.name.list$segments)
    last.seg.line <- fread.last(penalty_segments.csv, col.name.list$segments)
    first.cov.line <- fread.first(prob.cov.csv, col.name.list$data)
    last.cov.line <- fread.last(prob.cov.csv, col.name.list$data)
    penalty.loss <- fread(file=penalty_loss.tsv, col.names=col.name.list$loss)
    nrow.ok <- nrow(timing)==1 && nrow(penalty.loss)==1 &&
      nrow(first.seg.line)==1 && nrow(last.seg.line)==1 &&
      nrow(first.cov.line)==1 && nrow(last.cov.line)==1
    loss.segments.consistent <-
      first.seg.line$end-last.seg.line$start == penalty.loss$bases
    ## segments files are written by decoding/backtracking after
    ## dynamic progamming, so it is normal that the first line of the
    ## segment file is actually the last segment in terms of position
    ## in the data sequence.
    start.ok <- first.cov.line$start == last.seg.line$start
    end.ok <- last.cov.line$end == first.seg.line$end
    nrow.ok && loss.segments.consistent && start.ok && end.ok
  }, error=function(e){
    FALSE
  })
  if(!already.computed){
    seconds <- system.time({
      result <- geodesicFPOP_file(prob.cov.csv, penalty.str, db.file)
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
    file=penalty_segments.csv,
    col.names=col.name.list$segments),
    start)
  L <- list(
    segments=penalty.segs,
    loss=data.table(
      penalty.loss,
      timing[, list(megabytes, seconds)]))
  L
}

geodesicFPOP_file <- function(csv.file, pen.str, db.file=NULL){
  if(!(
    is.character(csv.file) &&
    length(csv.file)==1 &&
    file.exists(csv.file)
  )){
    stop(
      "csv.file=", csv.file,
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
  norm.file <- normalizePath(csv.file, mustWork=TRUE)
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
  prefix <- paste0(csv.file, "_penalty=", pen.str)
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
