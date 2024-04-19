#include "geodesicFPOP.h"
#include <Rcpp.h>

// [[Rcpp::export]]
void geodesicFPOP_interface
(const std::string in_file_str,
 const std::string penalty_str,
 const std::string temp_file_str){
  int status = geodesicFPOP
    (in_file_str.c_str(),
     penalty_str.c_str(),
     temp_file_str.c_str());
  if(status==ERROR_PENALTY_NOT_FINITE){
    Rcpp::stop("penalty=%s but must be finite", penalty_str);
  }
  if(status==ERROR_PENALTY_NEGATIVE){
    Rcpp::stop("penalty=%s must be non-negative", penalty_str);
  }
  if(status==ERROR_UNABLE_TO_OPEN_BEDGRAPH){
    Rcpp::stop("unable to open input file for reading %s", in_file_str);
  }
  if(status==ERROR_NOT_ENOUGH_COLUMNS){
    Rcpp::stop("each line of input data file %s should have exactly three space-separated columns: start, end, radians", in_file_str);
  }
  if(status==ERROR_INCONSISTENT_CHROMSTART_CHROMEND){
    Rcpp::stop("there should be no gaps (columns 1-2) in input data file %s", in_file_str);
  }
  if(status==ERROR_WRITING_COST_FUNCTIONS){
    Rcpp::stop("unable to write to cost function database file %s", temp_file_str);
  }
  if(status==ERROR_WRITING_LOSS_OUTPUT){
    Rcpp::stop("unable to write to loss output file %s_penalty=%s_loss.tsv",
	  in_file_str, penalty_str);
  }
  if(status==ERROR_WRITING_SEGMENTS_OUTPUT){
    Rcpp::stop("unable to write to segments output file %s_penalty=%s_segments.bed",
	  in_file_str, penalty_str);
  }
  if(status==ERROR_NO_DATA){
    Rcpp::stop("input file %s contains no data", in_file_str);
  }
  if(status==ERROR_PENALTY_NOT_NUMERIC){
    Rcpp::stop
      ("penalty string '%s' is not numeric; it should be convertible to double",
       penalty_str);
  }
  if(status != 0){
    Rcpp::stop("error code %d", status);
  }
}
