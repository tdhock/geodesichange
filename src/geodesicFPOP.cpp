#include <iomanip> //for setprecision.
#include <fstream> //for ifstream etc.
#include <exception>//for std::exception
#include <stdexcept>//for std::invalid_argument
#include <R.h> // Rprintf

#include "geodesicFPOP.h"

#include <list>
#include <math.h>
#include <stdio.h>
#include <R.h>

#define NEWTON_EPSILON 1e-12
#define NEWTON_STEPS 100
#define PREV_NOT_SET (-3)

#define ABS(x) ((x)<0 ? -(x) : (x))

LinearPiece::LinearPiece
(double li, double co, double m, double M, int i, double prev){
  Linear = li;
  Constant = co;
  min_angle_param = m;
  max_angle_param = M;
  data_i = i;
  prev_angle_param = prev;
}

LinearPiece::LinearPiece
(double li, double co, double m, double M){
  Linear = li;
  Constant = co;
  min_angle_param = m;
  max_angle_param = M;
  data_i = PREV_NOT_SET;
  prev_angle_param = INFINITY;
}

LinearPiece::LinearPiece(){
}

double LinearPiece::Loss(double angle_param){
  return Linear*angle_param + Constant;
}

PiecewiseLinearLossFun::PiecewiseLinearLossFun(){
  weight = 1;
}

void PiecewiseLinearLossFun::set_to_min_of_one
(PiecewiseLinearLossFun *input, int verbose){
  double best_loss = INFINITY,
    best_angle_param = INFINITY,
    prev_angle_param = INFINITY;
  int data_i;
  input->Minimize(&best_loss, &best_angle_param, &data_i, &prev_angle_param);
  piece_list.clear();
  piece_list.emplace_front(0, best_loss, 0, 2*PI, PREV_NOT_SET, best_angle_param);
}

void PiecewiseLinearLossFun::push_sum_pieces
(PiecewiseLinearLossFun *fun1,
 PiecewiseLinearLossFun *fun2,
 L1LossPieceList::iterator it1,
 L1LossPieceList::iterator it2,
 int verbose){
  double min_angle_param =
    (it1->min_angle_param < it2->min_angle_param) ?
    it2->min_angle_param : it1->min_angle_param;
  double max_angle_param =
    (it1->max_angle_param < it2->max_angle_param) ?
    it1->max_angle_param : it2->max_angle_param;
  enlarge_last_or_emplace
    (it1->Linear+it2->Linear, it1->Constant+it2->Constant,
     min_angle_param, max_angle_param, it2->data_i, it2->prev_angle_param);
}

void PiecewiseLinearLossFun::push_min_pieces
(PiecewiseLinearLossFun *fun1,
 PiecewiseLinearLossFun *fun2,
 L1LossPieceList::iterator it1,
 L1LossPieceList::iterator it2,
 int verbose){
  double min_angle_param =
    (it1->min_angle_param < it2->min_angle_param) ?
    it2->min_angle_param : it1->min_angle_param;
  double max_angle_param =
    (it1->max_angle_param < it2->max_angle_param) ?
    it1->max_angle_param : it2->max_angle_param;
  double diff_Constant = it1->Constant-it2->Constant;
  double diff_Linear = it1->Linear-it2->Linear;
  double root_angle_param = -diff_Constant/diff_Linear;
  if(diff_Linear == 0){
    if(diff_Constant < 0){
      // f1-f2<0 => f1<f2
      enlarge_last_or_emplace
	(it1, min_angle_param, max_angle_param);
    }else{
      enlarge_last_or_emplace
	(it2, min_angle_param, max_angle_param);
    }
  }else if(root_angle_param <= min_angle_param){
    if(diff_Linear < 0){
      // f1-f2<0 => f1<f2
      enlarge_last_or_emplace
	(it1, min_angle_param, max_angle_param);
    }else{
      // f1-f2>0 => f1>f2
      enlarge_last_or_emplace
	(it2, min_angle_param, max_angle_param);
    }
  }else if(root_angle_param >= max_angle_param){
    if(diff_Linear < 0){
      // f1-f2>0 => f1>f2
      enlarge_last_or_emplace
	(it2, min_angle_param, max_angle_param);
    }else{
      // f1-f2<0 => f1<f2
      enlarge_last_or_emplace
	(it1, min_angle_param, max_angle_param);
    }
  }else{
    //it1 intersects it2 between min and max -> add two pieces.
    if(diff_Linear < 0){
      //f1-f2>0 => f1>f2 before.
      enlarge_last_or_emplace
	(it2, min_angle_param, root_angle_param);
      enlarge_last_or_emplace
	(it1, root_angle_param, max_angle_param);
    }else{
      enlarge_last_or_emplace
	(it1, min_angle_param, root_angle_param);
      enlarge_last_or_emplace
	(it2, root_angle_param, max_angle_param);
    }      
  }
}

void PiecewiseLinearLossFun::set_to_min_of_two
(PiecewiseLinearLossFun *fun1,
 PiecewiseLinearLossFun *fun2,
 int verbose){
  while_piece_pairs(fun1, fun2, &PiecewiseLinearLossFun::push_min_pieces, verbose);
}

void PiecewiseLinearLossFun::set_to_sum_of
(PiecewiseLinearLossFun *fun1,
 PiecewiseLinearLossFun *fun2,
 int verbose){
  while_piece_pairs(fun1, fun2, &PiecewiseLinearLossFun::push_sum_pieces, verbose);
}

void PiecewiseLinearLossFun::while_piece_pairs
(PiecewiseLinearLossFun *fun1,
 PiecewiseLinearLossFun *fun2,
 push_fun_ptr push_pieces,
 int verbose){
  L1LossPieceList::iterator
    it1 = fun1->piece_list.begin(),
    it2 = fun2->piece_list.begin();
  piece_list.clear();
  while(it1 != fun1->piece_list.end() &&
	it2 != fun2->piece_list.end()){
    (this->*push_pieces)(fun1, fun2, it1, it2, verbose);
    double last_max_angle_param = piece_list.back().max_angle_param;
    if(it1->max_angle_param == last_max_angle_param){
      it1++;
    }
    if(it2->max_angle_param == last_max_angle_param){
      it2++;
    }
  }
}

void PiecewiseLinearLossFun::add(double Constant){
  L1LossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->Constant += Constant;
  }
}

void PiecewiseLinearLossFun::multiply(double x){
  L1LossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->Linear *= x;
    it->Constant *= x;
  }
}

void PiecewiseLinearLossFun::set_prev_seg_end(int prev_seg_end){
  L1LossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->data_i = prev_seg_end;
  }
}

void PiecewiseLinearLossFun::findMean
(double angle_param, int *seg_end, double *prev_angle_param){
  L1LossPieceList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_angle_param <= angle_param && angle_param <= it->max_angle_param){
      *seg_end = it->data_i;
      *prev_angle_param = it->prev_angle_param;
      return;
    }
  }
}

void PiecewiseLinearLossFun::print(){
  L1LossPieceList::iterator it;
  Rprintf("%5s %5s %5s %5s %5s %s\n",
	  "Linear", "Constant",
	  "min_angle_param", "max_angle_param",
	  "data_i", "prev_angle_param");
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->print();
  }
}

void LinearPiece::print(){
  Rprintf("%.20e %.20e %15f %15f %15f %d\n",
	 Linear, Constant,
	 min_angle_param, max_angle_param,
	 prev_angle_param, data_i);
}

void PiecewiseLinearLossFun::emplace_piece
(double Linear, double Constant,
 double min_angle_param, double max_angle_param){
  piece_list.emplace_back
    (Linear*weight, Constant*weight,
     min_angle_param, max_angle_param);
}

void PiecewiseLinearLossFun::emplace_piece
(double Linear, double Constant,
 double min_angle_param, double max_angle_param,
 int data_i, double prev_angle_param){
  piece_list.emplace_back
    (Linear*weight, Constant*weight,
     min_angle_param, max_angle_param,
     data_i, prev_angle_param);
}

void PiecewiseLinearLossFun::enlarge_last_or_emplace
(double Linear, double Constant,
 double min_angle_param, double max_angle_param,
 int data_i, double prev_angle_param){
  L1LossPieceList::iterator it=piece_list.end();
  if(it!=piece_list.begin()){
    it--;
    if(it->Linear == Linear && it->Constant == Constant &&
       it->data_i==data_i && it->prev_angle_param==prev_angle_param){
      it->max_angle_param = max_angle_param;
      return;
    }
  }
  emplace_piece
    (Linear*weight, Constant*weight,
     min_angle_param, max_angle_param,
     data_i, prev_angle_param);
}

void PiecewiseLinearLossFun::enlarge_last_or_emplace
(double Linear, double Constant,
 double min_angle_param, double max_angle_param){
  enlarge_last_or_emplace
    (Linear, Constant, min_angle_param, max_angle_param,
     PREV_NOT_SET, INFINITY);
}

void PiecewiseLinearLossFun::enlarge_last_or_emplace
(L1LossPieceList::iterator it,
 double min_angle_param, double max_angle_param){
  enlarge_last_or_emplace
    (it->Linear, it->Constant,
     min_angle_param, max_angle_param,
     it->data_i, it->prev_angle_param);
}

void PiecewiseLinearLossFun::init
(double angle, double weight_){
  weight = weight_;
  piece_list.clear();
  if(angle == 0){
    emplace_piece(1, 0, 0, PI);
    emplace_piece(-1, 2*PI, PI, 2*PI);
  }else if(angle < PI){
    emplace_piece(-1, angle, 0, angle);
    emplace_piece(1, -angle, angle, angle+PI);
    emplace_piece(-1, (2*PI+angle), angle+PI, 2*PI);
  }else if(angle == PI){
    emplace_piece(-1, PI, 0, PI);
    emplace_piece(1, -PI, PI, 2*PI);
  }else{
    emplace_piece(1, 2*PI-angle, 0, angle-PI);
    emplace_piece(-1, angle, angle-PI, angle);
    emplace_piece(1, -angle, angle, 2*PI);
  }
}

void PiecewiseLinearLossFun::Minimize
(double *best_loss,
 double *best_angle_param,
 int *data_i,
 double *prev_angle_param){
  for
    (L1LossPieceList::iterator it = piece_list.begin();
     it != piece_list.end();
     it ++){
    double it_angle_param =
      (it->Linear < 0) ? it->max_angle_param : it->min_angle_param;
    double it_loss = it->Loss(it_angle_param);
    if(it_loss < *best_loss){
      *best_loss = it_loss;
      *best_angle_param = it_angle_param;
      *data_i = it->data_i;
      *prev_angle_param = it->prev_angle_param;
    }
  }
}

int PiecewiseFunSize(const PiecewiseLinearLossFun&fun){
  int sizeof_piece = 2*sizeof(double) + sizeof(int);
  return sizeof_piece*fun.piece_list.size() +
    sizeof(int)*2; // n_pieces and chromEnd.
}

void PiecewiseFunCopy(void *dest, const PiecewiseLinearLossFun&fun){
  char *p = (char*)dest;
  int n_pieces = fun.piece_list.size();
  memcpy(p, &n_pieces, sizeof(int));
  p += sizeof(int);
  memcpy(p, &(fun.chromEnd), sizeof(int));
  p += sizeof(int);
  for(L1LossPieceList::const_iterator it = fun.piece_list.begin();
      it != fun.piece_list.end(); it++){
    memcpy(p, &(it->max_angle_param), sizeof(double));
    p += sizeof(double);
    memcpy(p, &(it->data_i), sizeof(int));
    p += sizeof(int);
    memcpy(p, &(it->prev_angle_param), sizeof(double));
    p += sizeof(double);
  }
}

void PiecewiseFunRestore(PiecewiseLinearLossFun&fun, const void *src){
  int n_pieces;
  char *p = (char*)src;
  LinearPiece piece;
  memcpy(&n_pieces, p, sizeof(int));
  p += sizeof(int);
  memcpy(&(fun.chromEnd), p, sizeof(int));
  p += sizeof(int);
  double min_angle_param = -INFINITY;
  for(int piece_i=0; piece_i < n_pieces; piece_i++){
    piece.min_angle_param = min_angle_param;
    memcpy(&(piece.max_angle_param), p, sizeof(double));
    p += sizeof(double);
    memcpy(&(piece.data_i), p, sizeof(int));
    p += sizeof(int);
    memcpy(&(piece.prev_angle_param), p, sizeof(double));
    p += sizeof(double);
    fun.piece_list.push_back(piece);
    min_angle_param = piece.max_angle_param;
  }
}

class UndefinedReadException : public std::exception {
  const char * what() const throw(){
    return "Attempt to read from undefined position";
  }
};

class AlreadyWrittenException : public std::exception {
  const char * what() const throw(){
    return "Attempt to write from already defined position";
  }
};

class WriteFailedException : public std::exception {
  const char * what() const throw(){
    return "Attempt to write to file failed";
  }
};

class DiskVector {
public:
  std::fstream db; //fstream supports both input and output.
  std::streampos beginning;
  int n_entries;
  void init(const char *filename, int N){
    n_entries = N;
    db.open(filename, std::ios::binary|std::ios::in|std::ios::out|std::ios::trunc);
    // reserve the first n_entries for streampos objects that will
    // tell us where to look for the data.
    beginning = db.tellp();
    for(int i=0; i<n_entries;i++){
      write_or_exception((char*)&beginning, sizeof(std::streampos));
    }
  }
  void write_or_exception(char * p, int size){
    db.write(p, size);
    if(db.fail()){
      throw WriteFailedException();
    }
  }
  void seek_element(int element){
    db.seekp(sizeof(std::streampos)*element, std::ios::beg);
  }
  std::streampos get_element_position(int element){
    seek_element(element);
    std::streampos pos;
    db.read((char*)&pos, sizeof(std::streampos));
    return pos;
  }
  PiecewiseLinearLossFun read(int element){
    std::streampos pos = get_element_position(element);
    if(pos == beginning){
      throw UndefinedReadException();
    }
    db.seekp(pos);
    int size;
    db.read((char*)&size, sizeof(int));
    void * buffer = malloc(size);
    db.read((char*)buffer, size);
    PiecewiseLinearLossFun fun;
    PiecewiseFunRestore(fun, buffer);
    free(buffer);
    return fun;
  }
  void write(int element, PiecewiseLinearLossFun fun){
    std::streampos pos = get_element_position(element);
    if(pos != beginning){
      throw AlreadyWrittenException();
    }
    // serialize at end of file.
    db.seekp(0, std::ios::end);
    pos = db.tellp();//save pos for later.
    // first write size of fun.
    int size = PiecewiseFunSize(fun);
    write_or_exception((char*)&size, sizeof(int));
    // then write the data itself.
    void * buffer = malloc(size);
    PiecewiseFunCopy(buffer, fun);
    write_or_exception((char*)buffer, size);
    free(buffer);
    // write position.
    seek_element(element);
    write_or_exception((char*)&pos, sizeof(std::streampos));
  }
};

int geodesicFPOP
(const char *bedGraph_file_name,
 const char *penalty_str,
 const char *db_file_name){
  bool penalty_is_Inf = strcmp(penalty_str, "Inf") == 0;
  double penalty;
  try{
    penalty = std::stod(penalty_str);
  }catch(const std::invalid_argument& e){
    return ERROR_PENALTY_NOT_NUMERIC;
  }
  // Handle penalty error cases before opening files.
  if(penalty_is_Inf){
    //ok.
  }else if(!std::isfinite(penalty)){
    return ERROR_PENALTY_NOT_FINITE;
  }else if(penalty < 0){
    return ERROR_PENALTY_NEGATIVE;
  }
  std::ifstream bedGraph_file(bedGraph_file_name);
  if(!bedGraph_file.is_open()){
    return ERROR_UNABLE_TO_OPEN_BEDGRAPH;
  }
  std::string line;
  int chromStart, chromEnd, items, line_i=0;
  double angle;
  char chrom[100];
  char extra[100] = "";
  double cum_weight_i = 0.0, cum_weight_prev_i=-1.0;
  int data_i = 0;
  double weight;
  int first_chromStart=-1, prev_chromEnd=-1;
  while(std::getline(bedGraph_file, line)){
    line_i++;
    items = sscanf
      (line.c_str(),
       "%s %d %d %lf%s\n",
       chrom, &chromStart, &chromEnd, &angle, extra);
    //Rprintf("%s %d %d %f%s\n", chrom, chromStart, chromEnd, angle, extra);
    if(items < 4){
      Rprintf("problem: %d items on line %d\n", items, line_i);
      return ERROR_NOT_ENOUGH_COLUMNS;
    }
    if(0 < strlen(extra)){
      return ERROR_NON_INTEGER_DATA;
    }
    weight = chromEnd-chromStart;
    cum_weight_i += weight;
    if(line_i == 1){
      first_chromStart = chromStart;
    }else{
      if(chromStart != prev_chromEnd){
	return ERROR_INCONSISTENT_CHROMSTART_CHROMEND;
      }
    }
    prev_chromEnd = chromEnd;
  }
  int data_count = line_i;
  if(data_count==0){
    return ERROR_NO_DATA;
  }
  double best_cost = INFINITY,
    best_angle_param = INFINITY,
    prev_angle_param = INFINITY;
  // open segments and loss files for writing.
  std::string penalty_prefix = bedGraph_file_name;
  penalty_prefix += "_penalty=";
  penalty_prefix += penalty_str;
  std::string segments_file_name = penalty_prefix + "_segments.bed";
  std::string loss_file_name = penalty_prefix + "_loss.tsv";
  std::ofstream segments_file, loss_file; // ofstream supports output only.
  // Opening both files here is fine even if we error exit, because
  // "any open file is automatically closed when the ofstream object
  // is destroyed."
  // http://www.cplusplus.com/reference/fstream/ofstream/close/
  loss_file.open(loss_file_name.c_str());
  segments_file.open(segments_file_name.c_str());
  bedGraph_file.clear();
  bedGraph_file.seekg(0, std::ios::beg);
  DiskVector cost_model_mat;
  try{
    cost_model_mat.init(db_file_name, data_count);
  }catch(WriteFailedException& e){
    return ERROR_WRITING_COST_FUNCTIONS;
  }
  PiecewiseLinearLossFun dist_fun_i, cost_up_to_i, cost_up_to_prev, cost_of_change, min_term;
  int verbose=0;
  cum_weight_i = 0;
  double total_intervals = 0.0, max_intervals = 0.0;
  while(std::getline(bedGraph_file, line)){
    items = sscanf(line.c_str(), "%*s\t%d\t%d\t%lf\n", &chromStart, &chromEnd, &angle);
    //Rprintf("data_i=%d start=%d end=%d angle=%f\n", data_i, chromStart, chromEnd, angle);
    weight = chromEnd-chromStart;
    cum_weight_i += weight;
    dist_fun_i.init(angle, weight);
    if(data_i==0){
      cost_up_to_i = dist_fun_i;
    }else{
      if(penalty_is_Inf){
	min_term = cost_up_to_prev;
      }else{
	cost_of_change.set_to_min_of_one(&cost_up_to_prev, verbose);
	// V_t(m) = (gamma_t + w_{1:t-1} * M_t(m))/w_{1:t}, where
	// M_t(m) = min{
	//   V_{t-1}(m),
	//   Vbar_{t-1} + penalty/w_{1:t-1}
	// in other words, we need to divide the penalty by the previous cumsum,
	// and add that to the min-less-ified function, before applying the min-env
	cost_of_change.set_prev_seg_end(data_i-1);
	cost_of_change.add(penalty/cum_weight_prev_i);
	if(penalty==0){
	  min_term = cost_of_change;
	}else{
	  min_term.set_to_min_of_two(&cost_of_change, &cost_up_to_prev, verbose);
	}
      }
      min_term.multiply(cum_weight_prev_i);
      cost_up_to_i.set_to_sum_of(&dist_fun_i, &min_term, verbose);
    }
    cost_up_to_i.multiply(1/cum_weight_i);
    cum_weight_prev_i = cum_weight_i;
    total_intervals += cost_up_to_i.piece_list.size();
    if(max_intervals < cost_up_to_i.piece_list.size()){
      max_intervals = cost_up_to_i.piece_list.size();
    }
    cost_up_to_prev = cost_up_to_i;
    cost_up_to_i.chromEnd = chromEnd;
    try{
      cost_model_mat.write(data_i, cost_up_to_i);
    }catch(WriteFailedException& e){
      return ERROR_WRITING_COST_FUNCTIONS;
    }
    data_i++;
  }//while(can read line in text file)
  // Decoding the cost_model_vec, and writing to the output matrices.
  int prev_seg_end = -10;
  cost_up_to_i.Minimize
    (&best_cost, &best_angle_param,
     &prev_seg_end, &prev_angle_param);
  //Rprintf("param=%f end_i=%d chromEnd=%d\n", best_angle_param, prev_seg_end, cost_up_to_i.chromEnd);
  prev_chromEnd = cost_up_to_i.chromEnd;
  line_i=1;
  while(0 <= prev_seg_end){
    line_i++;
    cost_up_to_i = cost_model_mat.read(prev_seg_end);
    segments_file << chrom << "\t" << cost_up_to_i.chromEnd << "\t" << prev_chromEnd << "\tUNUSED\t" << best_angle_param << "\n";
    prev_chromEnd = cost_up_to_i.chromEnd;
    best_angle_param = prev_angle_param;
    cost_up_to_i.findMean
      (best_angle_param, &prev_seg_end, &prev_angle_param);
    //Rprintf("param=%f end=%d chromEnd=%d\n", best_angle_param, prev_seg_end, up_cost.chromEnd);
  }//for(data_i
  segments_file << chrom << "\t" << first_chromStart << "\t" << prev_chromEnd << "\tUNUSED\t" << best_angle_param << "\n";
  double total_penalty = (line_i==1) ? 0 : penalty*(line_i-1);
  loss_file << std::setprecision(20) << penalty << //penalty constant
    "\t" << line_i << //segments
    "\t" << (int)cum_weight_i << //total weight
    "\t" << data_count << //bedGraph_lines
    "\t" << best_cost << //mean penalized cost
    "\t" << best_cost*cum_weight_i-total_penalty << //total un-penalized cost
    "\t" << total_intervals/data_count << //mean intervals
    "\t" << max_intervals <<
    "\n";
  if(loss_file.fail()){
    return ERROR_WRITING_LOSS_OUTPUT;
  }
  if(segments_file.fail()){
    return ERROR_WRITING_SEGMENTS_OUTPUT;
  }
  return 0;
}

