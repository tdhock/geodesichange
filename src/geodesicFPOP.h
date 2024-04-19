#define ERROR_PENALTY_NOT_FINITE 1
#define ERROR_PENALTY_NEGATIVE 2
#define ERROR_UNABLE_TO_OPEN_BEDGRAPH 3
#define ERROR_NOT_ENOUGH_COLUMNS 4
#define ERROR_NON_INTEGER_DATA 5
#define ERROR_INCONSISTENT_CHROMSTART_CHROMEND 6
#define ERROR_WRITING_COST_FUNCTIONS 7
#define ERROR_WRITING_LOSS_OUTPUT 8
#define ERROR_NO_DATA 9
#define ERROR_PENALTY_NOT_NUMERIC 10
#define ERROR_WRITING_SEGMENTS_OUTPUT 11

#include <list>

int geodesicFPOP(const char *, const char *, const char *, const char *);

class LinearPiece {
 public:
  double Linear;
  double Constant;
  double min_angle_param;
  double max_angle_param;
  int data_i;
  double prev_angle_param;
  LinearPiece();
  double Loss(double);
  LinearPiece(double li, double co, double m, double M);
  LinearPiece
    (double li, double co, double m, double M, int i, double);
  void print();
};

typedef std::list<LinearPiece> L1LossPieceList;

class PiecewiseLinearLossFun;

typedef void (PiecewiseLinearLossFun::*push_fun_ptr)
(PiecewiseLinearLossFun*,
 PiecewiseLinearLossFun*,
 L1LossPieceList::iterator,
 L1LossPieceList::iterator,
 int);

class PiecewiseLinearLossFun {
 public:
  L1LossPieceList piece_list;
  int chromEnd;
  double weight;
  PiecewiseLinearLossFun();
  void push_sum_pieces(PiecewiseLinearLossFun*, PiecewiseLinearLossFun*, L1LossPieceList::iterator, L1LossPieceList::iterator, int);
  void push_min_pieces(PiecewiseLinearLossFun*, PiecewiseLinearLossFun*, L1LossPieceList::iterator, L1LossPieceList::iterator, int);
  void while_piece_pairs(PiecewiseLinearLossFun*, PiecewiseLinearLossFun*, push_fun_ptr, int);
  void emplace_piece(double,double,double,double);
  void emplace_piece(double,double,double,double,int,double);
  void enlarge_last_or_emplace(L1LossPieceList::iterator,double,double);
  void enlarge_last_or_emplace(double,double,double,double);
  void enlarge_last_or_emplace(double,double,double,double,int,double);
  void init(double, double);
  void set_to_min_of_one(PiecewiseLinearLossFun *, int);
  void set_to_min_of_two(PiecewiseLinearLossFun *, PiecewiseLinearLossFun *, int);
  void set_to_sum_of(PiecewiseLinearLossFun *, PiecewiseLinearLossFun *, int);
  void add(double Constant);
  void multiply(double);
  void print();
  void set_prev_seg_end(int prev_seg_end);
  void findMean(double mean, int *seg_end, double *prev_angle_param);
  double findCost(double mean);
  void Minimize
    (double *best_cost,
     double *best_mean,
     int *data_i,
     double *prev_angle_param);
};


