#ifndef __SCIP_SEPA_UTILITY_H__
#define __SCIP_SEPA_UTILITY_H__

#include <list>
#include <set>
#include <tuple>
#include <vector>
#include <algorithm>
#include <ctime>
#include "objscip/objscip.h"
#include "scip/cons_and.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"

using namespace scip;
using namespace std;

#ifdef DEBUG_TIME
static SCIP_Real memory_time;
static SCIP_Real compute_time;
static int ct_zero;
static int ct_full;
static SCIP_Real sort_time, coef_time, loop_time, inner_loop_time;
static time_t t1;
static time_t t2;
static time_t t3;
static time_t t4; 
static time_t t5;
static time_t t6; 
static time_t t8;
static time_t t9; 
static time_t t10; 
static time_t t11; 
#endif

/*
 * multilinear function
 */

// f(x) = f_1(x) - f_2(x) + b \le t, t is obj_var
// we consider S = {(x,t):f_1(x) - f_2(x) + b - t \le 0}, and S^c= {(x,t):f_1(x) - f_2(x) + b - t \ge 0}
// f_1 all negative weights (submodular), f_2 all negative weights (submodular).

class MLF{
public:
    // original vars
	vector<int> probvartypes; // problem varibale types, 0: continuous, 1: f1, 2:f2, 3: both
	vector<SCIP_VAR *> probvars; // problem variables x
	SCIP_Real constant;  // b
	int size; 


	// linear parts
	vector<pair<int, SCIP_Real>> linear_term;	// linear terms: index and coefficient (including obj_var)
	int ct_ind;


	// mulitlinear terms 
	vector<pair<list<int>, SCIP_Real>> mlf1_term;	// f_1: multilinear terms, negative weights
	vector<tuple< unsigned short int, unsigned short int, SCIP_Real>> mlf1_termval; // sum, degree, coefficient
	vector<list<int>> f1termofvar;
	vector<bool> mlf1_termone;


	vector<pair<list<int>, SCIP_Real>> mlf2_term;	// f_2: multilinear terms, negative weights
	vector<tuple< unsigned short int, unsigned short int, SCIP_Real>> mlf2_termval; // sum, degree, coefficient
	vector<list<int>> f2termofvar;
	vector<bool> mlf2_termone;


	// lifting vars (note that liftvars have different index system, which is matched with multilinear term index.)
	vector<SCIP_VAR *> liftf1vars;
	int f1size;
	vector<SCIP_VAR *> liftf2vars;
	int f2size;
	vector<SCIP_VAR *> extendedf1vars; // the first part is probvars variables, and the second part is lift1vars
	vector<SCIP_VAR *> extendedf2vars; // the first part is probavars variables, and the second part is lift2vars
	
	// recording variable
	SCIP_Real val; // temp val of the whole function
	SCIP_Real linearobj, linearray; // sum of obj and ray of linear terms; 
	SCIP_Real linearf1obj, linearf1ray; // sum of obj and ray of linearized f1 term; 

	// working variables
	vector<SCIP_Real> sol; // temp solution values 
	vector<SCIP_Real> solf1, solf2; // temp solution values 
	vector<SCIP_Real> ray; // temp ray values 
	vector<SCIP_Real> linearize1; // linearize all the probvars
	vector<SCIP_Real> linearize2; // linearize all the probvars
	vector<SCIP_Real> w_indicator; // binarize all the probvars
	vector<pair<SCIP_Real, int>> w_var1linds; // given value and index of f1 
	vector<pair<SCIP_Real, int>> w_var2linds; // given value and index of f2

	vector<int> set1, set2; // the set of indices in f1, f2
	set<int> uniqeset1, uniqeset2;
	
	// tmp reference
	SCIP_Bool typef;
	int fsize = f1size;
	SCIP_Real prevOracle = 0;
	vector<int> *  fset = NULL;
	vector<SCIP_Real> *  linearize = NULL;
	vector<SCIP_VAR *> * liftfvars = NULL;
	vector<pair<SCIP_Real, int>>  * w_varlinds = NULL;
	vector<pair<list<int>, SCIP_Real>> * mlf_term = NULL;
	vector<SCIP_Real> * solf = NULL;
	vector<SCIP_VAR *> * extendedfvars = NULL;

	bool fast_oracle;

   	MLF(){};


   	void init();

   // f_1(x) or f_2(x) : supported for warm start by prevOracle and mlf_termone
   SCIP_Real oracleSepValueFast(
      const list<int> & termofvar,
      SCIP_Bool mlftype
   );

   // f_1(x) or f_2(x) : supported for warm start by prevOracle and mlf_termone
   SCIP_Real oracleSepValue(
      const vector<SCIP_Real> & mlf_term_vals,
      SCIP_Bool mlftype
   );



   // evaluate the underestimator of submodular function with the linear underestimator
   SCIP_Real evalUnderSubWith(
      vector<pair<SCIP_Real, int>> & w_valinds,    /**< variable value*/
      SCIP_Bool mlftype
   );

   /**< consider: f_1(x) - f_2(x) + b - t <= 0, evaluate the ray with the gradient information, isfull: no linearization, t == 0, and compute the gradient */
   SCIP_Real evalRayWith(
      const SCIP_Real & t,
      SCIP_Real & gradient,
      SCIP_Bool isfull
   );

   inline SCIP_Real evalRayWithZero(
   ){
	   	SCIP_Real gradient;
		return  evalRayWith(0, gradient, true);	   
   };

	// w liftfi(x) = fi(x), if fi(x) - w liftfi(x) < 0, we first want to enforce  fi(x) - w liftfi(x) >= 0 by intersection cut, 
	// otherwise, we use outer approximation cut to enforce: fi(x)  - w liftfi(x) <= 0,
   SCIP_Real evalRayfWithZero(
		SCIP_Bool typef_
   );
	// w liftfi(x) = fi(x), if fi(x) - w liftfi(x) < 0, we first want to enforce  fi(x) - w liftfi(x) >= 0 by intersection cut, 
	// otherwise, we use outer approximation cut to enforce: fi(x)  - w liftfi(x) <= 0,
   SCIP_Real evalRayfWith(
      const SCIP_Real & t,
      SCIP_Real & gradient
   );
   
};

struct SCIP_SepaData {
   	int    maxrounds;             /**< maximal number of separation rounds per node (-1: unlimited) */
   	int    maxroundsroot;         /**< maximal number of separation rounds in the root node (-1: unlimited) */
   	SCIP_Real mincutviol;            /**< minimum required violation of a cut */
	int   BINSEARCH_MAXITERS = 500; /**< default iteration limit for binary search */
	SCIP_Bool disaggregate; /**< disaggregate the cutting planes */
	MLF mlf;
};

void createMLFCons(
   MLF & mlf,
   SCIP_Real sideval,
   vector<SCIP_Real> & linear_term_coeffs,
   vector<SCIP_VAR *>  & linear_term_vars,
   vector<SCIP_Real> & mlf_term_coeffs,
   vector<vector<SCIP_VAR *>> & mlf_term_vars,
   SCIP_Real & constant,
   unordered_map<SCIP_VAR* , vector<SCIP_VAR*>> & multilinear_map,
   SCIP_Bool l
);

SCIP_RETCODE detectMLFCons(
   SCIP*                 scip,             /**< SCIP data structure */
   MLF &   mlf
);
#endif
