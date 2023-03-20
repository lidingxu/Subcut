#include "utility.h"


void MLF::init(){
	fast_oracle  = true;
	size = probvars.size();
	w_indicator = vector<SCIP_Real>(size);
	linearize1 = vector<SCIP_Real>(size);
	linearize2 = vector<SCIP_Real>(size);
	sol = vector<SCIP_Real>(size);
	ray = vector<SCIP_Real>(size);

	//set<int> uniqeset1;
	mlf1_termval = vector<tuple<unsigned short int, unsigned short int, SCIP_Real>> (mlf1_term.size());
	f1termofvar = vector<list<int>>(size);
	int termf1_id = 0;
	for(auto & p: mlf1_term){
		for(int ind: p.first){
			if(uniqeset1.find(ind) == uniqeset1.end()){
				uniqeset1.insert(ind);
				set1.push_back(ind);
			}
			f1termofvar[ind].push_back(termf1_id);
		}
		mlf1_termval[termf1_id] = make_tuple(0, p.first.size(), p.second);
		termf1_id++;
	}
	sort(set1.begin(), set1.end());
	mlf1_termone= vector<bool>(mlf1_term.size());

	//set<int> uniqeset2;
	mlf2_termval = vector<tuple<unsigned short int, unsigned short int, SCIP_Real>> (mlf2_term.size());
	f2termofvar = vector<list<int>>(size);
	int termf2_id = 0;
	for(auto & p: mlf2_term){
		for(int ind: p.first){
			if(uniqeset2.find(ind) == uniqeset2.end()){
				uniqeset2.insert(ind);
				set2.push_back(ind);
			}
			f2termofvar[ind].push_back(termf2_id);
		}
		mlf2_termval[termf2_id] = make_tuple(0, p.first.size(), p.second);
		termf2_id++;;
	}
	sort(set2.begin(), set2.end());
	mlf2_termone= vector<bool>(mlf2_term.size());

	for(int ind = 0; ind < size; ind++){
		bool inset1 =  uniqeset1.find(ind) != uniqeset1.end();
		bool inset2 =  uniqeset2.find(ind) != uniqeset2.end();
		if(inset1 && inset2){
			probvartypes[ind] = 3;
		}
		else if(inset1){
			probvartypes[ind] = 1;
		}
		else if(inset2){
			probvartypes[ind] = 2;
		}
	}

	w_var1linds = vector<pair<SCIP_Real, int>>(set1.size());

	w_var2linds = vector<pair<SCIP_Real, int>>(set2.size());

	// find continuous var (obj_var) index
	for(int i = 0; i < size; i++){
		ct_ind = probvartypes[i] ?  ct_ind: i;
	}


	// lifting variables
	f1size =  liftf1vars.size();
	f2size =  liftf2vars.size();
	solf1 = vector<SCIP_Real>(f1size);
	solf2 = vector<SCIP_Real>(f2size);

	prevOracle = 0;
	// exnteded vars
	extendedf1vars = vector<SCIP_VAR *>(probvars);
	extendedf1vars.insert(extendedf1vars.end(), liftf1vars.begin(), liftf1vars.end());

	SCIPdebugMessage("%d %d %d\n", probvars.size(), liftf1vars.size(), liftf2vars.size());

	//SCIPdebugMessage("192\n");
	extendedf2vars = vector<SCIP_VAR *>(probvars);
	//SCIPdebugMessage("193\n");
	extendedf2vars.insert(extendedf2vars.end(), liftf2vars.begin(), liftf2vars.end());
	//SCIPdebugMessage("194\n");
}

// f_1(x) or f_2(x) : supported for warm start by prevOracle and mlf_termone
SCIP_Real MLF::oracleSepValueFast(
	const list<int> & termofvar,
	SCIP_Bool mlftype
){
	SCIP_Real result = prevOracle;
	// f1 or f2
	auto & mlf_termval = mlftype ? mlf1_termval: mlf2_termval;
	//SCIPdebugMessage("%d %d %d", mlftype, mlf1_term.size(), mlf2_term.size());
	// only consider zeros terms
	#ifdef DEBUG_TIME
	time(&t1);
	#endif
	for(int term_id: termofvar){
		auto & termval = mlf_termval[term_id];
		// boolean product
		#ifdef DEBUG_TIME
		time(&t3);
		#endif
		get<0>(termval) += 1;
		result = (get<0>(termval) != get<1>(termval)) ? result : (result + get<2>(termval) );
		#ifdef DEBUG_TIME
		time(&t4);
		inner_loop_time += difftime(t4,t3);
		#endif
	}
	#ifdef DEBUG_TIME
	time(&t2);
	loop_time += difftime(t2,t1);
	#endif
	prevOracle = result;
	return result;
}

// f_1(x) or f_2(x) : supported for warm start by prevOracle and mlf_termone
SCIP_Real MLF::oracleSepValue(
	const vector<SCIP_Real> & mlf_term_vals,
	SCIP_Bool mlftype
){
	SCIP_Real result = prevOracle;
	// f1 or f2
	auto & mlf_term = mlftype ? mlf1_term: mlf2_term;
	auto & mlf_termone = mlftype ? mlf1_termone: mlf2_termone;
	//SCIPdebugMessage("%d %d %d", mlftype, mlf1_term.size(), mlf2_term.size());
	int numterm = mlf_term.size();
	// only consider zeros terms
	#ifdef DEBUG_TIME
	time(&t1);
	#endif
	for(int i = 0; i < numterm; i++){
		if(mlf_termone[i]){
			continue;
		}
		auto & p = mlf_term[i];
		bool allone = true;
		// boolean product
		#ifdef DEBUG_TIME
		time(&t3);
		#endif
		for(int ind: p.first){
			if(mlf_term_vals[ind] < 0.5){
				allone = false;
				break;
			}
		}
		#ifdef DEBUG_TIME
		time(&t4);
		inner_loop_time += difftime(t4,t3);
		#endif
		mlf_termone[i] = allone;
		result += p.second * allone;
	}
	#ifdef DEBUG_TIME
	time(&t2);
	loop_time += difftime(t2,t1);
	#endif
	prevOracle = result;
	return result;
}



// evaluate the underestimator of submodular function with the linear underestimator
SCIP_Real MLF::evalUnderSubWith(
	vector<pair<SCIP_Real, int>> & w_valinds,    /**< variable value*/
	SCIP_Bool mlftype
){
	// sorting
	#ifdef DEBUG_TIME
	time(&t8);
	#endif
	sort(w_valinds.begin(), w_valinds.end(), [](const pair<SCIP_Real, int> & p1, const  pair<SCIP_Real, int> & p2){ return p1.first > p2.first || (p1.first == p2.first && p1.second > p2.second);});
	// fill indicator
	#ifdef DEBUG_TIME
	time(&t9);
	sort_time += difftime(t9, t8);
	#endif
	// get number of values
	int numval = w_valinds.size();
	// function values
	SCIP_Real prev_val = 0;
	SCIP_Real underobj = 0;
	prevOracle = 0;
	// get linearize term and set to 0
	auto & linearize = mlftype ? linearize1 : linearize2;
	fill(linearize.begin(), linearize.end(), 0.);
	if(fast_oracle)
	{
		auto & mlf_termval = mlftype ? mlf1_termval: mlf2_termval;
		auto & ftermofvar = mlftype ? f1termofvar: f2termofvar;
		int term_size = mlf_termval.size();
		for(int i = 0; i < term_size; i++){
			get<0>(mlf_termval[i]) = 0;
		}
		for(int i = 0; i < numval; i++){
			int id  = w_valinds[i].second;
			SCIP_Real current_val = oracleSepValueFast(ftermofvar[id], mlftype);
			SCIP_Real coeff = current_val - prev_val;
			linearize[id] = coeff;
			prev_val = current_val;
			underobj += coeff * w_valinds[i].first;
		}
	}
	else{
		fill(w_indicator.begin(), w_indicator.end(), 0.);
		auto & mlf_termone = mlftype ? mlf1_termone: mlf2_termone;
		fill(mlf_termone.begin(), mlf_termone.end(), 0.);
		for(int i = 0; i < numval; i++){
			int id  = w_valinds[i].second;
			w_indicator[id] = 1;
			SCIP_Real current_val = oracleSepValue(w_indicator, mlftype);
			SCIP_Real coeff = current_val - prev_val;
			linearize[id] = coeff;
			prev_val = current_val;
			underobj += coeff * w_valinds[i].first;
		}
	}
	return underobj;
}

/**< consider: f_1(x) - f_2(x) + b - t <= 0, evaluate the ray with the gradient information, isfull: no linearization, t == 0, and compute the gradient */
SCIP_Real MLF::evalRayWith(
	const SCIP_Real & t,
	SCIP_Real & gradient,
	SCIP_Bool isfull
){
	// start with constant
	SCIP_Real underobj = constant;
	// init gradient
	gradient = 0;
	SCIP_Real f1 = 0, f2 = 0;
	if(isfull){ // no linearization
		// set variables of f1, f2
		int  i = 0;
		for(int id: set1){
			w_var1linds[i] = make_pair(sol[id], id);
			i++;
		}
		i = 0;
		for(int id: set2){
			w_var2linds[i] = make_pair(sol[id], id);
			i++;
		}
		// compute f1, f2	
		f1 = evalUnderSubWith(w_var1linds, TRUE);
		f2 = evalUnderSubWith(w_var2linds, FALSE);		
		// compute the linear part
		// linearized f1
		linearf1obj = 0;	
		for(int id: set1){
			linearf1obj += linearize1[id] * (sol[id]);
		}
		// compute the linear part
		linearobj = 0;
		for(auto p: linear_term){
			SCIP_Real value =  sol[p.first] * p.second;
			underobj += value;
			linearobj += value;
		}
		underobj += f1 - f2; //evalUnderSubWith(w_var1linds, TRUE) - evalUnderSubWith(w_var2linds, FALSE);
		val = underobj;
		#ifdef SCIP_DEBUG_
		SCIPdebugMessage("%f,  %f\n", f1, f2);		
		#endif
	}
	else{
		// compute f2
		// sort and compute
		int i = 0;
		for(int id: set2){
			w_var2linds[i] = make_pair(sol[id] + t * ray[id], id);
			i++;
		}
		f2 = evalUnderSubWith(w_var2linds, FALSE);	
		// compute f2 gradient
		for(int id: set2){
			gradient -= linearize2[id] * ray[id];
			//SCIPdebugMessage("%f ", linearize2[id]);
		}
		// f1
		f1 = linearf1obj + t * linearf1ray;
		gradient += linearf1ray;
		// f1 - f2
		underobj += f1 - f2;
		// f1 -f2 + t
		underobj += linearobj +  t * linearray;
		gradient += linearray;
	}    
	//printf("%f %f\n", test_val, underobj);
	return underobj;	
}



// w liftfi(x) = fi(x), if fi(x) - w liftfi(x) < 0, we first want to enforce  fi(x) - w liftfi(x) >= 0 by intersection cut, 
// otherwise, we use outer approximation cut to enforce: fi(x)  - w liftfi(x) <= 0,
SCIP_Real MLF::evalRayfWithZero(
	SCIP_Bool typef_
){
	// set data
	typef = typef_;
	fset = typef ? &set1 : &set2;
	linearize = typef ? &linearize1 : &linearize2;
	liftfvars = typef ? &liftf1vars : &liftf2vars;
	w_varlinds = typef ? &w_var1linds : &w_var2linds;
	fsize = typef ? f1size : f2size;
	mlf_term = typef ? &mlf1_term : &mlf2_term;
	solf = typef ? &solf1 : &solf2;
	extendedfvars = typef ? &extendedf1vars : &extendedf2vars;
	int  i = 0;
	val = 0;
	// prepare
	for(int id: *fset){
		(*w_varlinds)[i] = make_pair(sol[id], id);
		i++;
	}
	// sorting and compute fi(x)
	SCIP_Real f = evalUnderSubWith(*w_varlinds, typef);
	// w liftfi(x)
	linearobj = 0;
	for(int i = 0; i < fsize; i++){
		linearobj += (*mlf_term)[i].second * (*solf)[i];
	}
	//  w liftfi(x) - fi(x)
	SCIP_Real underobj = f - linearobj;
	// tmp val
	val = underobj;
	return underobj;   
}

// w liftfi(x) = fi(x), if fi(x) - w liftfi(x) < 0, we first want to enforce  fi(x) - w liftfi(x) >= 0 by intersection cut, 
// otherwise, we use outer approximation cut to enforce: fi(x)  - w liftfi(x) <= 0,
SCIP_Real MLF::evalRayfWith(
	const SCIP_Real & t,
	SCIP_Real & gradient
){
	gradient = 0;
	int i = 0;
	// prepare
	for(int id: *fset){
		(*w_varlinds)[i] = make_pair(sol[id] + t * ray[id], id);
		i++;
	}
	// fi
	SCIP_Real f = evalUnderSubWith(*w_varlinds, typef);
	// compute fi gradient
	for(int id: *fset){
		gradient += (*linearize)[id] * ray[id];
	}
	// liftfi(w) - f(x)
	SCIP_Real underobj =  f -  (linearobj +  t * linearray);
	gradient -= linearray;
	return underobj;
}


// create Boolean mulilinear quadratic constraint
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
){
   printf("mlf to  creat\n");
   unordered_map<SCIP_VAR *, int> var2ind;
   int varind = 0;
   SCIP_Real multiply = l ? 1 : -1;
   mlf.constant = multiply * (sideval - constant);
   for(int i = 0; i < linear_term_coeffs.size(); i++){
      if(linear_term_coeffs[i] != 0){
         if(var2ind.find(linear_term_vars[i]) == var2ind.end()){
            var2ind[linear_term_vars[i]] = varind;
            varind++;
         }
         mlf.linear_term.push_back( make_pair(var2ind[linear_term_vars[i]], - multiply* linear_term_coeffs[i]));
      }
   }
   for(int i = 0; i < mlf_term_coeffs.size(); i++){
      if(mlf_term_coeffs[i] == 0){
         continue;
      }
      //SCIPdebugMessage("%d\n", mlf_term_vars[i].size());
      list<int> indlst;
      //SCIPdebugMessage("%d %d\n", indlst.size(), mlf_term_vars[i].size());
      for(auto & var: mlf_term_vars[i]){
         if( var2ind.find(var) == var2ind.end() ){
            var2ind[var] = varind;
            indlst.push_back(varind);
            varind++;
         }
         else{
            indlst.push_back(var2ind[var]);
         }
      }
      // match the lifting variable
      SCIP_VAR * liftvar = NULL;
      for(auto & ele: multilinear_map){
         if(ele.second.size() != mlf_term_vars[i].size()){
            continue;
         }
         auto p1 =  ele.second.begin();
         auto p2 = mlf_term_vars[i].begin();
         while(p1 != ele.second.end() && *p1 == *p2){
            p1++;
            p2++;
         }
         if(p1 == ele.second.end()){ // matched
            liftvar = ele.first;
            break;
         }
      }
      assert(liftvar != NULL);
      if(mlf_term_coeffs[i] * multiply  > 0){
         //p.quad1_term_coeffs.push_back(-multiple*quad_term_coeffs[i]);
         //SCIPdebugMessage("%d\n", indlst.size());
         mlf.liftf1vars.push_back(liftvar);
         mlf.mlf1_term.push_back(make_pair(indlst, -multiply*mlf_term_coeffs[i]));
      }
      if(mlf_term_coeffs[i] * multiply < 0){
         //p.quad2_term_coeffs.push_back(multiply * quad_term_coeffs[i]);
         //SCIPdebugMessage("%d\n", indlst.size());
         mlf.liftf2vars.push_back(liftvar);
         mlf.mlf2_term.push_back(make_pair(indlst, multiply*mlf_term_coeffs[i]) );
      }
   }
   mlf.probvars = vector<SCIP_VAR *> (var2ind.size());
   mlf.probvartypes = vector<int> (var2ind.size(), 0);
   for(auto & mlf_var: var2ind){
      //SCIPdebugMessage("%d\n", SCIPvarIsBinary( mlf_var.first));
      mlf.probvars[mlf_var.second] = mlf_var.first;
      mlf.probvartypes[mlf_var.second] = SCIPvarIsBinary(mlf_var.first);// && multilinear_map.find(mlf_var.first) == multilinear_map.end();
   }
   // verification
   SCIP_Bool isprint = false;
   if(isprint){
      //printf("\n f1: %d %d\n", p.quad1_term.size(), p.quad2_term.size());
      int ct = 0;
      printf("\n f1:");
      for(auto & p: mlf.mlf1_term){
         printf("%f",  p.second);
         for(auto var: p.first){
            printf("%s", SCIPvarGetName(mlf.probvars[var]));
         }
         printf(" ");
         //ct++;
      }
      printf("\n");
      printf("\n f2:");
      for(auto & p: mlf.mlf2_term){
         printf("%f",  p.second);
         for(auto var: p.first){
            printf("%s", SCIPvarGetName(mlf.probvars[var]));
         }
         printf(" ");
         //ct++;
      }
      printf("\n");
      printf("\n linear:");
      for(auto pr: mlf.linear_term){
         printf("%f%s ", pr.second, SCIPvarGetName(mlf.probvars[pr.first]) );
      }
      printf("\n");
   }
   printf("mlf to init\n");
   mlf.init();
   printf("mlf inited\n");
}

SCIP_RETCODE detectMLFCons(
   SCIP*                 scip,             /**< SCIP data structure */
   MLF &   mlf
){

   int ncons = SCIPgetNConss(scip);
   //SCIPdebugMessage("%d\n", ncons);


   /* check whether there are and constraints available */
   SCIP_CONSHDLR*  conshdlrand = SCIPfindConshdlr(scip, "and");
   SCIP_CONSHDLR* conshdlrknap = SCIPfindConshdlr(scip, "knapsack");
   SCIP_CONSHDLR* conshdlrlinear = SCIPfindConshdlr(scip, "linear");

   //SCIPdebugMessage("s0 %d %d\n", SCIPconshdlrGetNConss(conshdlrand), SCIPconshdlrGetNConssFound(conshdlrand));
   if( conshdlrand == NULL || SCIPconshdlrGetNConss(conshdlrand) == 0 )
      return SCIP_OKAY;
   
   /* find all Boolean bilinear terms from and constraints  */
   unordered_map<SCIP_VAR* , vector<SCIP_VAR*>> multilinear_map;
   //SCIPdebugMessage("s1 %d\n", SCIPconshdlrGetNConss(conshdlrand));
   for(int c = 0; c < SCIPconshdlrGetNConss(conshdlrand); ++c )
   {
      SCIP_CONS* cons;

      cons = SCIPconshdlrGetConss(conshdlrand)[c];
      assert(cons != NULL);
      int nvarsand = SCIPgetNVarsAnd(scip, cons);
      SCIP_VAR ** vars = SCIPgetVarsAnd (scip, cons);
      SCIP_VAR * result_var = SCIPgetResultantAnd(scip, cons); 
      //SCIPdebugMessage("%d\n", nvarsand);	
      vector<SCIP_VAR *> list_vars(vars, vars + nvarsand);
      multilinear_map[result_var]= list_vars;
      //bilinear_map[result_var]= make_pair(vars[1], vars[0]);
   }   

   //SCIPdebugMessage("s2\n");
   // find all bilinear reformulation constraints

   //SCIP_CONSHDLR* conshdlrlinear = SCIPfindConshdlr(scip, "linear");
   //SCIPdebugMessage("s2 %d\n", SCIPconshdlrGetNConss(conshdlrlinear));
   for(int c = 0; c < SCIPconshdlrGetNConss(conshdlrlinear); ++c )
   {
      SCIP_CONS* cons;

	  //SCIPprintCons(scip, cons, NULL);

      cons = SCIPconshdlrGetConss(conshdlrlinear)[c];
      assert(cons != NULL);

      int nvars = SCIPgetNVarsLinear(scip, cons);
      SCIP_VAR ** vars = SCIPgetVarsLinear(scip, cons);
      SCIP_Real * vals = SCIPgetValsLinear(scip, cons);
      SCIP_Bool is_mlf = FALSE;

      vector<SCIP_Real> linear_term_coeffs;
      vector<SCIP_VAR *> linear_term_vars;
      vector<SCIP_Real> mlf_term_coeffs;
      vector<vector<SCIP_VAR*>> mlf_term_vars;
      SCIP_Real constant = 0;
      for(int i = 0; i < nvars; i++){
         if(multilinear_map.find(vars[i]) != multilinear_map.end()){
            mlf_term_vars.push_back(multilinear_map[vars[i]]);
            mlf_term_coeffs.push_back(vals[i]);
            is_mlf = TRUE;
         }
         else{
            linear_term_vars.push_back(vars[i]);
            linear_term_coeffs.push_back(vals[i]);    
            //SCIPdebugMessage("%f\n", vals[i]);        
         }
      }         //SCIPdebugMsg(scip,"\n");
      //SCIPprintCons(scip, cons, NULL);
      //SCIPdebugMsg(scip,"\n");
      //SCIPdebugMessage("%d\n", quad_term_vars.size());
      if( mlf_term_vars.size() != 0){         
         //SCIPdebugMsg(scip,"\n");
         //SCIPprintCons(scip, cons, NULL);
         //SCIPdebugMsg(scip,"\n");
         SCIP_Bool success = FALSE;
         SCIP_Real lhs = SCIPconsGetLhs(scip, cons, &success); 
         SCIP_Real rhs = SCIPconsGetRhs(scip, cons, &success); 
         // lhs <= f(x,t) <= rhs
		 printf("lhs\n");
         if(!SCIPisInfinity(scip, -lhs)){
            createMLFCons(mlf, lhs, linear_term_coeffs, linear_term_vars, mlf_term_coeffs, mlf_term_vars, constant, multilinear_map, TRUE);
         }
		 printf("rhs\n");
         if(!SCIPisInfinity(scip, rhs)){
            createMLFCons(mlf, rhs, linear_term_coeffs, linear_term_vars, mlf_term_coeffs, mlf_term_vars, constant, multilinear_map, FALSE);
         }
      }

   }  

   //SCIPABORT();
   return SCIP_OKAY;   
}