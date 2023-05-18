/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__NMDA
#define _nrn_initial _nrn_initial__NMDA
#define nrn_cur _nrn_cur__NMDA
#define _nrn_current _nrn_current__NMDA
#define nrn_jacob _nrn_jacob__NMDA
#define nrn_state _nrn_state__NMDA
#define _net_receive _net_receive__NMDA 
#define state state__NMDA 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define tau1 _p[0]
#define tau1_columnindex 0
#define tau2 _p[1]
#define tau2_columnindex 1
#define e _p[2]
#define e_columnindex 2
#define weight _p[3]
#define weight_columnindex 3
#define i _p[4]
#define i_columnindex 4
#define g _p[5]
#define g_columnindex 5
#define factor _p[6]
#define factor_columnindex 6
#define cond_scale _p[7]
#define cond_scale_columnindex 7
#define tp _p[8]
#define tp_columnindex 8
#define A _p[9]
#define A_columnindex 9
#define B _p[10]
#define B_columnindex 10
#define rescale_epsps _p[11]
#define rescale_epsps_columnindex 11
#define DA _p[12]
#define DA_columnindex 12
#define DB _p[13]
#define DB_columnindex 13
#define v _p[14]
#define v_columnindex 14
#define _g _p[15]
#define _g_columnindex 15
#define _tsav _p[16]
#define _tsav_columnindex 16
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_mgblock(void*);
 static double _hoc_table_t_table(void*);
 static double _hoc_t_table(void*);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "mgblock", _hoc_mgblock,
 "table_t_table", _hoc_table_t_table,
 "t_table", _hoc_t_table,
 0, 0
};
#define _f_mgblock _f_mgblock_NMDA
#define mgblock mgblock_NMDA
#define table_t_table table_t_table_NMDA
#define t_table t_table_NMDA
 extern double _f_mgblock( _threadargsprotocomma_ double );
 extern double mgblock( _threadargsprotocomma_ double );
 extern double table_t_table( );
 extern double t_table( _threadargsprotocomma_ double );
 
static void _check_mgblock(double*, Datum*, Datum*, NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, int _type) {
   _check_mgblock(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
#define mg mg_NMDA
 double mg = 1;
#define usetable usetable_NMDA
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_NMDA", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "mg_NMDA", "mM",
 "tau1", "ms",
 "tau2", "ms",
 "e", "mV",
 "weight", "uS",
 "A", "uS",
 "B", "uS",
 "i", "nA",
 "g", "uS",
 "factor", "1",
 "cond_scale", "1",
 "tp", "ms",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "mg_NMDA", &mg_NMDA,
 "usetable_NMDA", &usetable_NMDA,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"NMDA",
 "tau1",
 "tau2",
 "e",
 "weight",
 0,
 "i",
 "g",
 "factor",
 "cond_scale",
 "tp",
 0,
 "A",
 "B",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 17, _prop);
 	/*initialize range parameters*/
 	tau1 = 0.5;
 	tau2 = 44;
 	e = 0;
 	weight = 0.0005;
  }
 	_prop->param = _p;
 	_prop->param_size = 17;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _exp2synNMDA_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 17, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 NMDA exp2synNMDA.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_mgblock;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 
static void* _ptable_t_table = (void*)0;
 static double _n_mgblock(_threadargsprotocomma_ double _lv);
 static int _slist1[2], _dlist1[2];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   DA = - A / tau1 ;
   DB = - B / tau2 ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 DA = DA  / (1. - dt*( ( - 1.0 ) / tau1 )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau2 )) ;
  return 0;
}
 /*END CVODE*/
 static int state (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
    A = A + (1. - exp(dt*(( - 1.0 ) / tau1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - A) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tau2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - B) ;
   }
  return 0;
}
 static double _mfac_mgblock, _tmin_mgblock;
  static void _check_mgblock(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_mg;
  if (!usetable) {return;}
  if (_sav_mg != mg) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_mgblock =  - 140.0 ;
   _tmax =  80.0 ;
   _dx = (_tmax - _tmin_mgblock)/1000.; _mfac_mgblock = 1./_dx;
   for (_i=0, _x=_tmin_mgblock; _i < 1001; _x += _dx, _i++) {
    _t_mgblock[_i] = _f_mgblock(_p, _ppvar, _thread, _nt, _x);
   }
   _sav_mg = mg;
  }
 }

 double mgblock(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _lv) { 
#if 0
_check_mgblock(_p, _ppvar, _thread, _nt);
#endif
 return _n_mgblock(_p, _ppvar, _thread, _nt, _lv);
 }

 static double _n_mgblock(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_mgblock(_p, _ppvar, _thread, _nt, _lv); 
}
 _xi = _mfac_mgblock * (_lv - _tmin_mgblock);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_mgblock[0];
 }
 if (_xi >= 1000.) {
 return _t_mgblock[1000];
 }
 _i = (int) _xi;
 return _t_mgblock[_i] + (_xi - (double)_i)*(_t_mgblock[_i+1] - _t_mgblock[_i]);
 }

 
double _f_mgblock ( _threadargsprotocomma_ double _lv ) {
   double _lmgblock;
 _lmgblock = 1.0 / ( 1.0 + exp ( 0.062 * - _lv ) * ( mg / 3.57 ) ) ;
   
return _lmgblock;
 }
 
static double _hoc_mgblock(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 
#if 1
 _check_mgblock(_p, _ppvar, _thread, _nt);
#endif
 _r =  mgblock ( _p, _ppvar, _thread, _nt, *getarg(1) );
 return(_r);
}
 
double t_table ( _threadargsprotocomma_ double _lg ) {
 double _arg[1];
 _arg[0] = _lg;
 return hoc_func_table(_ptable_t_table, 1, _arg);
 }
/*  }
  */
 
static double _hoc_t_table(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  t_table ( _p, _ppvar, _thread, _nt, *getarg(1) );
 return(_r);
}
 double table_t_table ( ) {
	hoc_spec_table(&_ptable_t_table, 1);
	return 0.;
}
 
static double _hoc_table_t_table(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  table_t_table (  );
 return(_r);
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _thread = (Datum*)0; _nt = (NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   double _lt_tmp ;
 _lt_tmp = t_table ( _threadargscomma_ g / _args[0] ) ;
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A;
    double __primary = (_args[0] * factor * exp ( - _lt_tmp / tau1 )) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau1 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - __primary );
    A += __primary;
  } else {
 A = _args[0] * factor * exp ( - _lt_tmp / tau1 ) ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B;
    double __primary = (_args[0] * factor * exp ( - _lt_tmp / tau2 )) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau2 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - __primary );
    B += __primary;
  } else {
 B = _args[0] * factor * exp ( - _lt_tmp / tau2 ) ;
     }
 } }
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  A = A0;
  B = B0;
 {
   if ( tau1 / tau2 > .9999 ) {
     tau1 = .9999 * tau2 ;
     }
   A = 0.0 ;
   B = 0.0 ;
   tp = ( tau1 * tau2 ) / ( tau2 - tau1 ) * log ( tau2 / tau1 ) ;
   factor = - exp ( - tp / tau1 ) + exp ( - tp / tau2 ) ;
   factor = 1.0 / factor ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];

#if 0
 _check_mgblock(_p, _ppvar, _thread, _nt);
#endif
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g = B - A ;
   i = g * mgblock ( _threadargscomma_ v ) * ( v - e ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {   state(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = A_columnindex;  _dlist1[0] = DA_columnindex;
 _slist1[1] = B_columnindex;  _dlist1[1] = DB_columnindex;
   _t_mgblock = makevector(1001*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "exp2synNMDA.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "\n"
  "Combines NEURON'S Exp2Syn mechanism with 1) a previously reported\n"
  "voltage-dependence of Mg2+ block (Obtained from ModelDB model #150239)\n"
  "and 2) a new method to limit the total conductance to a maximum value to\n"
  "produce a computationally efficient NMDAR conductance.\n"
  "\n"
  "Explanation of the mechanism where multiple events (presynaptic\n"
  "activations) cause this single NMDA synapse to saturate at most to the\n"
  "supplied conductance, weight, (the blocking factor reduces the total\n"
  "conductance) rather than summate:\n"
  "\n"
  "Let's forget about the Mg2+ blocking factor for a moment and just look\n"
  "at the factor g=A-B inherited from the exp2syn synapse. The g trace\n"
  "rises and falls similar to an alpha function, however is formed by the\n"
  "difference of two decaying exponentials A and B.  Each event resets\n"
  "\"g\" to the rising curve before the peak (if it was already there it is\n"
  "unchanged, if it was after the peak it is moved to that earlier point\n"
  "that shares the same g value). Therefor at the time of an event the\n"
  "value of g remains constant while the slope g' will change\n"
  "from falling to rising if it is not already rising.\n"
  "\n"
  "Tom Morse\n"
  "\n"
  "Original comment:\n"
  "\n"
  "Two state kinetic scheme synapse described by rise time tau1,\n"
  "and decay time constant tau2. The normalized peak condunductance is 1.\n"
  "Decay time MUST be greater than rise time.\n"
  "\n"
  "The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is\n"
  " A = a*exp(-t/tau1) and\n"
  " G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))\n"
  "	where tau1 < tau2\n"
  "\n"
  "If tau2-tau1 -> 0 then we have a alphasynapse.\n"
  "and if tau1 -> 0 then we have just single exponential decay.\n"
  "\n"
  "The factor is evaluated in the\n"
  "initial block such that an event of weight 1 generates a\n"
  "peak conductance of 1.\n"
  "\n"
  "Because the solution is a sum of exponentials, the\n"
  "coupled equations can be solved as a pair of independent equations\n"
  "by the more efficient cnexp method.\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS NMDA : Exp2SynNMDA\n"
  "	RANGE tau1, tau2, e, i, weight, A, B\n"
  "	NONSPECIFIC_CURRENT i\n"
  "\n"
  "	RANGE g, tp, factor, cond_scale, offset\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(uS) = (microsiemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	tau1= 0.5 (ms) : 5 (ms) .1 (ms) <1e-9, tau2>\n"
  "	tau2 = 44 (ms) : 100 (ms) <tau1,1e9>\n"
  "	e=0	(mV)\n"
  "	mg=1    (mM)		: external magnesium concentration\n"
  "        weight = 0.5e-3 (uS) : reset by events\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	i (nA)\n"
  "	g (uS)\n"
  "	factor (1)\n"
  "        rescale_epsps (1)\n"
  "        cond_scale (1)\n"
  "        tp (ms)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	A (uS)\n"
  "	B (uS)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	if (tau1/tau2 > .9999) {\n"
  "		tau1 = .9999*tau2\n"
  "	}\n"
  "	A = 0\n"
  "	B = 0\n"
  "	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)\n"
  "	factor = -exp(-tp/tau1) + exp(-tp/tau2)\n"
  "	factor = 1/factor\n"
  "}\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	g = B - A\n"
  "	i = g*mgblock(v)*(v - e)\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "	A' = -A/tau1\n"
  "	B' = -B/tau2\n"
  "}\n"
  "\n"
  "FUNCTION mgblock(v(mV)) {\n"
  "	TABLE \n"
  "	DEPEND mg\n"
  "	FROM -140 TO 80 WITH 1000\n"
  "\n"
  "	: from Jahr & Stevens\n"
  "\n"
  "	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))\n"
  "}\n"
  "\n"
  ": For clarity the mg voltage dependent term is ommitted from the below \n"
  ": discussion however it of course is part of the overall conductance.\n"
  ": t_table returns the time on the rising curve of conductance that matches\n"
  ": the conductance (so that if the conductance is falling after the peak it\n"
  ": can rise to the peak again on the occurence of an event).\n"
  ": To recap:  if the conductance, g, at the time of an event was falling\n"
  ": the conductance is effectively moved to the same value however to\n"
  ": the point on the g=A-B curve thats rising.  If the g=A-B point\n"
  ": is rising when the event occurs then that point is looked up again\n"
  ": in the time table, t_table, so the overall effect is that events\n"
  ": occuring on the rising conductance curve have no effect (the g=A-B\n"
  ": continues to rise to the peak.\n"
  ": The function table values are set from hoc in set_t_table.hoc\n"
  "\n"
  "FUNCTION_TABLE t_table(g) (ms)  \n"
  "\n"
  "NET_RECEIVE(weight (uS)) {\n"
  "        LOCAL t_tmp\n"
  "        t_tmp = t_table(g/weight)\n"
  "        A = weight * factor * exp( -t_tmp/tau1)\n"
  "        B = weight * factor * exp( -t_tmp/tau2)\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  ;
#endif
