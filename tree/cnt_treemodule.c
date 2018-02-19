/* File: cnt_treemodule.c
 * This file is auto-generated with f2py (version:2).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * See http://cens.ioc.ee/projects/f2py2e/
 * Generation date: Mon Nov 27 19:34:30 2017
 * $Revision:$
 * $Date:$
 * Do not edit this file directly unless you know what you are doing!!!
 */

#ifdef __cplusplus
extern "C" {
#endif

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "Python.h"
#include <stdarg.h>
#include "fortranobject.h"
#include <string.h>
#include <math.h>

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *cnt_tree_error;
static PyObject *cnt_tree_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
typedef char * string;

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/
\
#define FAILNULL(p) do {                                            \
    if ((p) == NULL) {                                              \
        PyErr_SetString(PyExc_MemoryError, "NULL pointer found");   \
        goto capi_fail;                                             \
    }                                                               \
} while (0)

#define STRINGMALLOC(str,len)\
  if ((str = (string)malloc(sizeof(char)*(len+1))) == NULL) {\
    PyErr_SetString(PyExc_MemoryError, "out of memory");\
    goto capi_fail;\
  } else {\
    (str)[len] = '\0';\
  }

#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif

#define rank(var) var ## _Rank
#define shape(var,dim) var ## _Dims[dim]
#define old_rank(var) (PyArray_NDIM((PyArrayObject *)(capi_ ## var ## _tmp)))
#define old_shape(var,dim) PyArray_DIM(((PyArrayObject *)(capi_ ## var ## _tmp)),dim)
#define fshape(var,dim) shape(var,rank(var)-dim-1)
#define len(var) shape(var,0)
#define flen(var) fshape(var,0)
#define old_size(var) PyArray_SIZE((PyArrayObject *)(capi_ ## var ## _tmp))
/* #define index(i) capi_i ## i */
#define slen(var) capi_ ## var ## _len
#define size(var, ...) f2py_size((PyArrayObject *)(capi_ ## var ## _tmp), ## __VA_ARGS__, -1)

#define STRINGFREE(str) do {if (!(str == NULL)) free(str);} while (0)

#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
  PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
  fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif

#define CHECKSCALAR(check,tcheck,name,show,var)\
  if (!(check)) {\
    char errstring[256];\
    sprintf(errstring, "%s: "show, "("tcheck") failed for "name, var);\
    PyErr_SetString(cnt_tree_error,errstring);\
    /*goto capi_fail;*/\
  } else 
#define STRINGCOPYN(to,from,buf_size)                           \
    do {                                                        \
        int _m = (buf_size);                                    \
        char *_to = (to);                                       \
        char *_from = (from);                                   \
        FAILNULL(_to); FAILNULL(_from);                         \
        (void)strncpy(_to, _from, sizeof(char)*_m);             \
        _to[_m-1] = '\0';                                      \
        /* Padding with spaces instead of nulls */              \
        for (_m -= 2; _m >= 0 && _to[_m] == '\0'; _m--) {      \
            _to[_m] = ' ';                                      \
        }                                                       \
    } while (0)


/************************ See f2py2e/cfuncs.py: cfuncs ************************/
static int f2py_size(PyArrayObject* var, ...)
{
  npy_int sz = 0;
  npy_int dim;
  npy_int rank;
  va_list argp;
  va_start(argp, var);
  dim = va_arg(argp, npy_int);
  if (dim==-1)
    {
      sz = PyArray_SIZE(var);
    }
  else
    {
      rank = PyArray_NDIM(var);
      if (dim>=1 && dim<=rank)
        sz = PyArray_DIM(var, dim-1);
      else
        fprintf(stderr, "f2py_size: 2nd argument value=%d fails to satisfy 1<=value<=%d. Result will be 0.\n", dim, rank);
    }
  va_end(argp);
  return sz;
}

static int string_from_pyobj(string *str,int *len,const string inistr,PyObject *obj,const char *errmess) {
  PyArrayObject *arr = NULL;
  PyObject *tmp = NULL;
#ifdef DEBUGCFUNCS
fprintf(stderr,"string_from_pyobj(str='%s',len=%d,inistr='%s',obj=%p)\n",(char*)str,*len,(char *)inistr,obj);
#endif
  if (obj == Py_None) {
    if (*len == -1)
      *len = strlen(inistr); /* Will this cause problems? */
    STRINGMALLOC(*str,*len);
    STRINGCOPYN(*str,inistr,*len+1);
    return 1;
  }
  if (PyArray_Check(obj)) {
    if ((arr = (PyArrayObject *)obj) == NULL)
      goto capi_fail;
    if (!ISCONTIGUOUS(arr)) {
      PyErr_SetString(PyExc_ValueError,"array object is non-contiguous.");
      goto capi_fail;
    }
    if (*len == -1)
      *len = (PyArray_ITEMSIZE(arr))*PyArray_SIZE(arr);
    STRINGMALLOC(*str,*len);
    STRINGCOPYN(*str,PyArray_DATA(arr),*len+1);
    return 1;
  }
  if (PyString_Check(obj)) {
    tmp = obj;
    Py_INCREF(tmp);
  }
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(obj)) {
    tmp = PyUnicode_AsASCIIString(obj);
  }
  else {
    PyObject *tmp2;
    tmp2 = PyObject_Str(obj);
    if (tmp2) {
      tmp = PyUnicode_AsASCIIString(tmp2);
      Py_DECREF(tmp2);
    }
    else {
      tmp = NULL;
    }
  }
#else
  else {
    tmp = PyObject_Str(obj);
  }
#endif
  if (tmp == NULL) goto capi_fail;
  if (*len == -1)
    *len = PyString_GET_SIZE(tmp);
  STRINGMALLOC(*str,*len);
  STRINGCOPYN(*str,PyString_AS_STRING(tmp),*len+1);
  Py_DECREF(tmp);
  return 1;
capi_fail:
  Py_XDECREF(tmp);
  {
    PyObject* err = PyErr_Occurred();
    if (err==NULL) err = cnt_tree_error;
    PyErr_SetString(err,errmess);
  }
  return 0;
}

static int int_from_pyobj(int* v,PyObject *obj,const char *errmess) {
  PyObject* tmp = NULL;
  if (PyInt_Check(obj)) {
    *v = (int)PyInt_AS_LONG(obj);
    return 1;
  }
  tmp = PyNumber_Int(obj);
  if (tmp) {
    *v = PyInt_AS_LONG(tmp);
    Py_DECREF(tmp);
    return 1;
  }
  if (PyComplex_Check(obj))
    tmp = PyObject_GetAttrString(obj,"real");
  else if (PyString_Check(obj) || PyUnicode_Check(obj))
    /*pass*/;
  else if (PySequence_Check(obj))
    tmp = PySequence_GetItem(obj,0);
  if (tmp) {
    PyErr_Clear();
    if (int_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
    Py_DECREF(tmp);
  }
  {
    PyObject* err = PyErr_Occurred();
    if (err==NULL) err = cnt_tree_error;
    PyErr_SetString(err,errmess);
  }
  return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
extern void F_FUNC_US(count_tree,COUNT_TREE)(string,int*,int*,int*,int*,int*,size_t);
extern void F_FUNC_US(load_tree,LOAD_TREE)(string,int*,int*,int*,float*,int*,float*,float*,float*,float*,int*,int*,int*,int*,int*,size_t);
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/********************************* count_tree *********************************/
static char doc_f2py_rout_cnt_tree_count_tree[] = "\
count_tree(fn,n_halos_all,flist_index,slist_index,nsteps,big_run)\n\nWrapper for ``count_tree``.\
\n\nParameters\n----------\n"
"fn : input string(len=256)\n"
"n_halos_all : input int\n"
"flist_index : input int\n"
"slist_index : input int\n"
"nsteps : input int\n"
"big_run : input int";
/* extern void F_FUNC_US(count_tree,COUNT_TREE)(string,int*,int*,int*,int*,int*,size_t); */
static PyObject *f2py_rout_cnt_tree_count_tree(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(string,int*,int*,int*,int*,int*,size_t)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  string fn = NULL;
  int slen(fn);
  PyObject *fn_capi = Py_None;
  int n_halos_all = 0;
  PyObject *n_halos_all_capi = Py_None;
  int flist_index = 0;
  PyObject *flist_index_capi = Py_None;
  int slist_index = 0;
  PyObject *slist_index_capi = Py_None;
  int nsteps = 0;
  PyObject *nsteps_capi = Py_None;
  int big_run = 0;
  PyObject *big_run_capi = Py_None;
  static char *capi_kwlist[] = {"fn","n_halos_all","flist_index","slist_index","nsteps","big_run",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOOO:cnt_tree.count_tree",\
    capi_kwlist,&fn_capi,&n_halos_all_capi,&flist_index_capi,&slist_index_capi,&nsteps_capi,&big_run_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable fn */
  slen(fn) = 256;
  f2py_success = string_from_pyobj(&fn,&slen(fn),"",fn_capi,"string_from_pyobj failed in converting 1st argument `fn' of cnt_tree.count_tree to C string");
  if (f2py_success) {
  /* Processing variable flist_index */
    f2py_success = int_from_pyobj(&flist_index,flist_index_capi,"cnt_tree.count_tree() 3rd argument (flist_index) can't be converted to int");
  if (f2py_success) {
  /* Processing variable slist_index */
    f2py_success = int_from_pyobj(&slist_index,slist_index_capi,"cnt_tree.count_tree() 4th argument (slist_index) can't be converted to int");
  if (f2py_success) {
  /* Processing variable n_halos_all */
    f2py_success = int_from_pyobj(&n_halos_all,n_halos_all_capi,"cnt_tree.count_tree() 2nd argument (n_halos_all) can't be converted to int");
  if (f2py_success) {
  /* Processing variable nsteps */
    f2py_success = int_from_pyobj(&nsteps,nsteps_capi,"cnt_tree.count_tree() 5th argument (nsteps) can't be converted to int");
  if (f2py_success) {
  /* Processing variable big_run */
    big_run = (int)PyObject_IsTrue(big_run_capi);
    f2py_success = 1;
  if (f2py_success) {
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
        (*f2py_func)(fn,&n_halos_all,&flist_index,&slist_index,&nsteps,&big_run,slen(fn));
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("");
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  } /*if (f2py_success) of big_run*/
  /* End of cleaning variable big_run */
  } /*if (f2py_success) of nsteps*/
  /* End of cleaning variable nsteps */
  } /*if (f2py_success) of n_halos_all*/
  /* End of cleaning variable n_halos_all */
  } /*if (f2py_success) of slist_index*/
  /* End of cleaning variable slist_index */
  } /*if (f2py_success) of flist_index*/
  /* End of cleaning variable flist_index */
    STRINGFREE(fn);
  }  /*if (f2py_success) of fn*/
  /* End of cleaning variable fn */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/***************************** end of count_tree *****************************/

/********************************* load_tree *********************************/
static char doc_f2py_rout_cnt_tree_load_tree[] = "\
load_tree(fn,fatherID,fatherIDx,sonID,fatherMass,i_arr,f_arr,aexp_arr,omega_t_arr,age_univ_arr,big_run,[n_halos_all,n_all_fathers,n_all_sons,nsteps])\n\nWrapper for ``load_tree``.\
\n\nParameters\n----------\n"
"fn : input string(len=256)\n"
"fatherID : input rank-1 array('i') with bounds (n_all_fathers + 1)\n"
"fatherIDx : input rank-1 array('i') with bounds (n_all_fathers + 1)\n"
"sonID : input rank-1 array('i') with bounds (n_all_sons)\n"
"fatherMass : input rank-1 array('f') with bounds (n_all_fathers)\n"
"i_arr : input rank-2 array('i') with bounds (n_halos_all,15)\n"
"f_arr : input rank-2 array('f') with bounds (n_halos_all,25)\n"
"aexp_arr : input rank-1 array('f') with bounds (nsteps)\n"
"omega_t_arr : input rank-1 array('f') with bounds (nsteps)\n"
"age_univ_arr : input rank-1 array('f') with bounds (nsteps)\n"
"big_run : input int\n"
"\nOther Parameters\n----------------\n"
"n_halos_all : input int, optional\n    Default: shape(i_arr,0)\n"
"n_all_fathers : input int, optional\n    Default: (len(fatherID)-1)\n"
"n_all_sons : input int, optional\n    Default: len(sonID)\n"
"nsteps : input int, optional\n    Default: len(aexp_arr)";
/* extern void F_FUNC_US(load_tree,LOAD_TREE)(string,int*,int*,int*,float*,int*,float*,float*,float*,float*,int*,int*,int*,int*,int*,size_t); */
static PyObject *f2py_rout_cnt_tree_load_tree(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(string,int*,int*,int*,float*,int*,float*,float*,float*,float*,int*,int*,int*,int*,int*,size_t)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  string fn = NULL;
  int slen(fn);
  PyObject *fn_capi = Py_None;
  int *fatherID = NULL;
  npy_intp fatherID_Dims[1] = {-1};
  const int fatherID_Rank = 1;
  PyArrayObject *capi_fatherID_tmp = NULL;
  int capi_fatherID_intent = 0;
  PyObject *fatherID_capi = Py_None;
  int *fatherIDx = NULL;
  npy_intp fatherIDx_Dims[1] = {-1};
  const int fatherIDx_Rank = 1;
  PyArrayObject *capi_fatherIDx_tmp = NULL;
  int capi_fatherIDx_intent = 0;
  PyObject *fatherIDx_capi = Py_None;
  int *sonID = NULL;
  npy_intp sonID_Dims[1] = {-1};
  const int sonID_Rank = 1;
  PyArrayObject *capi_sonID_tmp = NULL;
  int capi_sonID_intent = 0;
  PyObject *sonID_capi = Py_None;
  float *fatherMass = NULL;
  npy_intp fatherMass_Dims[1] = {-1};
  const int fatherMass_Rank = 1;
  PyArrayObject *capi_fatherMass_tmp = NULL;
  int capi_fatherMass_intent = 0;
  PyObject *fatherMass_capi = Py_None;
  int *i_arr = NULL;
  npy_intp i_arr_Dims[2] = {-1, -1};
  const int i_arr_Rank = 2;
  PyArrayObject *capi_i_arr_tmp = NULL;
  int capi_i_arr_intent = 0;
  PyObject *i_arr_capi = Py_None;
  float *f_arr = NULL;
  npy_intp f_arr_Dims[2] = {-1, -1};
  const int f_arr_Rank = 2;
  PyArrayObject *capi_f_arr_tmp = NULL;
  int capi_f_arr_intent = 0;
  PyObject *f_arr_capi = Py_None;
  float *aexp_arr = NULL;
  npy_intp aexp_arr_Dims[1] = {-1};
  const int aexp_arr_Rank = 1;
  PyArrayObject *capi_aexp_arr_tmp = NULL;
  int capi_aexp_arr_intent = 0;
  PyObject *aexp_arr_capi = Py_None;
  float *omega_t_arr = NULL;
  npy_intp omega_t_arr_Dims[1] = {-1};
  const int omega_t_arr_Rank = 1;
  PyArrayObject *capi_omega_t_arr_tmp = NULL;
  int capi_omega_t_arr_intent = 0;
  PyObject *omega_t_arr_capi = Py_None;
  float *age_univ_arr = NULL;
  npy_intp age_univ_arr_Dims[1] = {-1};
  const int age_univ_arr_Rank = 1;
  PyArrayObject *capi_age_univ_arr_tmp = NULL;
  int capi_age_univ_arr_intent = 0;
  PyObject *age_univ_arr_capi = Py_None;
  int n_halos_all = 0;
  PyObject *n_halos_all_capi = Py_None;
  int n_all_fathers = 0;
  PyObject *n_all_fathers_capi = Py_None;
  int n_all_sons = 0;
  PyObject *n_all_sons_capi = Py_None;
  int big_run = 0;
  PyObject *big_run_capi = Py_None;
  int nsteps = 0;
  PyObject *nsteps_capi = Py_None;
  static char *capi_kwlist[] = {"fn","fatherID","fatherIDx","sonID","fatherMass","i_arr","f_arr","aexp_arr","omega_t_arr","age_univ_arr","big_run","n_halos_all","n_all_fathers","n_all_sons","nsteps",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOOOOOOOO|OOOO:cnt_tree.load_tree",\
    capi_kwlist,&fn_capi,&fatherID_capi,&fatherIDx_capi,&sonID_capi,&fatherMass_capi,&i_arr_capi,&f_arr_capi,&aexp_arr_capi,&omega_t_arr_capi,&age_univ_arr_capi,&big_run_capi,&n_halos_all_capi,&n_all_fathers_capi,&n_all_sons_capi,&nsteps_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable fn */
  slen(fn) = 256;
  f2py_success = string_from_pyobj(&fn,&slen(fn),"",fn_capi,"string_from_pyobj failed in converting 1st argument `fn' of cnt_tree.load_tree to C string");
  if (f2py_success) {
  /* Processing variable big_run */
    big_run = (int)PyObject_IsTrue(big_run_capi);
    f2py_success = 1;
  if (f2py_success) {
  /* Processing variable fatherID */
  ;
  capi_fatherID_intent |= F2PY_INTENT_IN;
  capi_fatherID_tmp = array_from_pyobj(NPY_INT,fatherID_Dims,fatherID_Rank,capi_fatherID_intent,fatherID_capi);
  if (capi_fatherID_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(cnt_tree_error,"failed in converting 2nd argument `fatherID' of cnt_tree.load_tree to C/Fortran array" );
  } else {
    fatherID = (int *)(PyArray_DATA(capi_fatherID_tmp));

  /* Processing variable sonID */
  ;
  capi_sonID_intent |= F2PY_INTENT_IN;
  capi_sonID_tmp = array_from_pyobj(NPY_INT,sonID_Dims,sonID_Rank,capi_sonID_intent,sonID_capi);
  if (capi_sonID_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(cnt_tree_error,"failed in converting 4th argument `sonID' of cnt_tree.load_tree to C/Fortran array" );
  } else {
    sonID = (int *)(PyArray_DATA(capi_sonID_tmp));

  /* Processing variable i_arr */
  i_arr_Dims[1]=15;
  capi_i_arr_intent |= F2PY_INTENT_IN;
  capi_i_arr_tmp = array_from_pyobj(NPY_INT,i_arr_Dims,i_arr_Rank,capi_i_arr_intent,i_arr_capi);
  if (capi_i_arr_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(cnt_tree_error,"failed in converting 6th argument `i_arr' of cnt_tree.load_tree to C/Fortran array" );
  } else {
    i_arr = (int *)(PyArray_DATA(capi_i_arr_tmp));

  /* Processing variable aexp_arr */
  ;
  capi_aexp_arr_intent |= F2PY_INTENT_IN;
  capi_aexp_arr_tmp = array_from_pyobj(NPY_FLOAT,aexp_arr_Dims,aexp_arr_Rank,capi_aexp_arr_intent,aexp_arr_capi);
  if (capi_aexp_arr_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(cnt_tree_error,"failed in converting 8th argument `aexp_arr' of cnt_tree.load_tree to C/Fortran array" );
  } else {
    aexp_arr = (float *)(PyArray_DATA(capi_aexp_arr_tmp));

  /* Processing variable n_all_sons */
  if (n_all_sons_capi == Py_None) n_all_sons = len(sonID); else
    f2py_success = int_from_pyobj(&n_all_sons,n_all_sons_capi,"cnt_tree.load_tree() 3rd keyword (n_all_sons) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(len(sonID)>=n_all_sons,"len(sonID)>=n_all_sons","3rd keyword n_all_sons","load_tree:n_all_sons=%d",n_all_sons) {
  /* Processing variable n_all_fathers */
  if (n_all_fathers_capi == Py_None) n_all_fathers = (len(fatherID)-1); else
    f2py_success = int_from_pyobj(&n_all_fathers,n_all_fathers_capi,"cnt_tree.load_tree() 2nd keyword (n_all_fathers) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR((len(fatherID)-1)>=n_all_fathers,"(len(fatherID)-1)>=n_all_fathers","2nd keyword n_all_fathers","load_tree:n_all_fathers=%d",n_all_fathers) {
  /* Processing variable n_halos_all */
  if (n_halos_all_capi == Py_None) n_halos_all = shape(i_arr,0); else
    f2py_success = int_from_pyobj(&n_halos_all,n_halos_all_capi,"cnt_tree.load_tree() 1st keyword (n_halos_all) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(shape(i_arr,0)==n_halos_all,"shape(i_arr,0)==n_halos_all","1st keyword n_halos_all","load_tree:n_halos_all=%d",n_halos_all) {
  /* Processing variable nsteps */
  if (nsteps_capi == Py_None) nsteps = len(aexp_arr); else
    f2py_success = int_from_pyobj(&nsteps,nsteps_capi,"cnt_tree.load_tree() 4th keyword (nsteps) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(len(aexp_arr)>=nsteps,"len(aexp_arr)>=nsteps","4th keyword nsteps","load_tree:nsteps=%d",nsteps) {
  /* Processing variable fatherIDx */
  fatherIDx_Dims[0]=n_all_fathers + 1;
  capi_fatherIDx_intent |= F2PY_INTENT_IN;
  capi_fatherIDx_tmp = array_from_pyobj(NPY_INT,fatherIDx_Dims,fatherIDx_Rank,capi_fatherIDx_intent,fatherIDx_capi);
  if (capi_fatherIDx_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(cnt_tree_error,"failed in converting 3rd argument `fatherIDx' of cnt_tree.load_tree to C/Fortran array" );
  } else {
    fatherIDx = (int *)(PyArray_DATA(capi_fatherIDx_tmp));

  /* Processing variable fatherMass */
  fatherMass_Dims[0]=n_all_fathers;
  capi_fatherMass_intent |= F2PY_INTENT_IN;
  capi_fatherMass_tmp = array_from_pyobj(NPY_FLOAT,fatherMass_Dims,fatherMass_Rank,capi_fatherMass_intent,fatherMass_capi);
  if (capi_fatherMass_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(cnt_tree_error,"failed in converting 5th argument `fatherMass' of cnt_tree.load_tree to C/Fortran array" );
  } else {
    fatherMass = (float *)(PyArray_DATA(capi_fatherMass_tmp));

  /* Processing variable f_arr */
  f_arr_Dims[0]=n_halos_all,f_arr_Dims[1]=25;
  capi_f_arr_intent |= F2PY_INTENT_IN;
  capi_f_arr_tmp = array_from_pyobj(NPY_FLOAT,f_arr_Dims,f_arr_Rank,capi_f_arr_intent,f_arr_capi);
  if (capi_f_arr_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(cnt_tree_error,"failed in converting 7th argument `f_arr' of cnt_tree.load_tree to C/Fortran array" );
  } else {
    f_arr = (float *)(PyArray_DATA(capi_f_arr_tmp));

  /* Processing variable omega_t_arr */
  omega_t_arr_Dims[0]=nsteps;
  capi_omega_t_arr_intent |= F2PY_INTENT_IN;
  capi_omega_t_arr_tmp = array_from_pyobj(NPY_FLOAT,omega_t_arr_Dims,omega_t_arr_Rank,capi_omega_t_arr_intent,omega_t_arr_capi);
  if (capi_omega_t_arr_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(cnt_tree_error,"failed in converting 9th argument `omega_t_arr' of cnt_tree.load_tree to C/Fortran array" );
  } else {
    omega_t_arr = (float *)(PyArray_DATA(capi_omega_t_arr_tmp));

  /* Processing variable age_univ_arr */
  age_univ_arr_Dims[0]=nsteps;
  capi_age_univ_arr_intent |= F2PY_INTENT_IN;
  capi_age_univ_arr_tmp = array_from_pyobj(NPY_FLOAT,age_univ_arr_Dims,age_univ_arr_Rank,capi_age_univ_arr_intent,age_univ_arr_capi);
  if (capi_age_univ_arr_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(cnt_tree_error,"failed in converting 10th argument `age_univ_arr' of cnt_tree.load_tree to C/Fortran array" );
  } else {
    age_univ_arr = (float *)(PyArray_DATA(capi_age_univ_arr_tmp));

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
        (*f2py_func)(fn,fatherID,fatherIDx,sonID,fatherMass,i_arr,f_arr,aexp_arr,omega_t_arr,age_univ_arr,&n_halos_all,&n_all_fathers,&n_all_sons,&big_run,&nsteps,slen(fn));
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("");
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  if((PyObject *)capi_age_univ_arr_tmp!=age_univ_arr_capi) {
    Py_XDECREF(capi_age_univ_arr_tmp); }
  }  /*if (capi_age_univ_arr_tmp == NULL) ... else of age_univ_arr*/
  /* End of cleaning variable age_univ_arr */
  if((PyObject *)capi_omega_t_arr_tmp!=omega_t_arr_capi) {
    Py_XDECREF(capi_omega_t_arr_tmp); }
  }  /*if (capi_omega_t_arr_tmp == NULL) ... else of omega_t_arr*/
  /* End of cleaning variable omega_t_arr */
  if((PyObject *)capi_f_arr_tmp!=f_arr_capi) {
    Py_XDECREF(capi_f_arr_tmp); }
  }  /*if (capi_f_arr_tmp == NULL) ... else of f_arr*/
  /* End of cleaning variable f_arr */
  if((PyObject *)capi_fatherMass_tmp!=fatherMass_capi) {
    Py_XDECREF(capi_fatherMass_tmp); }
  }  /*if (capi_fatherMass_tmp == NULL) ... else of fatherMass*/
  /* End of cleaning variable fatherMass */
  if((PyObject *)capi_fatherIDx_tmp!=fatherIDx_capi) {
    Py_XDECREF(capi_fatherIDx_tmp); }
  }  /*if (capi_fatherIDx_tmp == NULL) ... else of fatherIDx*/
  /* End of cleaning variable fatherIDx */
  } /*CHECKSCALAR(len(aexp_arr)>=nsteps)*/
  } /*if (f2py_success) of nsteps*/
  /* End of cleaning variable nsteps */
  } /*CHECKSCALAR(shape(i_arr,0)==n_halos_all)*/
  } /*if (f2py_success) of n_halos_all*/
  /* End of cleaning variable n_halos_all */
  } /*CHECKSCALAR((len(fatherID)-1)>=n_all_fathers)*/
  } /*if (f2py_success) of n_all_fathers*/
  /* End of cleaning variable n_all_fathers */
  } /*CHECKSCALAR(len(sonID)>=n_all_sons)*/
  } /*if (f2py_success) of n_all_sons*/
  /* End of cleaning variable n_all_sons */
  if((PyObject *)capi_aexp_arr_tmp!=aexp_arr_capi) {
    Py_XDECREF(capi_aexp_arr_tmp); }
  }  /*if (capi_aexp_arr_tmp == NULL) ... else of aexp_arr*/
  /* End of cleaning variable aexp_arr */
  if((PyObject *)capi_i_arr_tmp!=i_arr_capi) {
    Py_XDECREF(capi_i_arr_tmp); }
  }  /*if (capi_i_arr_tmp == NULL) ... else of i_arr*/
  /* End of cleaning variable i_arr */
  if((PyObject *)capi_sonID_tmp!=sonID_capi) {
    Py_XDECREF(capi_sonID_tmp); }
  }  /*if (capi_sonID_tmp == NULL) ... else of sonID*/
  /* End of cleaning variable sonID */
  if((PyObject *)capi_fatherID_tmp!=fatherID_capi) {
    Py_XDECREF(capi_fatherID_tmp); }
  }  /*if (capi_fatherID_tmp == NULL) ... else of fatherID*/
  /* End of cleaning variable fatherID */
  } /*if (f2py_success) of big_run*/
  /* End of cleaning variable big_run */
    STRINGFREE(fn);
  }  /*if (f2py_success) of fn*/
  /* End of cleaning variable fn */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/****************************** end of load_tree ******************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/
/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {
  {"count_tree",-1,{{-1}},0,(char *)F_FUNC_US(count_tree,COUNT_TREE),(f2py_init_func)f2py_rout_cnt_tree_count_tree,doc_f2py_rout_cnt_tree_count_tree},
  {"load_tree",-1,{{-1}},0,(char *)F_FUNC_US(load_tree,LOAD_TREE),(f2py_init_func)f2py_rout_cnt_tree_load_tree,doc_f2py_rout_cnt_tree_load_tree},

/*eof routine_defs*/
  {NULL}
};

static PyMethodDef f2py_module_methods[] = {

  {NULL,NULL}
};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "cnt_tree",
  NULL,
  -1,
  f2py_module_methods,
  NULL,
  NULL,
  NULL,
  NULL
};
#endif

#if PY_VERSION_HEX >= 0x03000000
#define RETVAL m
PyMODINIT_FUNC PyInit_cnt_tree(void) {
#else
#define RETVAL
PyMODINIT_FUNC initcnt_tree(void) {
#endif
  int i;
  PyObject *m,*d, *s;
#if PY_VERSION_HEX >= 0x03000000
  m = cnt_tree_module = PyModule_Create(&moduledef);
#else
  m = cnt_tree_module = Py_InitModule("cnt_tree", f2py_module_methods);
#endif
  Py_TYPE(&PyFortran_Type) = &PyType_Type;
  import_array();
  if (PyErr_Occurred())
    {PyErr_SetString(PyExc_ImportError, "can't initialize module cnt_tree (failed to import numpy)"); return RETVAL;}
  d = PyModule_GetDict(m);
  s = PyString_FromString("$Revision: $");
  PyDict_SetItemString(d, "__version__", s);
#if PY_VERSION_HEX >= 0x03000000
  s = PyUnicode_FromString(
#else
  s = PyString_FromString(
#endif
    "This module 'cnt_tree' is auto-generated with f2py (version:2).\nFunctions:\n"
"  count_tree(fn,n_halos_all,flist_index,slist_index,nsteps,big_run)\n"
"  load_tree(fn,fatherID,fatherIDx,sonID,fatherMass,i_arr,f_arr,aexp_arr,omega_t_arr,age_univ_arr,big_run,n_halos_all=shape(i_arr,0),n_all_fathers=(len(fatherID)-1),n_all_sons=len(sonID),nsteps=len(aexp_arr))\n"
".");
  PyDict_SetItemString(d, "__doc__", s);
  cnt_tree_error = PyErr_NewException ("cnt_tree.error", NULL, NULL);
  Py_DECREF(s);
  for(i=0;f2py_routine_defs[i].name!=NULL;i++)
    PyDict_SetItemString(d, f2py_routine_defs[i].name,PyFortranObject_NewAsAttr(&f2py_routine_defs[i]));


/*eof initf2pywraphooks*/
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
  if (! PyErr_Occurred())
    on_exit(f2py_report_on_exit,(void*)"cnt_tree");
#endif

  return RETVAL;
}
#ifdef __cplusplus
}
#endif
