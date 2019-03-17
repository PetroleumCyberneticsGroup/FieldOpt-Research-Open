/***********************************************************
ADGPRS, version 1.0, Copyright (c) 2010-2017 SUPRI-B
Author(s): Oleg Volkov          (ovolkov@stanford.edu)
***********************************************************/

#include "LibraryHandler.h"
#include "SNOPTLoader.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

static soHandle_t SNOPT_handle = NULL;

typedef void (*voidfun)(void);

 voidfun LSL_loadSym (soHandle_t h,
                     const char *symName, 
                     char *msgBuf, int msgLen);


typedef void (*snopta_t)(integer*, integer*, integer*,
                         integer*, integer*, doublereal*, 
                         integer*, char*, My_fp, integer*, 
                         integer*, integer*, integer*, 
                         doublereal*, integer*, integer*, 
                         integer*, integer*, doublereal*, 
                         doublereal*, char*, doublereal*, 
                         doublereal*, char*, doublereal*, 
                         integer*, doublereal*, doublereal*, 
                         integer*, doublereal*, integer*,
                         integer*, integer*, integer*, 
                         integer*, integer*, doublereal*, 
                         char*, integer*, integer*, integer*, 
                         doublereal*, integer*, char*, 
                         integer*, integer*, integer*, 
                         doublereal*, integer*, ftnlen, 
                         ftnlen, ftnlen, ftnlen, ftnlen);


typedef void (*sninit_t)(integer*, integer*, char*,
                         integer*, integer*, integer*,
                         doublereal*, integer*, ftnlen);


typedef void (*sngeti_t)(char*, integer*, integer*,
                         char*, integer*, integer*,
                         integer*, doublereal*, integer*,
                         ftnlen, ftnlen);


typedef void (*sngetr_t)(char*, doublereal*, integer*,
                         char*, integer*, integer*,
                         integer*, doublereal*, integer*,
                         ftnlen, ftnlen);


typedef void (*snset_t)(char*, integer*, integer*,
                        integer*, char*, integer*,
                        integer*, integer*,
                        doublereal*, integer*,
                        ftnlen, ftnlen);


typedef void (*sngetc_t)(char*, char*, integer*,
                         char*, integer*, integer*,
                         integer*, doublereal*, integer*,
                         ftnlen, ftnlen, ftnlen);


typedef void (*snseti_t)(char*, integer*, integer*,
                         integer*, integer*, char*,
                         integer*, integer*, integer*,
                         doublereal*, integer*, ftnlen,
                         ftnlen);


typedef void (*snsetr_t)(char*, doublereal*, integer *,
                         integer*, integer*, char*,
                         integer*, integer*, integer*,
                         doublereal*, integer*, ftnlen,
                         ftnlen);


typedef void (*snspec_t)(integer*, integer*, char*,
                         integer*, integer*, integer*,
                         doublereal*, integer*, ftnlen);


typedef void (*snmema_t)(integer*, integer*, integer*, 
                         integer*, integer*, integer*, 
                         integer*, integer*, integer*,
                         integer*, char*, integer*, 
                         integer*, integer*, doublereal*, 
                         integer*, ftnlen);


typedef void (*snjac_t)(integer*, integer*, integer*, 
                        My_fp, integer*, integer*, 
                        integer*, integer*, doublereal*, 
                        integer*, integer*, integer*, 
                        integer*, doublereal*, doublereal*, 
                        doublereal*, integer*, integer*,
                        integer*, char*, integer*, integer*, 
                        integer*, doublereal*, integer*, 
                        char*, integer*, integer*, integer*, 
                        doublereal*, integer*, ftnlen, 
                        ftnlen);


typedef void (*openfile_t)(integer*, char*, 
                           integer*, integer);


typedef void (*closefile_t)(integer*);

static snopta_t    func_snopta = NULL;
static sninit_t    func_sninit = NULL;
static sngeti_t    func_sngeti = NULL;
static sngetr_t    func_sngetr = NULL;
static sngetc_t    func_sngetc = NULL;
static snset_t     func_snset = NULL;
static snseti_t    func_snseti = NULL;
static snsetr_t    func_snsetr = NULL;
static snspec_t    func_snspec = NULL;
static snmema_t    func_snmema = NULL;
static snjac_t     func_snjac = NULL;
static openfile_t  func_openfile = NULL;
static closefile_t func_closefile = NULL;


void snopta_(integer *start, integer *nf, 
             integer *n, integer *nxname, 
             integer *nfname, doublereal *objadd, 
             integer *objrow, char *prob, 
             My_fp usrfun, integer *iafun, 
             integer *javar, integer *lena, 
             integer *nea, doublereal *a, 
             integer *igfun, integer *jgvar, 
             integer *leng, integer *neg, 
             doublereal *xlow, doublereal *xupp, 
             char *xnames, doublereal *flow, 
             doublereal *fupp, char *fnames, 
             doublereal *x, integer *xstate, 
             doublereal *xmul, doublereal *f, 
             integer *fstate, doublereal *fmul, 
             integer *inform__, integer *mincw, 
             integer *miniw, integer *minrw, 
             integer *ns, integer *ninf, 
             doublereal *sinf, char *cu, 
             integer *lencu, integer *iu,
             integer *leniu, doublereal *ru, 
             integer *lenru, char *cw, integer *lencw,
             integer *iw, integer *leniw, 
             doublereal *rw, integer *lenrw,
             ftnlen prob_len, ftnlen xnames_len, 
             ftnlen fnames_len, ftnlen cu_len,
             ftnlen cw_len) {

  if (func_snopta == NULL)
    exit(EXIT_FAILURE);

  func_snopta(start, nf, n, nxname, nfname, objadd, 
              objrow, prob, usrfun, iafun, javar, 
              lena, nea, a, igfun, jgvar, leng, neg, 
              xlow, xupp, xnames, flow, fupp, fnames, 
              x, xstate, xmul, f, fstate, fmul, inform__,
              mincw, miniw, minrw, ns, ninf, sinf, cu, 
              lencu, iu, leniu, ru, lenru, cw, lencw, 
              iw, leniw, rw, lenrw, prob_len, xnames_len, 
              fnames_len, cu_len, cw_len);
}


void sninit_(integer *iPrint, integer *iSumm, 
             char *cw, integer *lencw, integer *iw, 
             integer *leniw, doublereal *rw, 
             integer *lenrw, ftnlen cw_len ) {

  if (func_sninit == NULL)
    exit(EXIT_FAILURE);

  func_sninit(iPrint, iSumm, cw, lencw, 
              iw, leniw, rw, lenrw, cw_len);
}


void sngeti_(char *buffer, integer *ivalue, 
             integer *inform__, char *cw, 
             integer *lencw, integer *iw, 
             integer *leniw, doublereal *rw, 
             integer *lenrw, ftnlen buffer_len, 
             ftnlen cw_len) {

  if (func_sngeti == NULL)
    exit(EXIT_FAILURE);

  func_sngeti(buffer, ivalue, inform__, cw, lencw, iw,
              leniw, rw, lenrw, buffer_len, cw_len);
}


void sngetr_(char *buffer, doublereal *ivalue, 
             integer *inform__, char *cw, 
             integer *lencw, integer *iw, 
             integer *leniw, doublereal *rw, 
             integer *lenrw, ftnlen buffer_len, 
             ftnlen cw_len) {

  if (func_sngetr == NULL)
    exit(EXIT_FAILURE);

  func_sngetr(buffer, ivalue, inform__, cw, lencw, iw,
              leniw, rw, lenrw, buffer_len, cw_len);
}


void snset_(char *buffer, integer *iprint, 
            integer *isumm, integer *inform__, 
            char *cw, integer *lencw, integer *iw, 
            integer *leniw, doublereal *rw, 
            integer *lenrw, ftnlen buffer_len, 
            ftnlen cw_len) {

  if (func_snset == NULL)
    exit(EXIT_FAILURE);

  func_snset(buffer, iprint, isumm, inform__, cw, lencw,
             iw, leniw, rw, lenrw, buffer_len, cw_len);
}


void sngetc_(char *buffer, char *ivalue, 
             integer *inform__, char *cw, 
             integer *lencw, integer *iw, 
             integer *leniw, doublereal *rw, 
             integer *lenrw, ftnlen buffer_len, 
             ftnlen ivalue_len, ftnlen cw_len) {

  if (func_sngetc == NULL)
    exit(EXIT_FAILURE);

  func_sngetc(buffer, ivalue, inform__, cw, lencw, 
              iw, leniw, rw, lenrw, buffer_len, 
              ivalue_len, cw_len);
}


void snseti_(char *buffer, integer *ivalue, 
             integer *iprint, integer *isumm, 
             integer *inform__, char *cw, 
             integer *lencw, integer *iw, 
             integer *leniw, doublereal *rw, 
             integer *lenrw, ftnlen buffer_len,
            ftnlen cw_len) {

  if (func_snseti == NULL)
    exit(EXIT_FAILURE);

  func_snseti(buffer, ivalue, iprint, isumm, 
              inform__, cw, lencw, iw, leniw, 
              rw, lenrw, buffer_len, cw_len);
}


void snsetr_(char *buffer, doublereal *rvalue, 
             integer * iprint, integer *isumm, 
             integer *inform__, char *cw, 
             integer *lencw, integer *iw, 
             integer *leniw, doublereal *rw, 
             integer *lenrw, ftnlen buffer_len,
            ftnlen cw_len) {

  if (func_snsetr == NULL)
    exit(EXIT_FAILURE);

  func_snsetr(buffer, rvalue, iprint, isumm, 
              inform__, cw, lencw, iw, leniw, 
              rw, lenrw, buffer_len, cw_len);
}


void snspec_(integer *ispecs, integer *inform__, 
             char *cw, integer *lencw, integer *iw, 
             integer *leniw, doublereal *rw, 
             integer *lenrw, ftnlen cw_len) {

  if (func_snspec == NULL)
    exit(EXIT_FAILURE);

  func_snspec(ispecs, inform__, cw, lencw, 
              iw, leniw, rw, lenrw, cw_len);
}


void snmema_(integer *iexit, integer *nf, integer *n, 
             integer *nxname, integer *nfname, 
             integer *nea, integer *neg, integer *mincw, 
             integer *miniw, integer *minrw, char *cw, 
             integer *lencw, integer *iw, integer *leniw, 
             doublereal *rw, integer *lenrw, ftnlen cw_len)
{
  if (func_snmema == NULL)
    exit(EXIT_FAILURE);

  func_snmema(iexit, nf, n, nxname, nfname, nea, neg,
              mincw, miniw, minrw, cw, lencw, iw,
              leniw, rw, lenrw, cw_len);
}


void snjac_(integer *inform__, integer *nf, integer *n, 
            My_fp userfg, integer *iafun, integer *javar, 
            integer *lena, integer *nea, doublereal *a, 
            integer *igfun, integer *jgvar, integer *leng, 
            integer *neg, doublereal *x, doublereal *xlow, 
            doublereal *xupp, integer *mincw, integer *miniw,
            integer *minrw, char *cu, integer *lencu,
            integer *iu, integer *leniu, doublereal *ru,
            integer *lenru, char *cw, integer *lencw, 
            integer *iw, integer *leniw, doublereal *rw, 
            integer *lenrw, ftnlen cu_len, ftnlen cw_len ) {

  if (func_snjac == NULL)
    exit(EXIT_FAILURE);

  func_snjac(inform__, nf, n, userfg, iafun, javar, lena,
             nea, a, igfun, jgvar, leng, neg,
             x, xlow, xupp, mincw, miniw,
             minrw, cu, lencu, iu, leniu, ru,
             lenru, cw, lencw, iw, leniw, rw, lenrw,
             cu_len, cw_len);
}


void openfile_(integer* num, char* name, 
               integer* mode, integer prnt_len) {

  if (func_openfile == NULL)
    exit(EXIT_FAILURE);

  func_openfile(num, name, mode, prnt_len);
}


void closefile_(integer* num) {
  if (func_closefile == NULL)
    exit(EXIT_FAILURE);

  func_closefile(num);
}


int LSL_loadSNOPTLib(const char* libname, 
                     char* msgbuf, 
                     int msglen) {

  if (libname) {
    SNOPT_handle = LSL_loadLib(libname, msgbuf, msglen);
  } else  { /* try a default library name */
    SNOPT_handle = LSL_loadLib(SNOPTLIBNAME, msgbuf, msglen);
  }

  if (SNOPT_handle == NULL) {
    return EXIT_FAILURE;
  }

  func_snopta = (snopta_t)LSL_loadSym(SNOPT_handle, 
                                      "snopta_", 
                                      msgbuf, msglen);
  if (func_snopta == NULL) {
    return EXIT_FAILURE;
  }


  func_sninit = (sninit_t)LSL_loadSym(SNOPT_handle, 
                                      "sninit_", 
                                      msgbuf, msglen);

  if (func_sninit == NULL) {
    return EXIT_FAILURE;
  }


  func_sngeti = (sngeti_t)LSL_loadSym(SNOPT_handle, 
                                      "sngeti_", 
                                      msgbuf, msglen);

  if (func_sngeti == NULL) {
    return EXIT_FAILURE;
  }


  func_sngetr = (sngetr_t)LSL_loadSym(SNOPT_handle, 
                                      "sngetr_", 
                                      msgbuf, msglen);

  if (func_sngetr == NULL) {
    return EXIT_FAILURE;
  }


  func_sngetc = (sngetc_t)LSL_loadSym(SNOPT_handle, 
                                      "sngetc_", 
                                      msgbuf, msglen);

  if (func_sngetc == NULL) {
    return EXIT_FAILURE;
  }


  func_snset = (snset_t)LSL_loadSym(SNOPT_handle, 
                                    "snset_", 
                                    msgbuf, msglen);

  if (func_snset == NULL) {
    return EXIT_FAILURE;
  }


  func_snseti = (snseti_t)LSL_loadSym(SNOPT_handle, 
                                      "snseti_", 
                                      msgbuf, msglen);

  if (func_snseti == NULL) {
    return EXIT_FAILURE;
  }


  func_snsetr = (snsetr_t)LSL_loadSym(SNOPT_handle, 
                                      "snsetr_", 
                                      msgbuf, msglen);

  if (func_snsetr == NULL) {
    return EXIT_FAILURE;
  }


  func_snspec = (snspec_t)LSL_loadSym(SNOPT_handle, 
                                      "snspec_", 
                                      msgbuf, msglen);

  if (func_snspec == NULL) {
    return EXIT_FAILURE;
  }


  func_snmema = (snmema_t)LSL_loadSym(SNOPT_handle, 
                                      "snmema_", 
                                      msgbuf, msglen);

  if (func_snmema == NULL) {
    return EXIT_FAILURE;
  }


  func_snjac = (snjac_t)LSL_loadSym(SNOPT_handle, 
                                    "snjac_", 
                                    msgbuf, msglen);

  if (func_snjac == NULL) {
    return EXIT_FAILURE;
  }

// -------------------------------------------------------
  func_openfile = (openfile_t)LSL_loadSym(SNOPT_handle, 
                                          "openfile_", 
                                          msgbuf, msglen);

  if (func_openfile == NULL) {
    return EXIT_FAILURE;
  }

// -------------------------------------------------------
  func_closefile = (closefile_t)LSL_loadSym(SNOPT_handle, 
                                            "closefile_", 
                                            msgbuf, msglen);

  if (func_closefile == NULL) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


int LSL_unloadSNOPTLib() {

  int rc;

  if (SNOPT_handle == NULL) {
    return EXIT_SUCCESS;
  }


  rc = LSL_unloadLib(SNOPT_handle);
  SNOPT_handle = NULL;

  func_snopta = NULL;
  func_sninit = NULL;
  func_sngeti = NULL;
  func_sngetr = NULL;
  func_sngetc = NULL;
  func_snset = NULL;
  func_snseti = NULL;
  func_snsetr = NULL;
  func_snspec = NULL;
  func_snmema = NULL;
  func_snjac = NULL;
  func_openfile = NULL;
  func_closefile = NULL;

  return rc;
}


int LSL_isSNOPTLoaded() {
  return SNOPT_handle != NULL;
}


char* LSL_SNOPTLibraryName() {
  static char name[] = SNOPTLIBNAME;
  return name;
}