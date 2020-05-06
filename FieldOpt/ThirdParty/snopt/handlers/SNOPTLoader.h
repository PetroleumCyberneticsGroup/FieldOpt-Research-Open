/***********************************************************
ADGPRS, version 1.0, Copyright (c) 2010-2017 SUPRI-B
Author(s): Oleg Volkov          (ovolkov@stanford.edu)
***********************************************************/

#ifndef SNOPTLOADER_H_
#define SNOPTLOADER_H_

#include "../include/snopt/snopt.hpp"
//#include "uc_1ThirdParty/snopt/include/snopt/snopt.hpp"
// #include "uc_1ThirdParty/snopt/include/snopt-7.2.12/snopt.hpp"

#define LIBRARYSNOPT "libsnopt."
#define SNOPTLIBNAME LIBRARYSNOPT SHAREDLIBEXT

#ifdef __cplusplus
extern "C" {
#endif

// Tries to load a dynamically linked library with SNOPT.
// Return a failure if the library cannot be loaded or not
// all SNOPT symbols are found.
//
// @param libname The name under which the SNOPT lib can
// be found, or NULL to use a default name.
// @param msgbuf A buffer where we can store a failure
// message. Assumed to be NOT NULL!
// @param msglen Length of the message buffer.
// @return Zero on success, nonzero on failure.

int LSL_loadSNOPTLib(const char* libname,
                     char* msgbuf, int msglen);


// Unloads a loaded SNOPT library.
// @return Zero on success, nonzero on failure.
//
int LSL_unloadSNOPTLib();


// Indicates whether SNOPT lib has been successfully loaded.
// @return Zero if not loaded, nonzero if handle is loaded

int LSL_isSNOPTLoaded();


// Returns name of shared lib that should contain SNOPT
char* LSL_SNOPTLibraryName();

#ifdef __cplusplus
}
#endif

#endif /*SNOPTLOADER_H_*/