/* $Id: getopt.c,v 1.4 2003/10/20 15:03:42 chrisgreen Exp $ */
/**
 * Copyright (C) 2002 Martin Roesch <roesch@sourcefire.com><br/>
 **/

#include <cstdio>
#include <cstring>
#include "get_options.h"


// Public Static Member Variables /////////////////////////////////////////////

char* GetOpt::optarg = NULL;
int   GetOpt::optind = 1;
int   GetOpt::opterr = 1;


// Private Static Member Variable /////////////////////////////////////////////

char* GetOpt::pIndexPosition = NULL;


// Public Static Method ///////////////////////////////////////////////////////

int GetOpt::getopt(int argc, char *argv[], char *opstring) {
  char *pArgString = NULL;
  char *pOptString;

  if (pIndexPosition != NULL)
    if (*(++pIndexPosition))
      pArgString = pIndexPosition;

  if (pArgString == NULL) {
    if (optind >= argc) {
      pIndexPosition = NULL;
      return EOF;
    }
    pArgString = argv[optind++];
    if (('/' != *pArgString) && ('-' != *pArgString)) {
      --optind;
      optarg = NULL;
      pIndexPosition = NULL;
      return EOF;
    }
    if ((strcmp(pArgString, "-") == 0) || (strcmp(pArgString, "--") == 0)) {
      optarg = NULL;
      pIndexPosition = NULL;
      return EOF;
    }
    pArgString++;
  }

  if (':' == *pArgString) {
    return (opterr ? (int)'?' : (int)':');
  }
  else if ((pOptString = strchr(opstring, *pArgString)) == 0) {
    optarg = NULL;
    pIndexPosition = NULL;
    return (opterr ? (int)'?' : (int)*pArgString);
  }
  else {
    if (':' == ((char)*(pOptString+1))) {
      if ('\0' != ((char)*(pArgString+1)))
        optarg = &pArgString[1];
      else {
        if (optind < argc)
          optarg = argv[optind++];
        else {
          optarg = NULL;
          return (opterr ? (int)'?' : (int)*pArgString);
        }
      }
      pIndexPosition = NULL;
    }
    else {
      optarg = NULL;
      pIndexPosition = pArgString;
    }
    return (int)*pArgString;
  }
}

