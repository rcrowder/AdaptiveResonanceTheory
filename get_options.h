/* $Id: getopt.c,v 1.4 2003/10/20 15:03:42 chrisgreen Exp $ */
/**
 * Copyright (C) 2002 Martin Roesch <roesch@sourcefire.com><br/>
 **/


#ifndef GET_OPTIONS_H
#define GET_OPTIONS_H


/**
 * The GetOpt class is used to wrap an implementation of the unistd.h function
 * getopt() and the global variables associated with it into a namespace. This
 * implementation is portable to programming environments that don't have the
 * unistd.h header available to them (e.g. Visual C++).
 *
 * <p>
 * The getopt() function is a command line parser.  It returns the next option
 * character in argv that matches an option character in opstring.
 *
 * <p>
 * The argv argument points to an array of argc+1 elements containing argc
 * pointers to character strings followed by a null pointer.
 *
 * <p>
 * The opstring argument points to a string of option characters; if an option
 * character is followed by a colon, the option is expected to have an argument
 * that may or may not be separated from it by white space. The external
 * variable optarg is set to point to the start of the option argument on
 * return from getopt().
 *
 * <p>
 * The getopt() function places in optind the argv index of the next argument
 * to be processed.  The system initializes the external variable optind to 1
 * before the first call to getopt().
 *
 * <p>
 * When all options have been processed (that is, up to the first nonoption
 * argument), getopt() returns EOF.  The special option "--" may be used to
 * delimit the end of the options; EOF will be returned, and "--" will be
 * skipped.
 *
 * <p>
 * The getopt() function returns a question mark (?) when it encounters an
 * option character not included in opstring.  This error message can be
 * disabled by setting opterr to zero.  Otherwise, it returns the option
 * character that was detected.
 *
 * <p>
 * If the special option "--" is detected, or all options have been processed,
 * EOF is returned.
 *
 * <p>
 * Options are marked by either a minus sign (-) or a slash (/).
 *
 * <p>
 * No errors are defined.
 */
class GetOpt {

public:
  // Public Static Member Variables ///////////////////////////////////////////
  static char* optarg;
  static int   optind;
  static int   opterr;

  // Public Static Method /////////////////////////////////////////////////////
  static int getopt(int argc, char* argv[], char* opstring);

private:
  // Private Static Member Variable ///////////////////////////////////////////
  static char *pIndexPosition;
};


#endif

