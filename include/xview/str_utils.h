/*	@(#)str_utils.h 20.16 91/09/14 SMI */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_string_utils_DEFINED
#define xview_string_utils_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <xview/sun.h> 		/* defines TRUE, FALSE */
#include <xview/xv_c_types.h>

/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures
 ***********************************************************************
 */

/*
 * CharClass and CharAction are also defined in io_stream.h 
 */

#ifndef CHARCLASS
#define CHARCLASS
enum 	CharClass {Break, Sepr, Other};
#endif /* ~CHARCLASS */

#ifndef CHARACTION
#define CHARACTION
struct CharAction {
	Bool	stop;
	Bool	include;
};
#endif /* ~CHARACTION */

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

/*
 * PUBLIC Functions 
 */

#ifndef WHITESPACE
#define WHITESPACE
/*
 * xv_white_space is also defined in io_stream.h 
 * returns sepr for blanks, newlines, and tabs, other for everything else 
 */
EXTERN_FUNCTION (ENUM_TYPE(enum,CharClass) xv_white_space, (int c));
#endif

EXTERN_FUNCTION (int string_find, (char *s, char *target, Bool case_matters));
/* int string_find(s, target, case_matters)
 *	char *s, *target;
 *	Bool case_matters;
 * strfind searches one instance of a string for another.
 * If successful, returns the position in the string where the match began,
 * otherwise -1.
 * If case_matters = FALSE, 'a' will match with 'a' or 'A'.
 */

EXTERN_FUNCTION (Bool string_equal, (char *s1, char *s2, Bool case_matters));
/* Bool string_equal(s1, s2, case_matters)
 *	char *s1, *s2;
 *	Bool case_matters;
 * strequal compares two strings. 
 * If case_matters = FALSE, 'a' will match with 'a' or 'A'.
 * either s1 or s2 can be NULL without harm.
 */

EXTERN_FUNCTION (Bool xv_substrequal, (char *s1, int start1, char *s2, int start2, int n, Bool case_matters));
/* Bool xv_substrequal(s1, start1, s2, start2, n, case_matters)
 *	char *s1, *s2;
 *	int start1, start2, n;
 *	Bool case_matters;
 * xv_substrequal compares two substrings without having to construct them.
 * If case_matters = FALSE, 'a' will match with 'a' or 'A'.
 */

EXTERN_FUNCTION (char *string_get_token, (char *s, int *index, char *dest, ENUM_TYPE(enum,CharClass) (*charproc)(int)));
/*	char *s;
 *	int *index;
 *	char *dest;
 *	enum CharClass (*charproc) (char c);

 * string_get_token is used for tokenizing input, where more degree of
 * flexibility is required than simply delimiting tokens by white spaces
 * characters are divided into three classes, Break, Sepr, and Other.
 * separators (Sepr) serve to delimit a token. Leading separators are skipped.
 * think of separators as white space. Break characters delimit tokens, and
 * are themselves tokens. Thus, if a break character is the first thing seen
 * it is returned as the token. If any non-separator characters have been seen,
 * then they are returned as the token, and the break character will be the
 * returned as the result of the next call to get_token.
 * for example, if charproc returns Sepr for space, and Break for '(' and ')'
 * and Other for all alphabetic characters, then the string "now (is) the"
 * will yield five tokens consisting of "now" "(" "is" ")" and "the"

 * get_token stores the token that it constructs into dest,
 * which is also returned as its value.
 * index marks the current position in the string to "begin reading from"
 * it is updated so that the client program does not have to keep track of
 * how many characters have been read.

 * get_token returns NULL, rather than the empty string, corresponding to
 * the case where the token is empty
  */

EXTERN_FUNCTION (char *string_get_sequence, (char *s, int *index, char *dest, struct CharAction (*charproc)(void)));
/*	char *s;
 *	int *index;
 *	char *dest;
 *	struct CharAction (*charproc) ();

 * string_get_sequence is a more primitive tokenizer than get_token.
 * it takes a procedure which for each character specifies whether the
 * character is to terminate the sequence, and whether or not the
 * character is to be included in the sequence.
 * (If the character terminates the sequence, but is not included, then
 * it will be seen again on the next call.)
 * For example, having seen a \"\, to read to the matching \"\, call 
 * get_sequence with an action procedure that returns {TRUE, TRUE} for \"\
 * and  {FALSE, TRUE} for everything else. (If you want to detect the
 * case where a " is preceded by a \\, simply save the last character
 * and modify the procedure accordingly.

 * Note that gettoken can be defined in terms of get_sequence by
 * having Other characters return {FALSE, TRUE}, and also noticing whether
 * any have been seen yet, having Seprs return
 * {(seen_some_others ? TRUE : FALSE), FALSE}
 * and Break characters return {TRUE, (seen_some_others ? FALSE : TRUE)}

 * returns NULL for the empty sequence
 */

#endif /* ~xview_string_utils_DEFINED */
