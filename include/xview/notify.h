
/*      @(#)notify.h 20.30 91/09/14 SMI      */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef	xview_notify_DEFINED
#define	xview_notify_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <stdio.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <xview/base.h>
#ifdef SYSV_WAIT
#include <sys/rusage.h>
#endif 
#ifdef SYSV_UCONTEXT
#include <sys/ucontext.h>
#endif 

/* This is part of the hack to make fcntl work for SVR4. The problem
* is in short: because we are includeing the socket libarary, which
* has an unoverridable fcntl, this does not work anymore for SVR4.
* We have to redefine fcntl to xv_fcntl, and we have to replace the
* system call with a call to ther real fcntl.
* The argument to perror has not been changed.
* This hack works, as long as all the calls to fcntl have the same
* number of arguments. They do right now. [vmh - 7/31/90]
*/
#ifdef XV_USE_XVFCNTL
#define fcntl(a,b,c) xv_fcntl(a,b,c)
#endif 
/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

#define	NOTIFY_FUNC_NULL	((Notify_func)0)

/*      Macros to examine wait3/waitpid results provided for BSD rather
 *      than adding a lot of #ifdef'ed code in the .c's.  Note that AT&T-
 *      style macros expect the status value (not a pointer to it).
 *      Also provide a dummy rusage structure in the non-BSD case.
 */
#if  !(defined WTERMSIG) && !(defined SYSV_WAIT)
#define WTERMSIG(status)        ((status).w_termsig)
#endif
#if  !(defined WSTOPSIG) && !(defined SYSV_WAIT)
#define WSTOPSIG(status)        ((status).w_stopsig)
#endif
#if  !(defined WEXITSTATUS) && !(defined SYSV_WAIT)
#define WEXITSTATUS(status)     ((status).w_retcode)
#endif
#if  !(defined WCOREDUMP) && !(defined SYSV_WAIT)
#define WCOREDUMP(status)       ((status).w_coredump)
#endif

#ifdef POLL
#define notify_set_fd_func(nclient, func, poll_fd) \
                           ndet_set_fd_func(nclient, func, poll_fd, NTFY_FD)
#endif /* POLL */


/*
 * PRIVATE #defines 
 */

#define NOTIFY_CLIENT_NULL	((Notify_client)0)
#define	NOTIFY_COPY_NULL	((Notify_copy)0)
#define	NOTIFY_RELEASE_NULL	((Notify_release)0)
#define	NOTIFY_ARG_NULL		((Notify_arg)0)

/*
 * Mask bit generating macros 	(for prioritizer):
 */
#define	SIG_BIT(sig)		(1 << ((sig)-1))

/*
 ***********************************************************************
 *		Typedefs, enumerations, and structures
 ***********************************************************************
 */

/*
 * Opaque client handle.
 */
typedef	Xv_opaque Notify_client;

#ifndef _NOTIFY_MIN_SYMBOLS	/* Hack to reduce symbols in libraries. */

/*
 * Opaque client event.
 */
typedef	Xv_opaque Notify_event;

/*
 * Opaque client event optional argument.
 */
typedef	Xv_opaque Notify_arg;

/*
 * A pointer to a function returning a Notify_arg (used for client
 * event additional argument copying).
 */
typedef	Notify_arg (*Notify_copy)();

/*
 * A pointer to a function returning void (used for client
 * event additional argument storage releasing).
 */
typedef	void (*Notify_release)();

/*
 * For Debugging utility:
 */
typedef	enum notify_dump_type {
	NOTIFY_ALL=0,
	NOTIFY_DETECT=1,
	NOTIFY_DISPATCH=2
} Notify_dump_type;

#endif /* ~_NOTIFY_MIN_SYMBOLS */

/*
 * Client notification function return values for notifier to client calls.
 */
typedef enum notify_value {
	NOTIFY_DONE		= 0,	/* Handled notification */
	NOTIFY_IGNORED		= 1,	/* Did nothing about notification */
	NOTIFY_UNEXPECTED	= 2	/* Notification not expected */
} Notify_value;

/*
 * A pointer to a function returning a Notify_value.
 */
typedef	Notify_value (*Notify_func)();

/*
 * Error codes for client to notifier calls (returned when no other
 * return value or stored in notify_errno when other return value
 * indicates error condition).
 */
typedef enum notify_error {
	NOTIFY_OK		= 0,	/* Success */
	NOTIFY_UNKNOWN_CLIENT	= 1,	/* Client argument unknown to notifier */
	NOTIFY_NO_CONDITION	= 2,	/* Client not registered for given 
					 * condition 
					 */
	NOTIFY_BAD_ITIMER	= 3,	/* Itimer type unknown */
	NOTIFY_BAD_SIGNAL	= 4,	/* Signal number out of range */
	NOTIFY_NOT_STARTED	= 5,	/* Notify_stop called & notifier not 
					 * started 
					 */
	NOTIFY_DESTROY_VETOED	= 6,	/* Some client didn't want to die when 
					 * called notify_die(DESTROY_CHECKING)
					 */
	NOTIFY_INTERNAL_ERROR	= 7,	/* Something wrong in the notifier */
	NOTIFY_SRCH		= 8,	/* No such process */
	NOTIFY_BADF		= 9,	/* Bad file number */
	NOTIFY_NOMEM		= 10,	/* Not enough core */
	NOTIFY_INVAL		= 11,	/* Invalid argument */
	NOTIFY_FUNC_LIMIT	= 12	/* Too many interposition functions */
} Notify_error;


/*
 * Argument types
 */
typedef enum notify_signal_mode {
	NOTIFY_SYNC		= 0,
	NOTIFY_ASYNC		= 1
} Notify_signal_mode;

typedef enum notify_event_type {
	NOTIFY_SAFE		= 0,
	NOTIFY_IMMEDIATE	= 1
} Notify_event_type;

typedef enum destroy_status {
	DESTROY_PROCESS_DEATH	= 0,
	DESTROY_CHECKING	= 1,
	DESTROY_CLEANUP		= 2,
	DESTROY_SAVE_YOURSELF	= 3
} Destroy_status;

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

/*
 * PUBLIC variables 
 */

extern	struct itimerval NOTIFY_POLLING_ITIMER;	/* {{0,1},{0,1}} */

/*
 * PRIVATE variables 
 */

extern	Notify_error notify_errno;

#ifndef _NOTIFY_MIN_SYMBOLS
extern	struct itimerval NOTIFY_NO_ITIMER;	/* {{0,0},{0,0}} */
#endif	/* ~_NOTIFY_MIN_SYMBOLS */

/*
 * PUBLIC functions 
 */

#ifndef _NOTIFY_MIN_SYMBOLS

#ifdef SYSV_WAIT
EXTERN_FUNCTION (Notify_value   notify_default_wait3, (Notify_client nclient, int pid, int  *status, struct rusage *rusage));
#else
EXTERN_FUNCTION (Notify_value 	notify_default_wait3, (Notify_client nclient, int pid, union wait *status, struct rusage *rusage));
#endif 

EXTERN_FUNCTION (Notify_error 	notify_dispatch, (void));
EXTERN_FUNCTION (Notify_error	notify_do_dispatch, (void));
EXTERN_FUNCTION (Notify_error 	notify_itimer_value, (Notify_client nclient, int which, struct itimerval *value));
EXTERN_FUNCTION (Notify_value 	notify_next_destroy_func, (Notify_client nclient, Destroy_status status));
EXTERN_FUNCTION (Notify_value 	notify_next_event_func, (Notify_client nclient, Notify_event event, Notify_arg arg, Notify_event_type when));
EXTERN_FUNCTION (Notify_error 	notify_no_dispatch, (void));
EXTERN_FUNCTION (Notify_func 	notify_set_destroy_func, (Notify_client nclient, Notify_func func));
EXTERN_FUNCTION (Notify_func 	notify_set_exception_func, (Notify_client nclient, Notify_func func, int fd));
EXTERN_FUNCTION (Notify_func 	notify_set_input_func, (Notify_client nclient, Notify_func func, int fd));
EXTERN_FUNCTION (Notify_func 	notify_set_itimer_func, (Notify_client nclient, Notify_func func, int which, struct itimerval *value, struct itimerval *ovalue));
EXTERN_FUNCTION (Notify_func 	notify_set_output_func, (Notify_client nclient,	 Notify_func func, int fd));
EXTERN_FUNCTION (Notify_func 	notify_set_signal_func, (Notify_client nclient,	Notify_func func, int sig, Notify_signal_mode mode));
EXTERN_FUNCTION (Notify_func 	notify_set_wait3_func, 	(Notify_client nclient, Notify_func func, int pid));
EXTERN_FUNCTION (Notify_error 	notify_start, (void));
EXTERN_FUNCTION (Notify_error 	notify_stop, (void));
EXTERN_FUNCTION (Notify_error 	notify_veto_destroy, (Notify_client nclient));
EXTERN_FUNCTION (void 		notify_perror, (char *str));
EXTERN_FUNCTION (void 		notify_enable_rpc_svc, (int flag));

#endif	/* ~_NOTIFY_MIN_SYMBOLS */

/*
 * PRIVATE functions 
 */

EXTERN_FUNCTION (Notify_func 	notify_set_event_func, (Notify_client nclient, Notify_func func, Notify_event_type when));
EXTERN_FUNCTION (Notify_error 	notify_remove, (Notify_client nclient));

#ifndef _NOTIFY_MIN_SYMBOLS

EXTERN_FUNCTION (Notify_error 	notify_client, (Notify_client nclient));
EXTERN_FUNCTION (Notify_error 	notify_destroy, (Notify_client nclient, Destroy_status status));
EXTERN_FUNCTION (Notify_error 	notify_die, (Destroy_status status));
EXTERN_FUNCTION (Notify_error 	notify_event, (Notify_client nclient, Notify_event event, Notify_arg arg));
EXTERN_FUNCTION (Notify_error 	notify_exception, (Notify_client nclient, int fd));
EXTERN_FUNCTION (void 		notify_flush_pending, (Notify_client nclient));
EXTERN_FUNCTION (Notify_func 	notify_get_destroy_func, (Notify_client));
EXTERN_FUNCTION (Notify_func 	notify_get_event_func, (Notify_client nclient, Notify_event_type when));
EXTERN_FUNCTION (Notify_func 	notify_get_exception_func, (Notify_client nclient, int fd));
EXTERN_FUNCTION (Notify_func	notify_get_input_func, (Notify_client nclient, int fd));
EXTERN_FUNCTION (Notify_func	notify_get_itimer_func,	(Notify_client nclient, int which));
EXTERN_FUNCTION (Notify_func 	notify_get_output_func, (Notify_client nclient, int fd));
EXTERN_FUNCTION (Notify_func 	notify_get_prioritizer_func, (Notify_client nclient));
EXTERN_FUNCTION (Notify_func 	notify_get_scheduler_func, (void));
EXTERN_FUNCTION (int 		notify_get_signal_code, (void));

#ifndef SYSV_UCONTEXT
EXTERN_FUNCTION (struct sigcontext *notify_get_signal_context, (void));
#else 
EXTERN_FUNCTION (ucontext_t *notify_get_signal_context, (void));
#endif

EXTERN_FUNCTION (Notify_func 	notify_get_signal_func, (Notify_client nclient, int signal, Notify_signal_mode mode));
EXTERN_FUNCTION (Notify_func 	notify_get_wait3_func, (Notify_client nclient, int pid));
EXTERN_FUNCTION (Notify_error	notify_input, (Notify_client nclient, int fd));
EXTERN_FUNCTION (Notify_error	notify_interpose_destroy_func, (Notify_client nclient, Notify_func func));
EXTERN_FUNCTION (Notify_error 	notify_interpose_event_func, (Notify_client nclient, Notify_func func, Notify_event_type when));
EXTERN_FUNCTION (Notify_error	notify_interpose_exception_func, (Notify_client nclient, Notify_func func, int fd));
EXTERN_FUNCTION (Notify_error	notify_interpose_input_func, (Notify_client nclient, Notify_func func, int fd));
EXTERN_FUNCTION (Notify_error	notify_interpose_itimer_func, (Notify_client nclient, Notify_func func, int which));
EXTERN_FUNCTION (Notify_error 	notify_interpose_output_func, (Notify_client nclient, Notify_func func, int fd));
EXTERN_FUNCTION (Notify_error 	notify_interpose_signal_func, (Notify_client nclient, Notify_func func, int signal, Notify_signal_mode mode));
EXTERN_FUNCTION (Notify_error	notify_interpose_wait3_func, (Notify_client nclient, Notify_func func, int pid));
EXTERN_FUNCTION (Notify_error 	notify_itimer, 	(Notify_client nclient, int which));
EXTERN_FUNCTION (Notify_value 	notify_next_exception_func, (Notify_client nclient, int fd));
EXTERN_FUNCTION (Notify_value 	notify_next_input_func, (Notify_client nclient, int fd));
EXTERN_FUNCTION (Notify_value 	notify_next_itimer_func, (Notify_client nclient, int which));
EXTERN_FUNCTION (Notify_value	notify_next_output_func, (Notify_client nclient, int fd));
EXTERN_FUNCTION (Notify_value 	notify_next_signal_func, (Notify_client nclient, int signal, Notify_signal_mode mode));

#ifdef SYSV_WAIT
EXTERN_FUNCTION (Notify_value 	notify_next_wait3_func, (Notify_client nclient, int pid, int  *status, struct rusage *rusage));
#else
EXTERN_FUNCTION (Notify_value 	notify_next_wait3_func, (Notify_client nclient, int pid, union wait *status, struct rusage *rusage));
#endif 

EXTERN_FUNCTION (Notify_value	notify_nop, (void));
EXTERN_FUNCTION (Notify_error 	notify_output, (Notify_client nclient, int fd));
EXTERN_FUNCTION (Notify_error	notify_post_destroy, (Notify_client nclient, Destroy_status status, Notify_event_type type));
EXTERN_FUNCTION (Notify_error 	notify_post_event, (Notify_client nclient, Notify_event event, 	Notify_event_type when_hint));

/* vmh - 10/15/90: one argument was missing */
EXTERN_FUNCTION (Notify_error 	notify_post_event_and_arg, (Notify_client nclient, Notify_event event, Notify_event_type when_hint, Notify_arg arg, Notify_copy copy_func, Notify_release release_func));

EXTERN_FUNCTION (Notify_error 	notify_remove_destroy_func, (Notify_client nclient, Notify_func func));
EXTERN_FUNCTION (Notify_error 	notify_remove_event_func, (Notify_client nclient, Notify_func func, Notify_event_type when));
EXTERN_FUNCTION (Notify_error 	notify_remove_exception_func, (Notify_client nclient, Notify_func func, int fd));
EXTERN_FUNCTION (Notify_error	notify_remove_input_func, (Notify_client nclient, Notify_func func, int fd));
EXTERN_FUNCTION (Notify_error 	notify_remove_itimer_func, (Notify_client nclient, Notify_func func, int which));
EXTERN_FUNCTION (Notify_error	notify_remove_output_func, (Notify_client nclient, Notify_func func, int fd));
EXTERN_FUNCTION (Notify_error 	notify_remove_signal_func, (Notify_client nclient, Notify_func func, int signal, Notify_signal_mode mode));
EXTERN_FUNCTION (Notify_error 	notify_remove_wait3_func, (Notify_client nclient, Notify_func func, int pid));
EXTERN_FUNCTION (Notify_func 	notify_set_prioritizer_func, (Notify_client nclient, Notify_func func));
EXTERN_FUNCTION (Notify_func 	notify_set_scheduler_func, (Notify_func nclient));
EXTERN_FUNCTION (Notify_error 	notify_signal, (Notify_client nclient, int sig));
EXTERN_FUNCTION (Notify_error 	notify_wait3, (Notify_client nclient));

extern	Notify_error	notify_errno;

/*
 * FD manipulation functions
 */

EXTERN_FUNCTION (int 		ntfy_fd_cmp_and, (fd_set *a, fd_set *b));
EXTERN_FUNCTION (int 		ntfy_fd_cmp_or, (fd_set *a, fd_set *b));
EXTERN_FUNCTION (int 		ntfy_fd_anyset, (fd_set *a));
EXTERN_FUNCTION (fd_set *	ntfy_fd_cpy_or, (fd_set *a, fd_set *b));
EXTERN_FUNCTION (fd_set *	ntfy_fd_cpy_and, (fd_set *a, fd_set *b));
EXTERN_FUNCTION (fd_set *	ntfy_fd_cpy_xor, (fd_set *a, fd_set *b));

/*
 * Debugging Utility 
 */

EXTERN_FUNCTION (void 		notify_dump, (Notify_client nclient, Notify_dump_type type, FILE * file));

#endif /* ~_NOTIFY_MIN_SYMBOLS */

#endif	/* ~xview_notify_DEFINED */
