#ifndef LWP_PROTO_H
#define LWP_PROTO_H


#include <stdarg.h>
#include <lwp/lwp.h>
#include <lwp/stackdep.h>
#include <lwp/lwpmachdep.h>

int lwp_checkstkset(thread_t tid, caddr_t limit) ;
int lwp_stkcswset(thread_t tid, caddr_t limit) ;
int lwp_setstkcache(int minstksz, int numstks) ;
stkalign_t *lwp_newstk(void) ;
stkalign_t *lwp_datastk(caddr_t data, int size, caddr_t * addr) ;
int lwp_yield(thread_t tid) ;
int lwp_sleep(struct timeval *timeout) ;
int lwp_resched(int prio) ;
int lwp_setpri(thread_t tid, int prio) ;
int lwp_suspend(thread_t tid) ;
int lwp_resume(thread_t tid) ;
int lwp_join(thread_t tid) ;
int lwp_create(thread_t *tid, void (*func)(int tid, void *parm), int prio, 
							 int flags, stkalign_t *stack, int nargs, ...) ;
int lwp_destroy(thread_t tid) ;
void pod_setexit(int status) ;
int pod_getexit(int status) ;
void pod_exit(int status) ;
int pod_getmaxpri(void) ;
int pod_getmaxsize(void) ;
int pod_setmaxpri(int maxprio) ;

int lwp_enumerate(thread_t vec[], int maxsize) ;
int lwp_ping(thread_t tid) ;
int lwp_getregs(thread_t tid, machstate_t *machstate) ;
int lwp_setregs(thread_t tid, machstate_t *machstate) ;
int lwp_getstate(thread_t tid, statvec_t *statvec) ;
int lwp_self(thread_t *tid) ;

#endif
