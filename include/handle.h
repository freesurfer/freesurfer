#ifndef HANDLE_H
#define HANDLE_H

typedef unsigned long PTR_HANDLE ;

PTR_HANDLE HandleAlloc(void *ptr) ;
void       HandleFree(PTR_HANDLE handle) ;
void       *HandleToPtr(PTR_HANDLE handle) ;
int        HandleOk(PTR_HANDLE handle) ;


#endif
