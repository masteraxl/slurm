/* p2p locks.  This file includes code for picking a strategy,
 * whether generic (e.g., the SysV SEMOPS) or machine specific.
 * First we choose a method, and then, in a second step, generate the
 * appropriate declarations.
 */

/* Choose method for providing locks.  The decision will be made by
 * combining preferences (defined here) with available capabilities (should be,
 * and for the most part are, defined by configure as HAVE_, al).  Some of the
 * HAVE_'s may be defined in p2p_special.h, at least temporarily.
 */

#if !defined(LOCKS_PICKED) 

/* Prefer vendor-specific locks on IRIX */
#if defined(MPI_IRIX)
#    define PREFER_USLOCKS
#    define PREFER_SPINLOCKS	/* prefer spinlocks to uslock semaphores */
#endif

/* Prefer vendor-specific locks on HPUX */
#if defined(MPI_HPUX)
#    define PREFER_HPLOCKS
#endif

/* prefer assembly language locks to SEMOPs on SX-4 */
#if defined(MPI_SX_4)
#    define PREFER_TSLOCKS
#endif

/* Prefer MUTEX locks (currently only on Solaris) to SEMOPS */
#if defined(HAVE_MUTEX_INIT)
#    define PREFER_MUTEX
#endif

/* Prefer MSEMAPHORES (currently at least on AIX), except on HP*/
#if defined(HAVE_MSEM_INIT)
#    if !defined(HP_UX)
#        define PREFER_MSEM
#    endif
#endif

/* Now we convert our preferences and our capabilites into choices of what
 * to actually use.  Some cumbersomeness is generated by checking that we
 * do not pick two different schemes.
 */

#ifdef LOCKS_PICKED
#    undef LOCKS_PICKED		/* start clean, no final decisions so far */
#endif

#if defined(HAVE_USLOCKS) && defined(PREFER_USLOCKS)
#    if defined(LOCKS_PICKED)
         'Oops - trying to use two different locking schemes'
#    else
#        define USE_USLOCKS
#        define LOCKS_PICKED
#    endif
#endif

#if defined(HAVE_HPLOCKS) && defined(PREFER_HPLOCKS)
#    if defined(LOCKS_PICKED)
         'Oops - trying to use two different locking schemes'
#    else
#        define USE_HPLOCKS
#        define LOCKS_PICKED
#    endif
#endif

#if defined(HAVE_TSLOCKS) && defined(PREFER_TSLOCKS)
#    if defined(LOCKS_PICKED)
         'Oops - trying to use two different locking schemes'
#    else
#        define USE_TSLOCKS
#        define LOCKS_PICKED
#    endif
#endif

#if defined(HAVE_MSEM_INIT) && defined(PREFER_MSEM)
#    if defined(LOCKS_PICKED)
         'Oops - trying to use two different locking schemes'
#    else
#        define USE_MSEM
#        define LOCKS_PICKED
#    endif
#endif

#if defined(HAVE_MUTEX_INIT) && defined(PREFER_MUTEX)
#    if defined(LOCKS_PICKED)
         'Oops - trying to use two different locking schemes'
#    else
#        define USE_MUTEX
#        define LOCKS_PICKED
#    endif
#endif

/* Default: SysV semaphores, unfortunately */
#if !defined(LOCKS_PICKED) && defined(HAVE_SEMOP) && defined(HAVE_SEMGET)
#    define USE_SEMOP
#    define LOCKS_PICKED
#endif

#if !defined(LOCKS_PICKED) && defined(HAVE_SEM_POST) 
#    define USE_POSIX_SEM
#    define LOCKS_PICKED
#endif

/* Check that we picked at least one type of locks.  We have already
 * checked that we did not pick more than one.
 */
#if !defined(LOCKS_PICKED)
#error  Choke - no locks picked
#endif

#endif /* !defined(LOCKS_PICKED) */
/* At this point we should have generated exactly one USE_xxx.
 * Now we generate the appropriate includes and declarations.
 */

#if defined(USE_USLOCKS)
#    include <ulocks.h>
     /* usema_t is defined as "void" in /usr/include/ulocks.h */
     typedef usema_t           *p2p_lock_t;

#    if defined(PREFER_SPINLOCKS)
#        define p2p_lock_init(l)   (*(l)) = usnewlock(p2p_sgi_usptr)
#        define p2p_lock(l)     ussetlock(*(l))
#        define p2p_unlock(l)   usunsetlock(*(l))
#        define p2p_lock_name   "uslocks - spinlocks"
#    elif defined(PREFER_SEMAPHORES)
#        define p2p_lock_init(l)   (*(l)) = usnewsema(p2p_sgi_usptr,1)  
#        define p2p_lock(l)     uspsema(*(l))
#        define p2p_unlock(l)   usvsema(*(l))
#        define p2p_lock_name   "uslocks - semaphores"
#    else
         'Oops - no uslocks'
#    endif
#    define p2p_lock_free(l)
#endif

#if defined(USE_HPLOCKS)
    /* HPUX uses a special lock/unlock set (see mem.c file) given to us
       by Dan Golan of Convex */
    typedef int p2p_lock_t[4];
#    define p2p_lock_init(l)   { *((int*)(l)) = 1; }
#    define p2p_lock(l)        MPID_SHMEM__acquire_lock(l)
#    define p2p_unlock(l)      MPID_SHMEM__release_lock(l)
#    define p2p_lock_name      "HPUX assembly language locks"
#    define p2p_lock_free(l)
#endif

#if defined(USE_TSLOCKS)
    /* Pass the routines a (p2p_lock_t *) */
    typedef long p2p_lock_t;
#    define p2p_lock_init(l)   tslock_init(l)
#    define p2p_lock(l)        tslock(l)
#    define p2p_unlock(l)      tsunlock(l)
#    define p2p_lock_name      "tslocks"
#    define p2p_lock_free(l)
#endif

#if defined(USE_MSEM)
#include <sys/mman.h>
    /* Problem - some systems (e.g., HP) use struct msemaphore, others use just
                 msemaphore */
#    if defined(MSEMAPHORE_IS_STRUCT)   /* HP-style */
        typedef struct msemaphore  MPID_msemaphore;
#    else
        typedef msemaphore         MPID_msemaphore;
#    endif

    /* Place each lock on its own cache line.  We probably really want two
       lock types - one on its own line, and one within another structure */
    typedef struct { MPID_msemaphore lock; 
	            char pad[MPID_CACHE_LINE_SIZE - sizeof(MPID_msemaphore)]; }
                               p2p_lock_t;
    /* An alternative to this is to ALLOCATE these on separate cache lines,
       but this is simpler for now */

#    define p2p_lock_init(l)   msem_init(&(l)->lock, MSEM_UNLOCKED)
#    define p2p_lock(l)        msem_lock(&(l)->lock, 0)
#    define p2p_unlock(l )     msem_unlock(&(l)->lock, 0)
#    define p2p_lock_name      "msem_lock"

    /* Non-cache-line-separated locks:
    typedef msemaphore         p2p_lock_t;
#    define p2p_lock_init(l)   msem_init(l, MSEM_UNLOCKED)
#    define p2p_lock(l)        msem_lock(l, 0)
#    define p2p_unlock(l )     msem_unlock(l, 0)
    */
#    define p2p_lock_free(l)
#endif

#if defined(USE_MUTEX)
/*   Only known system is Solaris */
#    include <sys/systeminfo.h>
#    include <sys/processor.h>
#    include <sys/procset.h>
#    include <synch.h>
     typedef mutex_t           MPID_msemaphore;
     typedef mutex_t           p2p_lock_t;
#    define p2p_lock_init(l)   mutex_init(l,USYNC_PROCESS,(void *)NULL)
#    define p2p_lock(l)        mutex_lock(l)
#    define p2p_unlock(l)      mutex_unlock(l)
#    define p2p_lock_name      "mutex_lock"
#    define p2p_lock_free(l)
#endif

#if defined(USE_SEMOP)
#    include <sys/types.h>
#    include <sys/ipc.h>
#    include <sys/sem.h>
#    define INCLUDED_SYS_SEM
     typedef struct { int semid;  int semnum; }   p2p_lock_t;
#    define p2p_lock_init(l)  semop_init(l)
#    define p2p_lock(l)       semop_lock(l)
#    define p2p_unlock(l)     semop_unlock(l)
#    define p2p_lock_name     "semop_lock"
void semop_init (p2p_lock_t *);
void semop_lock (p2p_lock_t *);
void semop_unlock (p2p_lock_t *);
int MD_init_semop( void );
int MD_init_sysv_semop(void);
void MD_remove_sysv_sipc( void );
#    define p2p_lock_free(l)
#endif

#if defined(USE_POSIX_SEM)
#include <?>
#define INCLUDED_POSIX_SEM
typedef struct { sem_t mutex; } p2p_lock_t;
#define p2p_lock_init(l) sem_init( &((l).mutex), 0, 1 )
/* Sem_wait is a slow system call and can be interrupted */
#define p2p_lock(l)      \
    while (sem_wait( &((l).mutex) == -1 && errno == EAGAIN );
#define p2p_unlock(l)    sem_post( &((l).mutex) )
#define p2p_lock_name    "Posix sem"
#define p2p_lock_free(l) sem_destroy( &((l).mutex )
#endif

/* A few odds and ends */

/* Putting addresses of locks on separate cache lines */
/* This is probably NOT a good idea, since the standard locks are
 * retrieved as addresses by the usnewlock routine.
 * However, this DOES let the ADDRESSES of the locks reside on 
 * different cachelines.
 */
/* Note currently inactivated */
#if 0 && defined(MPID_CACHE_ALIGN)
#  if defined(MPI_IRIX)
     typedef usema_t *MD_lock_t;
     typedef union { MD_lock_t lock; PAD(MPID_CACHE_LINE_SIZE) } p2p_lock_t;
#    define p2p_lock_init(l)   (*l).lock = usnewlock(p2p_sgi_usptr)
#    define p2p_lock(l)        MD_lock((*l).lock)
#    define p2p_unlock(l)      MD_unlock((*l).lock)
#    define MD_lock(l)         ussetlock((l))
#    define MD_unlock(l)       usunsetlock((l))
#  endif
#endif

/* 
 * p2p_write_sync forces writes to be written to cache-coherent memory.
 * Some processors have special, assembly-language instructions for this.
 * Otherwise, you can usually do a lock/unlock.
 *
 * This has been made specific to MPID; this should really use
 * a lock initialized in p2p_init, with a p2p_write_sync_init()
 * macro used to do the lock initialization.
 */
#define p2p_write_sync() {\
     p2p_lock( &MPID_shmem->globlock );\
     p2p_unlock( &MPID_shmem->globlock );\
     }
/*
 * Define some condition variable operations.  To allow both in-process
 * spin-on-value and more general OS wait-for-change, these operations
 * contain both the local condition being changed (e.g., assigning to a 
 * pointer) and a condition variable.
 *
 * These are not yet used in the code.  They may eventually be used in the
 * MPID_SHMEM_ReadControl and _SendControl routines.
 */

#ifdef USE_SEM_CONDVAR
#define p2p_condvar_init(l) sem_init(&l)
#define p2p_condvar_free(l) sem_free(&l)
#define p2p_condvar_decl(l) sem_t l
#define p2p_condvar_post(exp,l) exp;sem_post(&l)
#define p2p_condvar_wait(boolexp,l) sem_wait(&l)
#endif

#ifndef p2p_condvar_init
     /* No condition variables.  Instead, use only the expressions */
#define p2p_condvar_init(l)
#define p2p_condvar_free(l)
#define p2p_condvar_decl(l)
#define p2p_condvar_post(exp,condvar) exp
     /* The wait needs a better "do while spinning" operation.
        Instead, should the API be 
	  p2p_condvar_wait(boolexp,condvar,optwhilewait)
	where optwhilewait is what is executed while spinning?
	Then when using condition variables, we could still do
	something like
	  if (!exp) {optwhilewait;} sem_wait(condvar)
	which would allow a implementation to use the "optwhilewait"
	to flush any pending stuff before entering the blocking call
     */
#define p2p_condvar_wait(boolexp,condvar) while (!(boolexp)) 
#endif