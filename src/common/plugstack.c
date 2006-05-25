/*****************************************************************************\
 *  plugstack.c -- stackable plugin architecture for node job kontrol (SPANK)
 *****************************************************************************
 *  Copyright (C) 2005 The Regents of the University of California.
 *  Produced at Lawrence Livermore National Laboratory (cf, DISCLAIMER).
 *  UCRL-CODE-217948.
 *
 *  This file is part of SLURM, a resource management program.
 *  For details, see <http://www.llnl.gov/linux/slurm/>.
 *
 *  SLURM is free software; you can redistribute it and/or modify it under
 *  the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  SLURM is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with SLURM; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
\*****************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string.h>

#include "src/common/plugin.h"
#include "src/common/xmalloc.h"
#include "src/common/xstring.h"
#include "src/common/xassert.h"
#include "src/common/safeopen.h"
#include "src/common/strlcpy.h"
#include "src/common/read_config.h"
#include "src/common/plugstack.h"
#include "src/common/optz.h"

#include "src/slurmd/slurmstepd/slurmstepd_job.h"

#include <slurm/spank.h>

#define REQUIRED "required"
#define OPTIONAL "optional"

struct spank_plugin_operations {
    spank_f * init;
    spank_f * user_init;
    spank_f * user_task_init;
    spank_f * task_post_fork;
    spank_f * task_exit;
    spank_f * exit;
};

const int n_spank_syms = 6;
const char * spank_syms [] = {
    "slurm_spank_init",
    "slurm_spank_user_init",
    "slurm_spank_task_init",
    "slurm_spank_task_post_fork",
    "slurm_spank_task_exit",
    "slurm_spank_exit"
};

struct spank_plugin {
    const char *name;
    char *fq_path;
    plugin_handle_t plugin;
    bool required;
    int ac;
    char **argv;
    struct spank_plugin_operations ops;
    struct spank_option *opts;
};

/*
 *  SPANK Plugin options 
 */
struct spank_plugin_opt {
    struct spank_option * opt;       /* Copy of plugin option info           */
    struct spank_plugin * plugin;    /* Link back to plugin structure        */
    int                   optval;    /* Globally unique value                */
    int                   found:1;   /* 1 if option was found, 0 otherwise   */
    int                   disabled:1;/* 1 if option is cached but disabled   */
    char *                optarg;    /* Option argument.                     */
};

/*
 *  Initial value for global optvals for SPANK plugin options
 */
static int spank_optval = 0xfff;

/*
 *  Cache of options provided by spank plugins
 */
static List option_cache = NULL;


/*
 *  SPANK handle for plugins
 */
struct spank_handle {
#   define SPANK_MAGIC 0x00a5a500
    int magic;
    slurmd_job_t *job;
    slurmd_task_info_t *task;
};


/*
 *  SPANK plugin hook types:
 */
typedef enum step_fn {
    STEP_INIT = 0,
    STEP_USER_INIT,
    STEP_USER_TASK_INIT,
    STEP_TASK_POST_FORK,
    STEP_TASK_EXIT,
    STEP_EXIT
} step_fn_t;

/*
 *  SPANK plugins stack
 */
static List spank_stack = NULL;

static pthread_mutex_t spank_mutex = PTHREAD_MUTEX_INITIALIZER;


/*
 *  Forward declarations
 */
static int _spank_plugin_options_cache (struct spank_plugin *p);


static void _argv_append (char ***argv, int ac, const char *newarg)
{
    *argv = xrealloc (*argv, (++ac + 1) * sizeof (char *));
    (*argv)[ac] = NULL;
    (*argv)[ac-1] = xstrdup (newarg);
    return;
}

static int
_plugin_stack_parse_line (char *line, char **plugin, int *acp, char ***argv,
    bool *required)
{
    int                 ac;
    const char *        separators = " \t\n";
    char *              path;
    char *              option;
    char *              s;
    char **             av;
    char *              sp;

    *plugin = NULL;
    *argv = NULL;
    *acp = 0;

    /* Nullify any comments
     */
    if ((s = strchr (line, '#')))
        *s = '\0';

    if (!(option = strtok_r (line, separators, &sp)))
        return 0;

    if (strncmp (option, REQUIRED, strlen (option)) == 0) {
        *required = true;
    } else if (strncmp (option, OPTIONAL, strlen (option)) == 0) {
        *required = false;
    } else {
        error ("spank: Invalid option \"%s\". Must be either %s or %s",
            option, REQUIRED, OPTIONAL);
        return (-1);
    }

    if (!(path = strtok_r (NULL, separators, &sp)))
        return (-1);

    ac = 0;
    av = NULL;

    while ((s = strtok_r (NULL, separators, &sp)))
        _argv_append (&av, ac++, s);

    *plugin = xstrdup (path);
    *argv  = av;
    *acp   = ac;

    return (0);
}

static struct spank_plugin *
_spank_plugin_create (char *path, int ac, char **av, bool required)
{
    struct spank_plugin *plugin;
    plugin_handle_t p;
    struct spank_plugin_operations ops;

    if (!(p = plugin_load_from_file (path))) 
        return NULL;

    if (plugin_get_syms (p, n_spank_syms, spank_syms, (void **) &ops) == 0) {
        error ("spank: \"%s\" exports 0 symbols\n", path);
        return NULL;
    }

    plugin = xmalloc (sizeof (struct spank_plugin));

    plugin->fq_path =  path; /* fq_path is xstrdup'd in *process_line */
    plugin->plugin =   p;
    plugin->name =     plugin_get_name (p); /* no need to dup */
    plugin->required = required;
    plugin->ac =       ac;
    plugin->argv =     av;
    plugin->ops =      ops;

    plugin->opts =     plugin_get_sym (p, "spank_options");

    return (plugin);
}

void _spank_plugin_destroy (struct spank_plugin *sp)
{
    if (sp == NULL)
        return;

    xfree (sp->fq_path);

    /* No need to free "name" it was defined within plugin */
    sp->name = NULL;

    plugin_unload (sp->plugin);
    sp->plugin = NULL;
    if (sp->argv) {
        int i;
        for (i = 0; sp->argv[i]; i++) 
            xfree (sp->argv[i]);
        xfree (sp->argv);
    }
    xfree (sp);
    return;
}

static int
_spank_stack_process_line (const char *file, int line, char *buf,
        struct spank_plugin **plugin)
{
    char **argv;
    int    ac;
    char  *path;
    bool   required;

    struct spank_plugin *p;

    *plugin = NULL;

    if (_plugin_stack_parse_line (buf, &path, &ac, &argv, &required) < 0) {
        error ("spank: %s: %d: Invalid line. Ignoring.", file, line);
        return (0);
    }

    if (path == NULL) /* No plugin listed on this line */
        return (0);

    if (!(p = _spank_plugin_create (path, ac, argv, required))) {
        error ("spank: %s: %d: Failed to load %s plugin from %s. %s",
               file, line, required ? "required" : "optional",
               path, required ? "Aborting." : "Ignoring.");
        return (required ? -1 : 0);
    }

    *plugin = p;

    return (0);
}

static int
_spank_stack_create (const char *path, List *listp)
{
    int   line;
    char  buf [4096];
    FILE *fp;

    *listp = NULL;

    verbose ("spank: opening plugin stack %s\n", path);

    if (!(fp = safeopen (path, "r", SAFEOPEN_NOCREATE)))
        return -1;

    line = 1;
    while (fgets (buf, sizeof (buf), fp)) {
        struct spank_plugin *p;

        if (_spank_stack_process_line (path, line, buf, &p) < 0)
            goto fail_immediately;

        if (p == NULL)
            continue;

        if (*listp == NULL)
            *listp = list_create ((ListDelF) _spank_plugin_destroy);

        verbose ("spank: loaded plugin %s\n", xbasename (p->fq_path));
        list_append (*listp, p);

        _spank_plugin_options_cache (p);

        line++;
    }

    return (0);

    fail_immediately:
    if (*listp != NULL) {
        list_destroy (*listp);
        *listp = NULL;
    }
    return (-1);
}

static int
_spank_handle_init (struct spank_handle *spank, slurmd_job_t *job, int taskid)
{
    memset (spank, 0, sizeof (*spank));
    spank->magic = SPANK_MAGIC;
    spank->job   = job;
    if (taskid >= 0)
        spank->task = job->task [taskid];
    return (0);
}

static const char *
_step_fn_name (step_fn_t type)
{
    switch (type) {
        case STEP_INIT:
            return ("init");
        case STEP_USER_INIT:
            return ("user_init");
        case STEP_USER_TASK_INIT:
            return ("task_init");
        case STEP_TASK_POST_FORK:
            return ("task_post_fork");
    case STEP_TASK_EXIT:
            return ("task_exit");
        case STEP_EXIT:
            return ("exit");
    }

    /* NOTREACHED */
    return ("unknown");
}

static int
_do_call_stack (step_fn_t type, slurmd_job_t *job, int taskid)
{
    int                   rc = 0;
    ListIterator          i;
    struct spank_plugin * sp;
    struct spank_handle   spank[1];
    const char *          fn_name;

    if (!spank_stack)
        return (0);

    if (_spank_handle_init (spank, job, taskid) < 0) {
        error ("spank: Failed to initialize handle for plugins");
        return (-1);
    }

    fn_name = _step_fn_name (type);

    i = list_iterator_create (spank_stack);
    while ((sp = list_next (i))) {
        const char * name = xbasename (sp->fq_path);

        switch (type) {
        case STEP_INIT:
            if (sp->ops.init) {
                rc = (*sp->ops.init) (spank, sp->ac, sp->argv);
                debug2 ("spank: %s: %s = %d\n", name, fn_name, rc);
            }
            break;
        case STEP_USER_INIT:
            if (sp->ops.user_init) {
                rc = (*sp->ops.user_init) (spank, sp->ac, sp->argv);
                debug2 ("spank: %s: %s = %d\n", name, fn_name, rc);
            }
            break;
        case STEP_USER_TASK_INIT:
            if (sp->ops.user_task_init) {
                rc = (*sp->ops.user_task_init) (spank, sp->ac, sp->argv);
                debug2 ("spank: %s: %s = %d\n", name, fn_name, rc);
            }
            break;
         case STEP_TASK_POST_FORK:
            if (sp->ops.task_post_fork) {
                rc = (*sp->ops.task_post_fork) (spank, sp->ac, sp->argv);
                debug2 ("spank: %s: %s = %d\n", name, fn_name, rc);
            }
            break;
         case STEP_TASK_EXIT:
            if (sp->ops.task_exit) {
                rc = (*sp->ops.task_exit) (spank, sp->ac, sp->argv);
                debug2 ("spank: %s: %s = %d", name, fn_name, rc);
            }
            break;
        case STEP_EXIT:
            if (sp->ops.exit) {
                rc = (*sp->ops.exit) (spank, sp->ac, sp->argv);
                debug2 ("spank: %s: %s = %d\n", name, fn_name, rc);
            }
            break;
       }

        if ((rc < 0) && sp->required) {
            error ("spank: required plugin %s: %s() failed with rc=%d", 
                    name, fn_name, rc);
            break;
    } else
            rc = 0;
    }

    list_iterator_destroy (i);

    return (rc);
}

int spank_load (void)
{
    slurm_ctl_conf_t *conf = slurm_conf_lock ();
    const char *path = conf->plugstack;
    slurm_conf_unlock ();

    if (_spank_stack_create (path, &spank_stack) < 0) {
        if (errno == ENOENT) /* No error if spank config doesn't exist */
            return (0);
        error ("spank: failed to create plugin stack: %m");
        return (-1);
    }

    return (0);
}

int spank_unload (void)
{
    if (spank_stack)
        list_destroy (spank_stack);
    if (option_cache)
        list_destroy (option_cache);

    return (0);
}

int spank_init (slurmd_job_t *job)
{
    if ((spank_stack == NULL) && (spank_load () < 0))
        return (-1);

    return (_do_call_stack (STEP_INIT, job, -1));
}

int spank_user (slurmd_job_t *job)
{
    return (_do_call_stack (STEP_USER_INIT, job, -1));
}

int spank_user_task (slurmd_job_t *job, int taskid)
{
    return (_do_call_stack (STEP_USER_TASK_INIT, job, taskid));
}

int spank_task_post_fork (slurmd_job_t * job, int taskid)
{
    return (_do_call_stack (STEP_TASK_POST_FORK, job, taskid));
}

int spank_task_exit (slurmd_job_t *job, int taskid)
{
    return (_do_call_stack (STEP_TASK_EXIT, job, taskid));
}

int spank_fini (slurmd_job_t *job)
{
    int rc = _do_call_stack (STEP_EXIT, job, -1);

    spank_unload ();

    return (rc);
}

/*
 *  SPANK options functions
 */

static int
_spank_next_option_val (void)
{
    int optval;
    slurm_mutex_lock (&spank_mutex);
    optval = spank_optval++;
    slurm_mutex_unlock (&spank_mutex);
    return (optval);
}

static struct spank_plugin_opt *
_spank_plugin_opt_create (struct spank_plugin *p, struct spank_option *opt,
        int disabled)
{
    struct spank_plugin_opt *spopt = xmalloc (sizeof (*spopt));
    spopt->opt =    opt;
    spopt->plugin = p;
    spopt->optval = _spank_next_option_val ();
    spopt->found =  0;
    spopt->optarg = NULL;

    spopt->disabled = disabled;

    return (spopt);
}

void
_spank_plugin_opt_destroy (struct spank_plugin_opt *spopt)
{
    xfree (spopt->optarg);
    xfree (spopt);
}

static int 
_opt_by_val (struct spank_plugin_opt *opt, int *optvalp)
{
    return (opt->optval == *optvalp);
}

static int 
_opt_by_name (struct spank_plugin_opt *opt, char *optname)
{
    return (strcmp (opt->opt->name, optname) == 0);
}

static int
_spank_plugin_options_cache (struct spank_plugin *p)
{
    int disabled = 0;
    struct spank_option *opt = p->opts;

    if ((opt == NULL) || opt->name == NULL)
        return (0);

    if (!option_cache) {
        option_cache = list_create ((ListDelF) _spank_plugin_opt_destroy);
    } 

    for (; opt && opt->name != NULL; opt++) {
        struct spank_plugin_opt *spopt;
        
        spopt = list_find_first (option_cache, (ListFindF) _opt_by_name, 
                                opt->name);
        if (spopt) {
            struct spank_plugin *q = spopt->plugin;
            info ("spank: option \"%s\" provided by both %s and %s", 
                   opt->name, xbasename (p->fq_path), xbasename (q->fq_path));
            /*
             *  Disable this option, but still cache it, in case
             *    options are loaded in a different order on the remote
             *    side.
             */
            disabled = 1;
        }

        if ((strlen (opt->name) > SPANK_OPTION_MAXLEN)) {
            error ("spank: option \"%s\" provided by %s too long. Ignoring.",
                   opt->name, p->name); 
            continue;
        }

        verbose ("SPANK: appending plugin option \"%s\"\n", opt->name);
        list_append (option_cache, _spank_plugin_opt_create (p, opt, disabled));
    }

    return (0);
}

struct option *
spank_option_table_create (const struct option *orig)
{
    struct spank_plugin_opt *spopt;
    struct option *opts = NULL;
    ListIterator i = NULL;

    opts = optz_create ();

    /*
     *  Start with original options:
     */
    if ((orig != NULL) && (optz_append (&opts, orig) < 0)) {
        optz_destroy (opts);
        return (NULL);
    }

    if (option_cache == NULL || (list_count (option_cache) == 0)) 
        return (opts);

    i = list_iterator_create (option_cache);
    while ((spopt = list_next (i))) {
        struct option opt;

        if (spopt->disabled)
            continue;

        opt.name = spopt->opt->name;
        opt.has_arg = spopt->opt->has_arg;
        opt.flag = NULL;
        opt.val = spopt->optval;

        if (optz_add (&opts, &opt) < 0) {
            if (errno == EEXIST) {
                error ("Ingoring conflicting option \"%s\" in plugin \"%s\"", 
                        opt.name, spopt->plugin->name);
            }
            else {
                error ("Unable to add option \"%s\" from plugin \"%s\"",
                        opt.name, spopt->plugin->name);
            }
            spopt->disabled = 1;
            continue;
        }

    }

    list_iterator_destroy (i);

    return (opts);
}

void 
spank_option_table_destroy (struct option *optz)
{
    optz_destroy (optz);
}

int 
spank_process_option (int optval, const char *arg)
{
    struct spank_plugin_opt *opt;
    int rc = 0;

    opt = list_find_first (option_cache, (ListFindF) _opt_by_val, &optval);

    if (!opt)
        return (-1);

    /*
     *  Call plugin callback if such a one exists
     */
    if (opt->opt->cb && (rc = ((*opt->opt->cb) (opt->opt->val, arg, 0))) < 0)
        return (rc);

    /*
     *  Set optarg and "found" so that option will be forwarded to remote side.
     */
    if (opt->opt->has_arg)
        opt->optarg = xstrdup (arg);
    opt->found = 1;

    return (0);
}

static void
_spank_opt_print (struct spank_option *opt, FILE *fp, int left_pad, int width)
{
    int    n;
    char * equals =    "";
    char * arginfo =   "";
    char   buf [81];


    if (opt->arginfo) {
        equals = "=";
        arginfo = opt->arginfo;
    }

    n = snprintf (buf, sizeof (buf), "%*s--%s%s%s", 
            left_pad, "", opt->name, equals, arginfo);

    if ((n < 0) || (n > sizeof (buf))) {
        const char trunc [] = "+";
        int len = strlen (trunc);
        char *p = buf + sizeof (buf) - len - 1;

        snprintf (p, len + 1, "%s", trunc);
    }

    if (n < width)
        fprintf (fp, "%-*s%s\n", width, buf, opt->usage); 
    else
        fprintf (fp, "\n%s\n%*s%s\n", buf, width, "", opt->usage);

    return;
}

int
spank_print_options (FILE *fp, int left_pad, int width)
{
    struct spank_plugin_opt *p;
    ListIterator i;

    if ((option_cache == NULL) || (list_count (option_cache) == 0))
        return (0);

    i = list_iterator_create (option_cache);
    while ((p = list_next (i))) {
        if (p->disabled)
            continue;
        _spank_opt_print (p->opt, fp, left_pad, width);
    }
    list_iterator_destroy (i);

    return (0);
}

#define OPT_TYPE_SPANK 0x4400

int 
spank_set_remote_options (job_options_t opts)
{
    struct spank_plugin_opt *p;
    ListIterator i;

    if ((option_cache == NULL) || (list_count (option_cache) == 0))
        return (0);

    i = list_iterator_create (option_cache);
    while ((p = list_next (i))) {
        char optstr [1024];

        if (!p->found)
            continue;

        snprintf (optstr, sizeof (optstr), "%s:%s", 
                  p->opt->name, p->plugin->name);

        job_options_append (opts, OPT_TYPE_SPANK, optstr, p->optarg);
    }
    list_iterator_destroy (i);
    return (0);
}

struct opt_find_args {
    const char *optname;
    const char *plugin_name;
};

static int _opt_find (struct spank_plugin_opt *p, struct opt_find_args *args)
{
    if (strcmp (p->plugin->name, args->plugin_name) != 0)
        return (0);
    if (strcmp (p->opt->name, args->optname) != 0)
        return (0);
    return (1);
}

static struct spank_plugin_opt *
_find_remote_option_by_name (const char *str)
{
    struct spank_plugin_opt *opt;
    struct opt_find_args args;
    char buf [256];
    char *name;

    if (strlcpy (buf, str, sizeof (buf)) >= sizeof (buf)) {
        error ("plugin option \"%s\" too big. Ignoring.", str);
        return (NULL);
    }

    if (!(name = strchr (buf, ':'))) {
        error ("Malformed plugin option \"%s\" recieved. Ignoring", str);
        return (NULL);
    }

    *(name++) = '\0';

    args.optname = buf;
    args.plugin_name = name;

    opt = list_find_first (option_cache, (ListFindF) _opt_find, &args);

    if (opt == NULL) {
        error ("warning: plugin \"%s\" option \"%s\" not found.", name, buf);
        return (NULL);
    }

    return (opt);
}

int
spank_get_remote_options (job_options_t opts)
{
    const struct job_option_info *j;

    job_options_iterator_reset (opts);
    while ((j = job_options_next (opts))) {
        struct spank_plugin_opt *opt;
        struct spank_option     *p;

        if (j->type != OPT_TYPE_SPANK)
            continue;

        if (!(opt = _find_remote_option_by_name (j->option)))
            continue;

        p = opt->opt;

        if (p->cb && (((*p->cb) (p->val, j->optarg, 1)) < 0)) {
            error ("spank: failed to process option %s=%s", 
                    p->name, j->optarg);
        }
    }

    return (0);
}



/*
 *  Global functions for SPANK plugins
 */
spank_err_t spank_get_item (spank_t spank, spank_item_t item, ...)
{
    int        *p2int;
    uint32_t   *p2uint32;
    uint16_t   *p2uint16;
    uid_t      *p2uid;
    gid_t      *p2gid;
    pid_t      *p2pid;
    char     ***p2argv;
    va_list     vargs;
    spank_err_t rc = ESPANK_SUCCESS;

    if ((spank == NULL) || (spank->magic != SPANK_MAGIC))
        return (ESPANK_BAD_ARG);

    va_start (vargs, item);
    switch (item) {
        case S_JOB_UID:
            p2uid = va_arg (vargs, uid_t *);
            *p2uid = spank->job->uid;
            break;
        case S_JOB_GID:
            p2gid = va_arg (vargs, gid_t *);
            *p2gid = spank->job->gid;
            break;
        case S_JOB_ID:
            p2uint32 = va_arg (vargs, uint32_t *);
            *p2uint32 = spank->job->jobid;
            break;
        case S_JOB_STEPID:
            p2uint32 = va_arg (vargs, uint32_t *);
            *p2uint32 = spank->job->stepid;
            break;
        case S_JOB_NNODES:
            p2uint32 = va_arg (vargs, uint32_t *);
            *p2uint32 = spank->job->nnodes;
            break;
        case S_JOB_NODEID:
            p2uint32 = va_arg (vargs, uint32_t *);
            *p2uint32 = spank->job->nodeid;
            break;
        case S_JOB_LOCAL_TASK_COUNT:
            p2uint32 = va_arg (vargs, uint32_t *);
            *p2uint32 = spank->job->ntasks;
            break;
        case S_JOB_TOTAL_TASK_COUNT:
            p2uint32 = va_arg (vargs, uint32_t *);
            *p2uint32 = spank->job->nprocs;
            break;
        case S_JOB_NCPUS:
            p2uint16 = va_arg (vargs, uint16_t *);
            *p2uint16 = spank->job->cpus;
            break;
        case S_JOB_ARGV:
            p2int = va_arg (vargs, int *);
            *p2int = spank->job->argc;
            p2argv = va_arg (vargs, char ***);
            *p2argv = spank->job->argv;
            break;
        case S_JOB_ENV:
            p2argv = va_arg (vargs, char ***);
            *p2argv = spank->job->env;
            break;
        case S_TASK_ID:
            p2int = va_arg (vargs, int *);
            if (!spank->task) {
                *p2int = -1;
                rc = ESPANK_NOT_TASK;
            } else {
                *p2int = spank->task->id;
            }
            break;
        case S_TASK_GLOBAL_ID:
            p2uint32 = va_arg (vargs, int *);
            if (!spank->task) {
                rc = ESPANK_NOT_TASK;
            } else {
                *p2uint32 = spank->task->gtid;
            }
            break;
        case S_TASK_EXIT_STATUS:
            p2int = va_arg (vargs, int *);
            if (!spank->task || !spank->task->exited) {
                rc = ESPANK_NOT_TASK;
            } else {
                *p2int = spank->task->estatus;
            }
            break;
        case S_TASK_PID:
            p2pid = va_arg (vargs, pid_t *);
            if (!spank->task) {
                rc = ESPANK_NOT_TASK;
                *p2pid = 0;
            } else {
                *p2pid = spank->task->pid;
            }
            break;
        default:
            rc = ESPANK_BAD_ARG;
            break;
    }
    va_end (vargs);
    return (rc);
}

spank_err_t spank_getenv (spank_t spank, const char *var, char *buf, int len)
{
    char *val;

    if ((spank == NULL) || (spank->magic != SPANK_MAGIC))
        return (ESPANK_BAD_ARG);

    if (spank->job == NULL)
        return (ESPANK_BAD_ARG);

    if (len < 0)
        return (ESPANK_BAD_ARG);

    if (!(val = getenvp (spank->job->env, var)))
        return (ESPANK_ENV_NOEXIST);

    if (strlcpy (buf, val, len) >= len)
        return (ESPANK_NOSPACE);

    return (ESPANK_SUCCESS);
}

spank_err_t spank_setenv (spank_t spank, const char *var, const char *val,
        int overwrite)
{
    if ((spank == NULL) || (spank->magic != SPANK_MAGIC))
        return (ESPANK_BAD_ARG);

    if (spank->job == NULL)
        return (ESPANK_BAD_ARG);

    if (getenvp (spank->job->env, var) && !overwrite)
        return (ESPANK_ENV_EXISTS);

    if (setenvf (&spank->job->env, var, "%s", val) < 0)
        return (ESPANK_ERROR);

    return (ESPANK_SUCCESS);
}

/*
 * vi: ts=4 sw=4 expandtab
 */
