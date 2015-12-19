/** [REDACTED] is crash-only software, so the ENSURE macros below call
 * abort rather than trying to recover.
 */
#pragma once

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

static inline void verbose_abort(const char *, const char *, const char *, const char *)  __attribute__ ((noreturn));
static inline void verbose_abort(const char *fl, const char *fn, const char *s, const char *extra) {
     write(2, fl, strlen(fl)); write(2, fn, strlen(fn)); write(2, s, strlen(s));
     if (extra != NULL) write(2, extra, strlen(extra));
     write(2, "\n", 1);
     abort(); _exit(-1);
}
#define STRINGIFY_PRIME(x) #x
#define STRINGIFY(x) STRINGIFY_PRIME(x)
#define ABORT(x) do { verbose_abort(__FILE__ ":" STRINGIFY(__LINE__) ": ", __func__, ": failed assertion: " x, NULL); } while (0)
#define ABORT_S(x,y) do { verbose_abort(__FILE__ ":" STRINGIFY(__LINE__) ": ", __func__, ": failed assertion: " x, y); } while (0)
#define ENSURE(x) do { errno = 0; if (!(x)) { if (errno) perror(__func__); ABORT(STRINGIFY(x)); } } while (0)
#define ENSURE_0(x) ENSURE(0 == (x))
#define ENSURE_NOT(v,x) ENSURE((v) != (x))
