/* Who calls memcpy, memmove and memset, how often, and with how many bytes.
 *
 *   gcc -shared -fPIC -O2 -o count_memops.so count_memops.c -ldl
 *   LD_PRELOAD=./count_memops.so ./model
 *
 * Worldbuilding frame: a COMPUTE measurement of the Vesper climate model on
 * this desktop. Nothing here is about the simulated planet.
 *
 * WHY AN INTERPOSER RATHER THAN A PROFILE. perf says a third of the run at T170
 * is inside a handful of addresses in libc and cannot say who called them.
 * Resolved against libc's debuginfo, those addresses are
 * __memset_avx512_unaligned_erms -- the model spends a third of its time
 * ZEROING memory, not copying it, which is why memset is interposed here and
 * why the file is no longer named for memcpy alone. libc here has no symbols
 * for the IFUNC-resolved implementations, AMD IBS needs system-wide sampling
 * that perf_event_paranoid forbids, and the notes record that neither dwarf nor
 * frame pointers unwind out of libc on this host. Counting at the call site
 * sidesteps all of that: the return address IS the caller, exactly, with no
 * unwinding and no sampling error.
 *
 * WHAT IT CANNOT SEE, and this bounds the conclusion. Only calls that go
 * through the PLT are interposed. gfortran emits inline move loops for many
 * array copies and calls memcpy only above a size threshold, so a small count
 * here does NOT mean there is little copying -- it means little copying goes
 * through this symbol. Read a large count as evidence; read a small one as
 * "look somewhere else", not as absolution.
 *
 * The counters are per-thread to avoid both a lock in the hot path and the
 * false sharing that a shared array would give. Buckets are keyed on the return
 * address, resolved afterwards with addr2line against the binary.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define NSLOT 4096u

struct site {
    void    *ret;
    uint64_t calls;
    uint64_t bytes;
};

/* Threadprivate so the hot path takes no lock. */
static __thread struct site tbl[NSLOT];
static __thread int registered;

/* One list of every thread's table, so the dump can total them. */
static struct site *tables[512];
static int ntables;
static __thread int myslot;

static void *(*real_memcpy)(void *, const void *, size_t);
static void *(*real_memmove)(void *, const void *, size_t);
static void *(*real_memset)(void *, int, size_t);

/* A byte loop, for the window before the real one is known.
 *
 * THE WINDOW IS REAL AND IT SEGFAULTS. The dynamic loader calls memcpy while it
 * is still wiring this library up, which is BEFORE the constructor runs, so a
 * bare `real_memcpy(...)` dereferences NULL. It is also order dependent: the
 * same binary interposed the same way ran sixty steps and then crashed on five.
 * dlsym cannot be called from inside the window either, because dlsym itself
 * copies, which would recurse. */
static void *slow_copy(void *d, const void *s, size_t n)
{
    unsigned char *dd = (unsigned char *)d;
    const unsigned char *ss = (const unsigned char *)s;
    if (dd < ss) { for (size_t i = 0; i < n; i++) dd[i] = ss[i]; }
    else { for (size_t i = n; i-- > 0; ) dd[i] = ss[i]; }
    return d;
}

__attribute__((constructor))
static void init(void)
{
    real_memcpy  = dlsym(RTLD_NEXT, "memcpy");
    real_memmove = dlsym(RTLD_NEXT, "memmove");
    real_memset  = dlsym(RTLD_NEXT, "memset");
}

static inline void note(void *ret, size_t n)
{
    if (!registered) {
        registered = 1;
        myslot = __atomic_fetch_add(&ntables, 1, __ATOMIC_RELAXED);
        if (myslot < 512) tables[myslot] = tbl;
    }
    /* Knuth multiplicative hash on the return address. Collisions merge two
     * call sites, which shows up as an address that addr2line resolves oddly;
     * with a few dozen live sites in 4096 slots that is unlikely. */
    unsigned h = (unsigned)(((uintptr_t)ret * 2654435761u) >> 11) & (NSLOT - 1u);
    for (unsigned i = 0; i < 64; i++) {
        unsigned s = (h + i) & (NSLOT - 1u);
        if (tbl[s].ret == ret || tbl[s].ret == NULL) {
            tbl[s].ret = ret;
            tbl[s].calls++;
            tbl[s].bytes += n;
            return;
        }
    }
}

void *memcpy(void *d, const void *s, size_t n)
{
    if (__builtin_expect(real_memcpy == NULL, 0)) return slow_copy(d, s, n);
    note(__builtin_return_address(0), n);
    return real_memcpy(d, s, n);
}

void *memmove(void *d, const void *s, size_t n)
{
    if (__builtin_expect(real_memmove == NULL, 0)) return slow_copy(d, s, n);
    note(__builtin_return_address(0), n);
    return real_memmove(d, s, n);
}

void *memset(void *d, int c, size_t n)
{
    if (__builtin_expect(real_memset == NULL, 0)) {
        unsigned char *dd = (unsigned char *)d;
        for (size_t i = 0; i < n; i++) dd[i] = (unsigned char)c;
        return d;
    }
    note(__builtin_return_address(0), n);
    return real_memset(d, c, n);
}

__attribute__((destructor))
static void dump(void)
{
    const char *path = getenv("MEMCPY_COUNT_OUT");
    FILE *f = path ? fopen(path, "w") : stderr;
    if (!f) return;

    /* Merge every thread's table on the return address. */
    struct site all[NSLOT];
    memset(all, 0, sizeof all);
    unsigned used = 0;
    int n = ntables < 512 ? ntables : 512;
    for (int t = 0; t < n; t++) {
        if (!tables[t]) continue;
        for (unsigned i = 0; i < NSLOT; i++) {
            if (!tables[t][i].ret) continue;
            unsigned j;
            for (j = 0; j < used; j++)
                if (all[j].ret == tables[t][i].ret) break;
            if (j == used) { all[used].ret = tables[t][i].ret; used++; }
            all[j].calls += tables[t][i].calls;
            all[j].bytes += tables[t][i].bytes;
        }
    }
    /* Heaviest by bytes first: the question is traffic, not call count. */
    for (unsigned a = 0; a + 1 < used; a++)
        for (unsigned b = a + 1; b < used; b++)
            if (all[b].bytes > all[a].bytes) {
                struct site tmp = all[a]; all[a] = all[b]; all[b] = tmp;
            }
    /* FILE OFFSETS, not runtime addresses. The model is a PIE, so a return
     * address means nothing to addr2line without the load base subtracted.
     * dladdr gives the base of whichever object the address fell in, which also
     * separates the model's own call sites from a library's. Resolved here at
     * exit rather than in note(), because dladdr walks the link map and has no
     * business in the hot path. */
    fprintf(f, "# threads %d, distinct call sites %u\n", n, used);
    fprintf(f, "# file_offset,calls,bytes,object,symbol\n");
    for (unsigned i = 0; i < used; i++) {
        Dl_info info;
        const char *obj = "?";
        const char *sym = "?";
        uintptr_t off = (uintptr_t)all[i].ret;
        if (dladdr(all[i].ret, &info) && info.dli_fbase) {
            off = (uintptr_t)all[i].ret - (uintptr_t)info.dli_fbase;
            if (info.dli_fname) obj = info.dli_fname;
            if (info.dli_sname) sym = info.dli_sname;
        }
        fprintf(f, "0x%lx,%llu,%llu,%s,%s\n", (unsigned long)off,
                (unsigned long long)all[i].calls,
                (unsigned long long)all[i].bytes, obj, sym);
    }
    if (path) fclose(f);
}
