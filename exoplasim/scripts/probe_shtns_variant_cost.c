/* What does the on-the-fly VARIANT choice cost, per transform call, per rung?
 *
 *   gcc -O2 -I<prefix>/include probe_shtns_variant_cost.c -o probe \
 *       -L<prefix>/lib -lshtns_omp -lfftw3_omp -lfftw3 -lm -fopenmp
 *   ./probe <NTRU> <NLAT> <NLON> <reps>
 *
 * Worldbuilding frame: a COMPUTE measurement of the transform library the
 * Vesper climate model calls. Nothing here is about the simulated planet.
 *
 * WHY THIS AND NOT A MODEL BENCH. CLIM-74's question is which on-the-fly
 * variant to pin per resolution. `shtns-viability.md` measures SHTns's own
 * kernels at 0.00% of model runtime at T170 with the wrappers at 0.74%
 * together, so a difference BETWEEN variants is a fraction of that, and a
 * model-level A/B whose self-scatter floor is 5% cannot resolve it -- that is
 * `docs/src/practice/failure-modes.md` class 34, an instrument too blunt for
 * the effect. This probe measures where the effect actually lives, with enough
 * repeats to say whether it is there at all.
 *
 * HOW A VARIANT IS FORCED. There is no public API for it. `config_load` reads
 * the plain-text `shtns_cfg` from the working directory before any timing and
 * suppresses the tuning pass, resolving algorithm NAMES through `sht_func`. So
 * the caller writes the file. This probe reports what the library says it
 * loaded, rather than assuming the request took: an unforced arm silently
 * running the default would compare fly2 against itself.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <shtns.h>

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

static int cmp(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

int main(int argc, char **argv)
{
    if (argc != 5) {
        fprintf(stderr, "usage: %s <NTRU> <NLAT> <NLON> <reps>\n", argv[0]);
        return 2;
    }
    int ntru = atoi(argv[1]), nlat = atoi(argv[2]), nlon = atoi(argv[3]);
    int reps = atoi(argv[4]);

    shtns_verbose(1);
    shtns_use_threads(1);
    shtns_cfg s = shtns_create(ntru, ntru, 1, sht_orthonormal);
    /* QUICK_INIT + LOAD_SAVE: loads an authored shtns_cfg if one matches this
     * grid, and falls back to the deterministic fly2 heuristic if not. */
    shtns_set_grid(s, sht_quick_init | SHT_PHI_CONTIGUOUS | SHT_LOAD_SAVE_CFG,
                   0.0, nlat, nlon);
    shtns_print_cfg(s);

    long nlm = s->nlm, nspat = s->nspat;
    cplx *Qlm = (cplx *)shtns_malloc(nlm * sizeof(cplx));
    cplx *Slm = (cplx *)shtns_malloc(nlm * sizeof(cplx));
    double *Qh = (double *)shtns_malloc(nspat * sizeof(double));
    double *Sh = (double *)shtns_malloc(nspat * sizeof(double));
    double *Th = (double *)shtns_malloc(nspat * sizeof(double));
    for (long i = 0; i < nlm; i++) { Qlm[i] = 0.0; Slm[i] = 0.0; }
    for (long i = 0; i < nspat; i++) { Qh[i] = 0.0; Sh[i] = 0.0; Th[i] = 0.0; }
    Qlm[2] = 1.0; Slm[2] = 1.0;

    /* The four calls shtnsmod.f90 actually makes: scalar synthesis and
     * analysis, and the vector pair behind sh_dv2uv and sh_uv2dv. */
    const char *name[4] = {"syn SH_to_spat", "ana spat_to_SH",
                           "vsy SHsphtor_to_spat", "van spat_to_SHsphtor"};
    double *t = (double *)malloc(reps * sizeof(double));

    for (int k = 0; k < 4; k++) {
        for (int r = 0; r < reps; r++) {
            double t0 = now();
            switch (k) {
                case 0: SH_to_spat(s, Qlm, Qh); break;
                case 1: spat_to_SH(s, Qh, Qlm); break;
                case 2: SHsphtor_to_spat(s, Slm, Qlm, Sh, Th); break;
                case 3: spat_to_SHsphtor(s, Sh, Th, Slm, Qlm); break;
            }
            t[r] = now() - t0;
        }
        qsort(t, reps, sizeof(double), cmp);
        double med = t[reps / 2];
        double q1 = t[reps / 4], q3 = t[(3 * reps) / 4];
        printf("RESULT %s nlat=%d median_us=%.3f iqr_us=%.3f reps=%d\n",
               name[k], nlat, med * 1e6, (q3 - q1) * 1e6, reps);
    }
    return 0;
}
