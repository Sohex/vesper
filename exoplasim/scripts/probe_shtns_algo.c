/* Which transform algorithm does SHTns actually select, per resolution?
 *
 *   gcc -O2 -I<prefix>/include probe_shtns_algo.c -o probe -L<prefix>/lib \
 *       -lshtns_omp -lfftw3_omp -lfftw3 -lm -fopenmp
 *   ./probe <NTRU> <NLAT> <NLON> <quick|timed|fly>
 *
 * Worldbuilding frame: a COMPUTE check on the Vesper climate model's transform
 * library. Nothing here is about the simulated planet.
 *
 * WHY. CLIM-74 was opened to measure what `SHT_QUICK_INIT` costs against
 * STORED LEGENDRE TABLES, on the reading that `nm` showing only `_fly` kernels
 * meant the library was recomputing what it could have streamed. The premise
 * does not survive reading the pinned revision: in SHTns 3.7.5 the matrix-based
 * algorithms are not merely unselected, they are ABSENT. `init_sht_array_func`
 * populates `SHT_ODD`, `SHT_SV` and the `SHT_FLY*`/`SHT_OMP*` slots and never
 * writes an `SHT_MEM` entry; `choose_best_sht` opens with
 * `const int on_the_fly_only = 1` and overrides `i0 = SHT_MEM` with `i0 =
 * SHT_SV`; `shtns_set_grid_auto` declares `const int on_the_fly = 1`.
 *
 * So the real choice is between ON-THE-FLY VARIANTS, which differ in unroll
 * depth (NWAY): `sht_quick_init` pins FLY2 by heuristic, and `sht_gauss` times
 * the available variants and keeps the winner. How many variants are AVAILABLE
 * is itself per-resolution -- `init_sht_array_func` caps `alg_lim` by
 * `nlat_2 / VSIZE2` -- which is why the pick is a per-rung question rather than
 * one answer for the ladder.
 *
 * This probe asks the library rather than the source. It prints, and does not
 * judge.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <shtns.h>

int main(int argc, char **argv)
{
    if (argc != 5) {
        fprintf(stderr, "usage: %s <NTRU> <NLAT> <NLON> <quick|timed|fly>\n", argv[0]);
        return 2;
    }
    int ntru = atoi(argv[1]), nlat = atoi(argv[2]), nlon = atoi(argv[3]);
    /* `fly` is SHT_GAUSS_FLY, which the task row called deterministic by
     * construction. It is not: shtns_set_grid_auto sets quick_init only for
     * quick_init, reg_fast and reg_poles, so gauss_fly falls through to
     * choose_best_sht exactly as sht_gauss does. This mode exists so that
     * claim can be tested rather than repeated. */
    int timed = (strcmp(argv[4], "timed") == 0);
    int fly = (strcmp(argv[4], "fly") == 0);

    /* Verbose 1 makes the library name the algorithm it settled on. That
     * printout IS the measurement here, so it goes to stdout as the library
     * writes it and is not reformatted. */
    shtns_verbose(1);
    shtns_use_threads(1);          /* the model's team is the parallelism */
    shtns_cfg s = shtns_create(ntru, ntru, 1, sht_orthonormal);
    enum shtns_type mode = fly ? sht_gauss_fly
                         : (timed ? sht_gauss : sht_quick_init);
    shtns_set_grid(s, mode | SHT_PHI_CONTIGUOUS, 0.0, nlat, nlon);
    printf("nlat %d nlon %d ntru %d mode %s\n", nlat, nlon, ntru, argv[4]);
    return 0;
}
