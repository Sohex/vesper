/* How many reals SHTns requires a spatial field to hold, against the NUGP the
 * model allocates.
 *
 *   gcc -O2 -I<prefix>/include probe_shtns_nspat.c -o probe -L<prefix>/lib \
 *       -lshtns_omp -lfftw3_omp -lfftw3 -lm -fopenmp
 *   ./probe <NTRU> <NLAT> <NLON>
 *
 * Worldbuilding frame: a COMPUTE check on the Vesper climate model's transform
 * library. Nothing here is about the simulated planet.
 *
 * WHY. `shtns.h` documents the spatial argument of `spat_to_SH` and
 * `SH_to_spat` as "a double array of size shtns->nspat", and separately says
 * the input is NOT GUARANTEED TO BE PRESERVED -- the library uses it as its own
 * FFT scratch. `shtnsmod.f90` hands those calls a local `zg(NUGP)`, where NUGP
 * is NLAT*NLON. If nspat exceeds NUGP the library reads and writes past the end
 * of a stack array, which a zero-initialised build hides and a sNaN-initialised
 * one traps on.
 *
 * It prints, and does not judge: the grid this model runs is a parameter, and
 * the answer wanted is the number, per resolution.
 */
#include <stdio.h>
#include <stdlib.h>
#include <shtns.h>

int main(int argc, char **argv)
{
    if (argc != 4) {
        fprintf(stderr, "usage: %s <NTRU> <NLAT> <NLON>\n", argv[0]);
        return 2;
    }
    int ntru = atoi(argv[1]), nlat = atoi(argv[2]), nlon = atoi(argv[3]);

    shtns_verbose(0);
    shtns_use_threads(1);
    /* The model's own configuration: shtnsmod.f90 shtns_setup. */
    shtns_cfg s = shtns_create(ntru, ntru, 1, sht_orthonormal);
    shtns_set_grid(s, sht_quick_init | SHT_PHI_CONTIGUOUS, 0.0, nlat, nlon);

    long nugp = (long)nlat * nlon;
    printf("NTRU %d  NLAT %d  NLON %d\n", ntru, nlat, nlon);
    printf("  NUGP  (what the model allocates) : %ld reals, %ld bytes\n",
           nugp, nugp * 8);
    printf("  nspat (what SHTns requires)      : %u reals, %lu bytes\n",
           s->nspat, (unsigned long)s->nspat * 8);
    printf("  overrun                          : %ld reals, %ld bytes\n",
           (long)s->nspat - nugp, ((long)s->nspat - nugp) * 8);
    return 0;
}
