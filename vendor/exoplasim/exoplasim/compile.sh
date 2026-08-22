#!/bin/bash
#
#               EXOPLASIM COMPILATION SCRIPT
#
#  Usage:
#       ./compile.sh -p 8 -r T21 -v 10 -n 16 -O mavx 
#
#
#  Options:
#
#      -p:   Set floating point precision. Argument should be either 4 or 8. Default: 4
#
#      -r:   Set resolution. Can be T21/T42/T63, or alternatively the number of latitudes.
#            Default: T21
#
#      -v:   Set number of vertical levels. Default: 10
#
#      -n:   Number of processing cores (MPI threads) to use. Default: 4
#
#      -t:   Number of years to use in most_plasim_run script. Default: 10
#
#      -O:   Specify additional compiler optimization flags, like -mavx (leave out prepended
#            hyphen).
#
#      -d:   Compile in debug mode (will produce line-number tracebacks on crash)
#
helptext=$(cat <<-END
               EXOPLASIM COMPILATION SCRIPT

  Usage:
       ./compile.sh -p 8 -r T21 -v 10 -n 16 -O mavx 


  Options:

      -p:   Set floating point precision. Argument should be either 4 or 8. Default: 4

      -r:   Set resolution. Can be T21/T42/T63, or alternatively the number of latitudes.
            Default: T21

      -v:   Set number of vertical levels. Default: 10

      -n:   Number of processing cores (MPI threads) to use. Default: 4

      -t:   Number of years to use in most_plasim_run script. Default: 10
      
      -f:   Force FFT991mod.f90

      -O:   Specify additional compiler optimization flags, like -mavx (leave out prepended
            hyphen).

      -d:   Compile in debug mode (will produce line-number tracebacks on crash)

      -m:   Compile with Mars routines
      
      -h:   Output this text
END
)


prec=4

# SHTns, from the project-local prefix exoplasim/scripts/build_shtns.sh builds.
# The include path goes into the compile flags because shtnsmod includes
# shtns.f03; the library goes on the link line. Both are empty and harmless if
# it has not been built, and shtnsmod is the only object that needs either --
# but the model links it unconditionally, so a run with nshtns=0 is unaffected
# and a run with nshtns=1 fails at build rather than at midnight.
# Three levels up: this script is at vendor/exoplasim/exoplasim/. It MUST be
# resolved here, at the top, because compile.sh cds into plasim/bld before
# it builds -- resolving it down there lands one directory short and the
# only symptom is a missing include much later.
SHTNS_PREFIX="$(cd "$(dirname "$0")/../../.." && pwd)/vendor/shtns-install"
if [ -f "$SHTNS_PREFIX/include/shtns.f03" ]; then
    export SHTNS_INC="-I$SHTNS_PREFIX/include"
    # -lgomp explicitly: libshtns_omp pulls fftw3_omp, which needs it, and the
    # MPI build does not compile with -fopenmp so it would not otherwise be on
    # the link line. The MPI build cannot USE SHTns -- it has no single address
    # space holding every latitude -- but it still links shtnsmod, so it still
    # has to resolve the symbols.
    export SHTNS_LIB="$(ls "$SHTNS_PREFIX"/lib/libshtns*.a | head -1) -lfftw3_omp -lfftw3 -lgomp -lm"
else
    echo "NOTE: no SHTns at $SHTNS_PREFIX; building without it."
    echo "      exoplasim/scripts/build_shtns.sh installs it. nshtns=1 needs it."
    export SHTNS_INC=""
    export SHTNS_LIB=""
fi

resolution="t21"
latitudes=32
longitudes=64
fftopt="fftmod"
levels=10
ncpus=4
parmode=""
debug=0
optimization=""
nopt=0
years=1
nmars=0

while getopts "p:r:v:n:O:t:jdhmu" opt; do
    case $opt in
        p)
            case $OPTARG in
                4)
                    prec=4 ;;
                8)
                    prec=8 ;;
                single)
                    prec=4 ;;
                double)
                    prec=8 ;;
                \?)
                    echo "INVALID PRECISION. Reverting to single precision."
                    ;;
            esac
            ;;
        r)
            case $OPTARG in #for >T63, recommend 30 levels
                T21)
                    resolution="t21"
                    latitudes=32
                    longitudes=64
                    ;;
                T42)
                    resolution="t42"
                    latitudes=64
                    longitudes=128
                    ;;
                T63)
                    resolution="t63"
                    latitudes=96
                    longitudes=192
                    fftopt="fft991mod"
                    ;;
                T85)
                    resolution="t85"
                    latitudes=128
                    longitudes=256
                    ;;
                T106)
                    resolution="t106"
                    latitudes=160
                    longitudes=320
                    fftopt="fft991mod"
                    ;;
                T127)
                    resolution="t127"
                    latitudes=192
                    longitudes=384
                    ;;
                T170)
                    resolution="t170"
                    latitudes=256
                    longitudes=512
                    ;;
                32)
                    resolution="t21"
                    latitudes=32
                    longitudes=64
                    ;;
                64)
                    resolution="t42"
                    latitudes=64
                    longitudes=128
                    ;;
                96)
                    resolution="t63"
                    latitudes=96
                    longitudes=192
                    fftopt="fft991mod"
                    ;;
                128)
                    resolution="t85"
                    latitudes=128
                    longitudes=256
                    ;;
                160)
                    resolution="t106"
                    latitudes=160
                    longitudes=320
                    fftopt="fft991mod"
                    ;;
                192)
                    resolution="t127"
                    latitudes=192
                    longitudes=384
                    ;;
                256)
                    resolution="t170"
                    latitudes=256
                    longitudes=512
                    ;;
                \?)
                    resolution="t21"
                    latitudes=32
                    longitudes=64
                    echo "INVALID RESOLUTION PASSED! Reverting to T21."
                    ;;
            esac ;;
        f)    
            fftopt="fft991mod"
            ;;
        v)
            levels=$OPTARG
            ;;
        n)
            ncpus=$OPTARG
            ;;
        j)
            # Threads instead of ranks: one process, NPRO OpenMP threads,
            # mpimod_omp in place of mpimod. The executable is named apart
            # from the MPI one because it is a different binary identity at
            # the same resolution and rank count, not a variant of it.
            parmode="omp"
            ;;
        u)
            # Unpaired latitudes: the stock contiguous layout. Required for
            # SHTns, which needs the grid in latitude order, and pointless
            # against it -- the paired layout exists to fold a mirror pair
            # together inside legmod, which SHTns replaces. A separate binary
            # identity because the grid layout genuinely differs.
            nopairlat=1
            ;;
        d)
            debug=1
            ;;
        t)
            years=$OPTARG
            ;;
        O)
            optimization="-"$OPTARG
            nopt=1
            ;;
        h)
            echo "$helptext"
            exit 0
            ;;
        m)
            echo "Compiling for Mars..."
            export PLAMOD=p_mars
            nmars=1
            ;;
        \?)
            echo "UNRECOGNIZED ARGUMENT PASSED: "$OPTARG
            ;;
    esac
done

suffix=""
[ "$parmode" = "omp" ] && suffix="_omp"
[ "${nopairlat:-0}" = 1 ] && suffix="${suffix}_np"
echo "PRODUCING: "$optimization" -r"$prec" -o most_plasim_"$resolution"_l"$levels"_p"$ncpus$suffix".x"
executable="most_plasim_"$resolution"_l"$levels"_p"$ncpus$suffix".x"

echo "Writing resmod.f90....."

# `plasim/bld` is a build directory and is not tracked, so on a fresh clone
# it does not exist. Without the guard below the `cd` fails, execution
# continues because this script does not set -e, and the `rm -rf *` on the
# next line then runs in the PACKAGE ROOT and deletes the source tree.
mkdir -p plasim/bld plasim/bin
cd plasim/bld/ || { echo "compile.sh: cannot enter plasim/bld" >&2; exit 1; }
rm -rf *


#       ! T85L30 on 16/32/64 processors
#       !parameter(NLAT_ATM = 128)
#       !parameter(NLEV_ATM = 30)
#       !!parameter(NPRO_ATM = 16)
#       !parameter(NPRO_ATM = 32)
#       !!parameter(NPRO_ATM = 64)
# 
# 
#       ! T127L30 on 16/32/64 processors
#       !parameter(NLAT_ATM = 192)
#       !parameter(NLEV_ATM = 30)
#       !!parameter(NPRO_ATM = 32)
#       !parameter(NPRO_ATM = 48) ! Does not work so well, 32 more efficient? \
# 
# 
# 
# 
#       ! T170L30 on 16/32/64 processors
#       !parameter(NLAT_ATM = 256)
#       !parameter(NLEV_ATM = 30)
#       !parameter(NPRO_ATM = 32)
#       !parameter(NPRO_ATM = 64) 


echo "      module resmod ! generated by compile.sh ">resmod.f90
echo "      parameter(NLAT_ATM = "$latitudes") ">>resmod.f90
echo "      parameter(NLEV_ATM = "$levels") ">>resmod.f90
echo "      parameter(NPRO_ATM = "$ncpus") ">>resmod.f90
echo "      end module resmod ">>resmod.f90
echo " ">>resmod.f90

rm plasim.x
rm ../bin/$executable
rm ../run/$executable
cp -p ../src/* .

# The build directory carries a marker for which parallel layer it was last
# built with, because the three share object and .mod names and a stale one
# links silently. Switching layers empties it.
if [ "$parmode" = "omp" ]
then
    [ ! -e OMP ] && rm -f *.o *.mod *.x MPI OMP
    touch OMP
    cp ../../most_compiler_omp compilerargs
elif [ "$ncpus" -gt 1 ]
then
    [ ! -e MPI ] && rm -f *.o *.mod *.x MPI OMP
    touch MPI
    cp ../../most_compiler_mpi compilerargs
else
    [ ! -e MPI ] && rm -f *.o *.mod *.x MPI OMP
    cp ../../most_compiler compilerargs
fi

# Precision gets a marker for the same reason, and it is the worst of the three
# to get wrong: -fdefault-real-8 changes the width of `real` in every
# declaration and every interface, so a real*4 object linked against real*8 ones
# does not fail to link -- it computes. The mismatch reaches you as numbers that
# are merely wrong, and it cost most of a session being read as a transform bug.
if [ ! -e "PREC$prec" ]
then
    rm -f *.o *.mod *.x PREC4 PREC8
fi
touch "PREC$prec"

# The pairing switch gets its own marker for the same reason the parallel layer
# does: it changes a parameter nearly every object sees, and the two settings
# share object and .mod names, so a stale one links silently.
if [ "${nopairlat:-0}" = 1 ]
then
    [ ! -e NOPAIR ] && rm -f *.o *.mod *.x
    touch NOPAIR
    sed -i.bak '3s/$/ -DNOPAIRLAT/' compilerargs && rm -f compilerargs.bak
else
    [ -e NOPAIR ] && rm -f *.o *.mod *.x NOPAIR
fi

dbgs=""
(($debug)) && dbgs='../../most_debug_options '

cp ../../most_precision_optionsx precisionargsx
if [ "$prec" -gt 4 ]
then
    sed -i.bak '1s/$/'$prec'/' precisionargsx
    rm -rf precisionargsx.bak
else
    echo "">precisionargsx
fi

(($nopt)) && sed -i.bak '3s/$/ '$optimization'/' compilerargs && rm -rf compilerargs.bak

cat compilerargs $dbgs../bld/precisionargsx make_plasim > makefile

echo "Writing makefile..."
echo ""

export OCEANCOUP=cpl_stub

export FFTMOD=$fftopt
#cat makefile

make -e
./most_snow_build$prec
./most_ice_build$prec
cp plasim.x ../bin/$executable
cp ../bin/$executable ../run/
cd ../../

(($nmars)) && cp plasim/dat/T"${resolution[@]:1}"_mars/* plasim/run/ || cp plasim/dat/T"${resolution[@]:1}"/* plasim/run/

namelist=example.nl
(($nmars)) && namelist=mars.nl
(($nmars)) && cp postprocessor/mars.nl plasim/run/

snamelist=snapshot.nl
(($nmars)) && snamelist=mars_snapshot.nl
(($nmars)) && cp postprocessor/mars_snapshot.nl plasim/run/

echo "#!/bin/bash ">plasim/run/most_plasim_run
echo "# run-script generated by compile.sh                       ">>plasim/run/most_plasim_run
echo "EXP=MOST    # Name your experiment here                    ">>plasim/run/most_plasim_run
echo "[ \$# ==1 ] && cd \$1                                      ">>plasim/run/most_plasim_run
echo "rm -f plasim_restart                                       ">>plasim/run/most_plasim_run
echo "rm -f Abort_Message                                        ">>plasim/run/most_plasim_run
echo "YEAR=0                                                     ">>plasim/run/most_plasim_run
echo "YEARS="$years"                                             ">>plasim/run/most_plasim_run
echo "while [ \$YEAR -lt \$YEARS ]                               ">>plasim/run/most_plasim_run
echo "do                                                         ">>plasim/run/most_plasim_run
echo "   YEAR=\`expr \$YEAR + 1\`                                ">>plasim/run/most_plasim_run
echo "   DATANAME=\`printf '%s.%05d' \$EXP \$YEAR\`              ">>plasim/run/most_plasim_run
echo "   SNAPNAME=\`printf '%s_SNAP.%05d' \$EXP \$YEAR\`         ">>plasim/run/most_plasim_run
echo "   DIAGNAME=\`printf '%s_DIAG.%05d' \$EXP \$YEAR\`         ">>plasim/run/most_plasim_run
echo "   RESTNAME=\`printf '%s_REST.%05d' \$EXP \$YEAR\`         ">>plasim/run/most_plasim_run
echo "   SNOWNAME=\`printf '%s_SNOW.%05d' \$EXP \$YEAR\`         ">>plasim/run/most_plasim_run
if [ "$parmode" = "omp" ]
then
   # One process. The thread count is compiled in, so the launcher supplies
   # only the two settings the team cannot run correctly without.
   #
   # OMP_STACKSIZE, because -fopenmp implies -frecursive and the model's large
   # local arrays become stack-allocated; 8 MB is not enough for them.
   #
   # OMP_PROC_BIND and OMP_PLACES, because libgomp's DEFAULT IS UNBOUND. Left
   # alone, threads land on arbitrary CPUs, the placement changes run to run,
   # and some threads share a physical core while whole cores sit idle -- which
   # on this processor also randomises which die a thread lands on, and the two
   # dies are not interchangeable at T127 and above. Bound this way a thread
   # takes core t, which is what Open MPI gives rank t, so the two layers are
   # measured on the same placement rather than on their defaults.
   # ulimit sizes the MASTER thread; OMP_STACKSIZE sizes the others. Both are
   # needed: at T127 the master alone overruns a 16 MB limit and segfaults.
   echo "   ulimit -s unlimited 2>/dev/null || ulimit -s 1048576         ">>plasim/run/most_plasim_run
   echo "   export OMP_STACKSIZE=\${OMP_STACKSIZE:-512M}                ">>plasim/run/most_plasim_run
   echo "   export OMP_PROC_BIND=\${OMP_PROC_BIND:-close}               ">>plasim/run/most_plasim_run
   echo "   export OMP_PLACES=\${OMP_PLACES:-cores}                     ">>plasim/run/most_plasim_run
   echo "   ./$executable                                        ">>plasim/run/most_plasim_run
elif [ "$ncpus" -gt 1 ]
then
   MPI_RUN=$(head -n 1 most_compiler_mpi | tr "=" "\n" | tail -1)
   echo "   $MPI_RUN -np $ncpus $executable                      ">>plasim/run/most_plasim_run
else
   echo "   ./$executable                                        ">>plasim/run/most_plasim_run
fi
echo "   [ -e Abort_Message ] && exit 1                          ">>plasim/run/most_plasim_run
echo "   [ -e plasim_output ] && mv plasim_output \$DATANAME     ">>plasim/run/most_plasim_run
echo "   [ -e plasim_snapshot ] && mv plasim_snapshot \$SNAPNAME ">>plasim/run/most_plasim_run
echo "   [[ -e burn7.x && -e "$namelist" && -e \$DATANAME ]] && ./burn7.x -n <"$namelist">burnout \$DATANAME \$DATANAME.nc ">>plasim/run/most_plasim_run
echo "   [[ -e burn7.x && -e "$snamelist" && -e \$SNAPNAME ]] && ./burn7.x -n <"$snamelist">snapout \$SNAPNAME \$SNAPNAME.nc ">>plasim/run/most_plasim_run
echo "   [ -e plasim_diag ] && mv plasim_diag \$DIAGNAME         ">>plasim/run/most_plasim_run
echo "   [ -e plasim_status ] && cp plasim_status plasim_restart ">>plasim/run/most_plasim_run
echo "   [ -e plasim_status ] && mv plasim_status \$RESTNAME     ">>plasim/run/most_plasim_run
echo "   [ -e restart_snow ] && mv restart_snow \$SNOWNAME       ">>plasim/run/most_plasim_run
echo "done                                                       ">>plasim/run/most_plasim_run

