#!/bin/bash

# The defaults
DO_MAKE=1  # measure the compilation time scaling
COMPILER_CONFIG="benchmarking"  # for the setup script it means "./compilers/../benchmarking/${COMPILER_CONFIG}.in"
N_PROC_LIST=""

usage() {
    echo "Usage:   $0 [options] [n_threads_1 [n_threads_2 ...]]"
    echo "         -h | --help : this message"
    echo "         -f | --fast : skip compilation test if possible"
    echo "         -c | --config file : use ./benchmarking/\${file}.in as the compiler config (default: ./benchmarking/benchmarking.in)"
    echo "default: $0 \$( seq number_of_logical_CPUs )"
    exit 1
}

j=1
skip=0
for i in "$@" ; do
    j=$(( $j + 1 ))
    [ $skip == 0 ] && case $i in
	"-f" | "--fast") DO_MAKE=0 ;;  # are we allowed to skip making and want to just run existing piernik
	"-h" | "--help") usage ;;
	"-c" | "--config") COMPILER_CONFIG=${*:$j:1} ; skip=1 ;;
	*) N_PROC_LIST="${N_PROC_LIST} $i"
    esac || skip=0
done

SCALE=${BIG:=1}

MEMG=$( LC_ALL=C free -m | awk '/Mem/ {print $2/1024.}' )
MEMM=$( LC_ALL=C free -m | awk '/Mem/ {print $2}' )
# alternatively (should give value closer to the amount of physical RAM installed):
#MEMG=0
#for mem in /sys/devices/system/memory/memory*; do
#    [[ "$(cat ${mem}/online)" == "1" ]] && MEMG=$(( MEMG+$((0x$(cat /sys/devices/system/memory/block_size_bytes)))))
#done
#MEMG=$( echo $MEMG | awk '{print $1 / 1024.^3 }' )

echo "## "$( cat /proc/cpuinfo  | grep "model name" | uniq )
echo "## Memory : $MEMG GB"

[ "$SCALE" != "1" ] && echo "# test domains are scaled by factor of $SCALE"

# create list of thread count to be tested
if [ ${#N_PROC_LIST} == 0 ] ; then
    N=$( awk 'BEGIN {c=0} /processor/ {if ($NF > c) c=$NF} END {print c+1}' /proc/cpuinfo )
    N_PROC_LIST=$( seq $N )
else
    N_PROC_LIST=$( echo $N_PROC_LIST | awk '{for (i=1; i<=NF; i++) printf("%d\n",1*$i); print ""}' | sort -n | uniq )
fi

for i in $N_PROC_LIST ; do
    if [ $i -le 0 ] ; then
	echo "N_PROC_LIST: ${N_PROC_LIST}"
	usage
    fi
done

TIMEFORMAT='real/user/sys/CPU%%: %1R %1U %1S %P%%'
MAKE=make

B_PROBLEM_LIST="sedov crtest maclaurin"
if [ $DO_MAKE == 0 ] ; then
    for p in $B_PROBLEM_LIST ; do
	[ -x runs/${p}_B_${p}/piernik ] || DO_MAKE=1
    done
fi

awkfor() {
    awk '{printf("%19s %7s %7s %7s %8s\n", $1, $2, $3, $4, $5)}'
}

# some cleanup
PROBLEM_LIST="maclaurin advection_test crtest jeans tearing sedov otvortex 2body"

if [ $DO_MAKE == 1 ] ; then
    touch Makefile
    # run correctly also when the clock is messed up
    rm -rf obj_B_*
    for i in $PROBLEM_LIST; do
	rm -rf runs/${i}_B_${i}
    done

    {
	# create object and run directories
	echo -n "Preparing objects                "
	SETUP_PARAMS="-c ../benchmarking/${COMPILER_CONFIG} -n --linkexe -d BENCHMARKING_HACK "
	( time for i in $PROBLEM_LIST; do
	./setup $i $SETUP_PARAMS -o "B_"$i > /dev/null
	done ) 2>&1 | awkfor

	get_n_problems() {
	    echo $1 $PROBLEM_LIST | awk '{for (i=0; i<$1; i++) printf("obj_B_%s ",$(2+i)); print ""}'
	}

	#
	# Benchmarking: make objects
	#

	OBJ_LIST=$( get_n_problems 1 )

	echo -n "Single-thread make object        "
	( time $MAKE $OBJ_LIST > /dev/null ) 2>&1 | awkfor
	$MAKE $OBJ_LIST CL=1 > /dev/null

	echo -n "Multi-thread make object         "
	( time $MAKE -j $OBJ_LIST > /dev/null ) 2>&1 | awkfor
	$MAKE $OBJ_LIST CL=1 > /dev/null

	OBJ_LIST=$( get_n_problems 2 )

	echo -n "Multi-thread make two objects    "
	( time $MAKE -j $OBJ_LIST > /dev/null ) 2>&1 | awkfor
	$MAKE $OBJ_LIST CL=1 > /dev/null

	OBJ_LIST=$( get_n_problems 4 )

	echo -n "Multi-thread make four objects   "
	( time $MAKE -j $OBJ_LIST > /dev/null ) 2>&1 | awkfor
	$MAKE $OBJ_LIST CL=1 > /dev/null

	OBJ_LIST=$( get_n_problems 8 )

	echo -n "Multi-thread make eight objects  "
	( time $MAKE -j $OBJ_LIST > /dev/null ) 2>&1 | awkfor
    }
    echo
else
    echo "## compilation skipped"
fi

#
# Trying to find best approximation of cuboid of given volume and integer edges.
# Call:
#     weak_cube number_of_work_units
# Tries to find a good decomposition of number_of_work_units into three factors.
#

weak_cube() {
i=$1

f1=`echo $i | awk '{eps=1+1e-10; x=$1; x3=(eps*x)**(1./3.); i1=1; for (i=x; i>int(x3); i--) if (x-i*int(x/i) == 0) i1=i; xx=x/i1; x2=(eps*xx)**(1./2.); i2=1; for (i=xx; i>int(x2); i--) if (xx-i*int(xx/i) == 0) i2=i; i3=xx/i2; print i1, i2, i3}'`

f2=`echo $i | awk '{eps=1+1e-10; x=$1; x3=(eps*x)**(1./3.); for (i=1; i<=int(x3); i++) if (x-i*int(x/i) == 0) i1=i; xx=x/i1; x2=(eps*xx)**(1./2.); for (i=1; i<=int(x2); i++) if (xx-i*int(xx/i) == 0) i2=i; i3=xx/i2; print i1, i2, i3}'`

if [ `factor $i | wc -w` -gt 4 ] ; then
    f3=`factor $i | awk '{if (NF <= 4) print $0; else if (NF==5) print $2*$3, $4, $5 ; else if (NF==6) print $2*$5, $3*$4, $6; else if (NF==7) print $2*$7, $3*$6, $4*$5; else {ii[0]=1; ii[1]=1; ii[2]=1; for (i=2; i<=NF; i++) ii[(i-2)%3] *= $i; print ii[0], ii[1], ii[2]} }'`
else
    f3=`factor $i | awk '{f1=1; f2=1; f3=1; if (NF>=2) f1=$2; if (NF>=3) f2=$3; if (NF==4) f3=$4; print f1, f2, f3}'`
fi

for f in "$f1" "$f2" "$f3" ; do
    echo $f
done | awk '{surf[NR]=$2*$3+$2*$1+$3*$1; f[NR]=$0} END {mn=4*'$i'; ff=-1; for (i=1;i<=NR;i++) if (surf[i]<mn) {mn=surf[i]; ff=i} print f[ff]}'
}

#
# Second approach to find best approximation of cuboid of given volume and integer edges.
# Call:
#     weak_cube2 number_of_work_units domain_size_1D block_size
# Tries to find a good decomposition of number_of_work_units*domain_size_1D^3 into three factors that allow for block_size^3 block decompositions.
# For some prime numbers of number_of_work_units choose a good enough approximation.
#

weak_cube2() {
[ $(( ( $2 / $3 ) * $3 - $2 )) == 0 ] || exit 1
echo $1 $2 $3 | awk 'function min(x,y) {return (x<y)?x:y} \
     function abs(x) {return (x>0)?x:-x} {\
     rsave = 10.; \
     isave = 0; \
     n = $1; \
     b = ($2/$3)^3; \
     nb = int((n*b)^(1./3.)); \
     for (i=1; i<=3*nb; i++) \
         for (j=1; j<=min(3*nb, n*b/i+1); j++) \
             for (k=1; k<=min(3*nb, (n*b)/(i*j)+1); k++) {\
                 ijk=i*j*k; \
		 r= ijk/(1.*n*b); \
		 if (abs(r-1.) < abs(rsave-1.)) { \
		    rsave = r; \
		    isave = 0; \
		 } \
		 if (r == rsave) {\
		    isave++;
		    ai[isave] = i; \
		    aj[isave] = j; \
		    ak[isave] = k; \
		 } \
	     } \
     if (isave > 0) \
        for (i=1; i<=isave; i++) \
	    print ai[i]*$3, aj[i]*$3, ak[i]*$3, ai[i]+aj[i]+ak[i], ai[i]*aj[i] + ai[i]*ak[i] + aj[i]*ak[i], rsave; \
     }' | sort -n -k 5 | head -n 1 | awk '{print $1, $2, $3}'
}

#
# Benchmarking: Piernik
#

for p in $B_PROBLEM_LIST ; do
    for t in flood strong weak ; do
	echo "## "$( date )
	echo "Benchmarking $p, $t scaling"
	(
	    RUNDIR=runs/${p}_B_${p}
	    cp  $( dirname $0 )"/problem.par."${p} ${RUNDIR}/problem.par
	    cd $RUNDIR
	    case $p in
		sedov)     echo "#Threads dWallClock1 dWallClock2 dWallClock3 dWallClock4 dWallClock5 dWallClock_Average";;
		crtest)    echo "#Threads MG_prepare MG_cycle Total_MG";;
		maclaurin) echo "#Threads MG_prepare MG_i-cycle MG_multipole MG_o-cycle Total_MG";;
	    esac
	    for i in $N_PROC_LIST ; do
                MPIRUN="mpirun"
                mpirun -np $i echo > /dev/null 2>&1 || MPIRUN="mpirun --use-hwthread-cpus"
                # OpenMPI refuses to run more jobs than CPU cores but with --use-hwthread-cpus its performance is poor on 2 threads
                max_mem=$( echo $MEMM $i | awk '{print int(0.90*$1*1024/$2)}' )
                rm *log 2> /dev/null
		if [ $t == flood ] ; then  # refresh the subdirectories for flood scaling runs
		    for j in $( seq $i ) ; do
			[ -e $j ] && rm -rf $j
			mkdir $j
			cp piernik problem.par $j
		    done
		fi

		case $p in
		    sedov)
		        BS=8  # BS=1 would perfectly work, but the awk search takes quite long then, O(BS^3)
			NX=$( echo 64 $SCALE | awk '{print int($1*$2)}') ;;
		    crtest)
		        BS=8
			NX=$( echo 32 $SCALE | awk '{print int($1*$2)}') ;;
		    maclaurin)
		        BS=32  # Has to match AMR_bsize
			if [ $t == flood ] ; then
			    NX=$( echo 64 $SCALE | awk '{print int($1*$2)}')
			else
			    case $t in
				weak)
				    NX=$( echo 64 $SCALE | awk '{print int($1*$2)}') ;;
				strong)
				    NX=$( echo 128 $SCALE | awk '{print int($1*$2)}') ;;
			    esac
			fi ;;
	        esac
		f=$( weak_cube2 $i $NX $BS )
		f1=$( echo $f | awk '{print $1}' )
		f2=$( echo $f | awk '{print $2}' )
		f3=$( echo $f | awk '{print $3}' )

		case $p in
		    sedov)
			if [ $t == flood ] ; then
			    for j in $( seq $i ) ; do
				cd $j
				rm *log 2> /dev/null
				./piernik -n '&BASE_DOMAIN n_d = 3*'$NX' /' > _stdout_ 2> _stderr_ &
				cd - > /dev/null
			    done
			    wait
			    sleep 1
			    for j in $( seq $i ) ; do
				( grep "dWallClock" $j/_stdout_ || echo "" ) | awk 'BEGIN {t=0; n=0; printf("%3d",'$i');} {if ($3 != 0) {printf("%7.2f ", $12); t+=$12; n++;} } END {printf("%7.3f\n", t/n)}'
			    done
			else
			    echo -n $i
			    case $t in
				weak)
				    # For cartesian decomposition there is degradation instead of improvement
				    # xm=$( echo $f1 $NX | awk '{print $1/(1.0 * $2)}')
				    # ym=$( echo $f2 $NX | awk '{print $1/(1.0 * $2)}')
				    # zm=$( echo $f3 $NX | awk '{print $1/(1.0 * $2)}')
				    # $MPIRUN -np $i ./piernik -n '&BASE_DOMAIN n_d = '$f1', '$f2', '$f3' xmin = -'$xm' xmax = '$xm' ymin = -'$ym' ymax = '$ym' zmin = -'$zm' zmax = '$zm' /' 2> /dev/null ;;
                                    $MPIRUN -np $i ./piernik -n '&BASE_DOMAIN n_d = '$(( $i * $NX ))', 2*'$NX' xmin = -'$(( $i * 1 ))' xmax = '$(( $i * 1 ))'/' 2> _stderr_weak_$i | tee weak_$i ;;
				strong)
				    $MPIRUN -np $i ./piernik -n '&BASE_DOMAIN n_d = 3*'$NX' /' 2> _stderr_strong_$i | tee strong_$i ;;
			    esac | grep "dWallClock" | awk 'BEGIN {t=0; n=0;} {if ($12 != 0.) {printf("%7.2f ", $12); t+=$12; n++;} } END {printf("%7.3f ", t/n)}'
			fi
			echo ;;
		    crtest)
			if [ $t == flood ] ; then
			    for j in $( seq $i ) ; do
				cd $j
				rm *log 2> /dev/null
				./piernik -n '&BASE_DOMAIN n_d = 3*'$NX' /' > _stdout_ 2> _stderr_ &
				cd - > /dev/null
			    done
			    wait
			    sleep 1
			    for j in $( seq $i ) ; do
				grep "C01cycles" $j/_stdout_ | awk '{if (NR==1) printf("%d %7.3f %7.3f ", '$i', $5, $8)}'
				awk '/Spent/ { printf("%s\n",$5) }' $j/*log
			    done
			else
			    echo -n $i
			    case $t in
				weak)
				    xm=$( echo $f1 $NX | awk '{print 512. * $1/(1.0 * $2)}')
				    ym=$( echo $f2 $NX | awk '{print 512. * $1/(1.0 * $2)}')
				    zm=$( echo $f3 $NX | awk '{print 512. * $1/(1.0 * $2)}')
				    $MPIRUN -np $i ./piernik -n '&BASE_DOMAIN n_d = '$f1', '$f2', '$f3' xmin = -'$xm' xmax = '$xm' ymin = -'$ym' ymax = '$ym' zmin = -'$zm' zmax = '$zm' /' 2> _stderr_weak_$i | tee weak_$i ;;
				strong)
				    $MPIRUN -np $i ./piernik -n '&BASE_DOMAIN n_d = 3*'$NX' /' 2> _stderr_strong_$i | tee strong_$i ;;
			    esac | grep "C01cycles" | awk '{if (NR==1) printf("%7.3f %7.3f ", $5, $8)}'
			    awk '/Spent/ { printf("%s ",$5) }' *log
			fi
			echo ;;
		    maclaurin)
			SKIP=0
			if [ $t == flood ] ; then
			    for j in $( seq $i ) ; do
				{
				    cd $j
				    rm *log 2> /dev/null
				    ./piernik -n '&BASE_DOMAIN n_d = 3*'$NX' / &MEMORY max_mem = '$max_mem'/' > _stdout_ 2> _stderr_ ;
			        } 2> /dev/null &
			    done
			    wait
			    sleep 1
			    [ $SKIP == 0 ] && for j in $( seq $i ) ; do
				grep cycles $j/_stdout_ | awk 'BEGIN {printf("%d", '$i');} {printf("%7.3f %7.3f ", $5, $8)}'
				grep -q cycles $j/_stdout_ || echo ""
				awk '/Spent/ { printf("%s\n",$5) }' $j/*log
			    done
			else
			    echo -n $i
			    case $t in
				weak)
				    xm=$( echo $f1 $NX | awk '{print 2. * $1/(1.0 * $2)}')
				    ym=$( echo $f2 $NX | awk '{print 2. * $1/(1.0 * $2)}')
				    zm=$( echo $f3 $NX | awk '{print 2. * $1/(1.0 * $2)}')
			            $MPIRUN -np $i ./piernik -n '&BASE_DOMAIN n_d = '$f1', '$f2', '$f3' xmin = -'$xm' xmax = '$xm' ymin = -'$ym' ymax = '$ym' zmin = -'$zm' zmax = '$zm' / &AMR bsize = 3*32 / &MEMORY max_mem = '$max_mem'/' 2> _stderr_weak_$i | tee weak_$i | grep cycles | awk '{printf("%7.3f %7.3f ", $5, $8)}' ;;
				strong)
				    $MPIRUN -np $i ./piernik -n '&BASE_DOMAIN n_d = 3*'$NX' / &AMR bsize = 3*32 / &MEMORY max_mem = '$max_mem'/' 2> _stderr_strong_$i | tee strong_$i | grep cycles | awk '{printf("%7.3f %7.3f ", $5, $8)}' ;;
			    esac
			    [ $SKIP == 0 ] && awk '/Spent/ { printf("%s ", $5) }' *log
			fi
			echo ;;
	        esac
	    done
	) | awk '{ if ($0 ~ "##") {print "## ", $0} else { if (substr($1,0,1) == "#") split($0,form); if (NF>0) {for (i=1;i<=NF;i++) printf("%-*s ",length(form[i]),$i); print ""; fflush()} } }'
	echo
    done
done
