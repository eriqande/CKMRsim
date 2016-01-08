#sh

# put default values here
# VAR=default



function usage {
      echo Syntax:
      echo "  $(basename $0)  Ped.out Def.in"
      echo
      echo "Ped.out is the mendel output file with the genotypes. Def.in is the mendel
input file with the number of alleles."
}

if [ $# -eq 0 ]; then
    usage;
    exit 1;
fi;

while getopts ":h" opt; do
    case $opt in
	h    )
	    usage
	    exit 1
	    ;;
	#m    )  VAR=$OPTARG;
	#    ;;
	\?   )
	    usage
	    exit  1
	    ;;
    esac
done

shift $((OPTIND-1));


# uncomment to test for right number of required args
if [ $# -ne 2 ]; then
    usage;
    exit 1;
fi

PED=$1
DEF=$2



(
  awk -F"," 'NF==3 {printf("num-alleles,%d\n", $3);}' $DEF;
  sed 's/:xyz//g; s/ *//g; s/\//,/g;' $PED
) | awk -F"," '
   $1=="num-alleles" {A[nL++] = $2; next}
   $2==1 || $2 == 2 {
    loc = 0;
    for(i=7;i<=NF;i+=2) {
      loc++; x=$i; y=$(i+1);
      if(x<y) {a=x; b=y}
      else {a=y; b=x}
      geno = 2 + (a - 1) * (A[loc] + 2) - (a * (a + 1) / 2) + (b - a)
      printf("%d %d %d %d\n", $1, $2, loc, geno);
      }
    }
'
