INDIR="$1"
SUFFIX=".yml"
SUFFIX_OUT=".cwl"

INITDIR=$(pwd)
cd $INDIR
INDIR=$(pwd)
CWLFILES=($( find -name *$SUFFIX))
cd $INITDIR

for FILE in ${CWLFILES[@]}
do
    echo " > converting: $FILE"
    cwl-upgrader $INDIR/$FILE > "$INDIR/${FILE%$SUFFIX}$SUFFIX_OUT"
done