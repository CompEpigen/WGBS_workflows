#!/bin/bash
if [ $# -ne 2 ]; then
	echo "This shell will split a large read file, *.fastq, *.csfasta or *.fasta *.qual into many small reads file."
	echo "Each file should have about one million reads."
	echo "The 1st parameter is the input file name; The 2nd parameter is the output file prefix."
	echo "A list of split files <prefix>.txt will be generated." 
	echo "You should be able to use PerM Ref.fa <prefix>.txt <Options> to do the mapping."
	echo "Note this shell split reads only by lines. Additional comments line may cuase reads shifts"
	echo "For SOLiD reads, the shell will also try to split the corresponding *.qual file"
	echo "Be very careful about the file format to make sure the quality scores and reads are matched."
fi
iFileN=$1
oFileNPrefix=$2

if test ! -s $iFileN; then
	echo "File $iFileN not found or empty"
	exit 1
fi

function countCommentLine {	
	echo $#
	if [ $# -gt 0 ]; then
		iFileN=$1
		echo $iFileN
		if [[ $# == 1 ]]; then
			echo "Lines starts with # are treated as comments"
			commentSymbol='#'
		else 
			commentSymbol=$2
		fi
	else
		echo "The 1st arg is the iFileN"
		echo "The 2nd arg is the comment symbol"
		exit 1
	fi
	
    let commentLineNo=0
	exec < $iFileN
	while read line
	do
		commentLine=$(echo $line | grep -P "^$commentSymbol")
		#echo $commentLine
		if [[ $commentLine != '' ]]; then
			let commentLineNo=commentLineNo+1
		else
			if [[ $commentLineNo -gt 0 ]]; then
				echo "The first "$commentLineNo" lines seems to be comment lines."
				echo "Try to use sed '1,$commentLineNo""d' <reads.csfasta> to delete the first $commentLineNo lines."
			fi
			return $commentLineNo 
		fi		
	done
	#echo $commentLineNo
	return 	$commentLineNo
}

function checkExtraLines {
	if [ $# -ne 2 ]; then
		echo "This shell coutn how many lines starts with the comment symbol in the beginning of the file."
		echo "The 1st arg are iFileN."
		echo "The 2nd symbol are."
	fi
 	iFileN=$1
	linesPerReads=$2
	# check head 
	firstLine=$(head -n 1 $iFileN)
	rightLine=$(echo $firstLine | grep -P "^[>@]") 
	if [[  $rightLine == '' ]]; then
		echo "Additional or missing lines in the head of $iFileN?"
		echo 'with extra line:' $firstLine
		countCommentLine $iFileN 
		return 0
	fi 
	# check tail
	lastLine=$(tail -n $linesPerReads $iFileN) 
	rightLine=$(echo $lastLine | head -n 1 | grep -P "^[>@]")
	if [[ $rightLine == '' ]]; then
		echo "Additional or missing lines  in the tail of $iFileN?"
		echo 'with extra lines:' 
		tail -n $linesPerReads $iFileN 
		return 0
	fi
	# check tail
}

readsSetFileN=$(basename $iFileN)
baseName=$(basename $readsSetFileN ".csfasta")
baseName=$(basename $baseName ".fasta")
baseName=$(basename $baseName ".fa")
baseName=$(basename $baseName ".fastq")
baseName=$(basename $baseName ".fq")
extName=$(echo $readsSetFileN | sed "s/$baseName//")
readsPerFile=1000000
if [[ $extName == '.fastq' || $extName == '.fq' ]];then
	linesPerRead=4
elif [[ $extName == '.fasta' || $extName == '.fa' || $extName == '.csfasta' ]];then
	linesPerRead=2
else
	echo "Non-recognized ext name" 
	exit 1
fi


iFileLineNo=$(wc -l $iFileN | cut -f 1 -d ' ' ) 
let MODULUS=$iFileLineNo%$linesPerRead
if [ $MODULUS -ne 0 ]; then
	checkExtraLines $iFileN $linesPerRead
	echo "Remove extra lines before split."
	exit 1
#else 
	#checkExtraLines $iFileN
fi

linesPerFile=$(expr $readsPerFile \* $linesPerRead)
#echo $iFileLineNo '-le' $linesPerFile
if [ $iFileLineNo -le $linesPerFile ]; then
   echo "The input file has less than $readsPerFile reads."
   echo "No need to split the file."
   echo "Use PerM <Ref.fa> $iFileN [options] dirrectly."
   exit 1
fi 

# split the file # Should check the lines
split -l $linesPerFile $iFileN $oFileNPrefix'_' 
splitFiles=$(ls | grep -P ^$oFileNPrefix'_'[a-z]+[a-z]$)
for splitFile in $splitFiles;
do mv $splitFile $splitFile$extName; done

# make a split reads file list
listFile=$oFileNPrefix".txt"
ls | grep -P ^$oFileNPrefix'_'[a-z]+[a-z]$extName$ > $listFile
noReadsFiles=$(ls | grep -P ^$oFileNPrefix'_'[a-z]+[a-z]$extName$ | wc -l)
echo "Sucessfully split the reads into $noReadsFiles files"
echo "and make the read list $listFile."

# split the file name
qualFileN=''
if test -s $baseName'.qual'; then
	qualExtFileN='.qual'
	qualFileN=$baseName'.qual'
elif test -s $baseName'.QUAL'; then
	qualExtFileN='.QUAL'
	qualFileN=$baseName'.QUAL'
elif test -s $baseName'_QV.QUAL'; then
	qualExtFileN='_QV.QUAL'
	qualFileN=$baseName'_QV.QUAL'
elif test -s $baseName'_QV.qual'; then
	qualExtFileN='_QV.qual'
	qualFileN=$baseName'_QV.qual'
fi

if [[ $extName == '.csfasta' ]];then
	if test -s $qualFileN; then
		qFileLineNo=$(wc -l $qualFileN | cut -f 1 -d ' ')
		echo $qFileLineNo
		if [ $iFileLineNo -eq $qFileLineNo ]; then
			firstReadN=$(head -n 1 $iFileN | cut -f 1 -d ' ' | cut -f 1 -d ',')
			firstReadQN=$(head -n 1 $qualFileN | cut -f 1 -d ' ' | cut -f 1 -d ',')
			if [ $firstReadN == $firstReadQN ]; then
				split -l $linesPerFile $qualFileN $oFileNPrefix'_' 
				splitFiles=$(ls | grep -P ^$oFileNPrefix"_"[a-z]+[a-z]$)
				for splitFile in $splitFiles;
				do mv ./$splitFile $splitFile$qualExtFileN; done
				exit
			else
				echo "$firstReadN != $firstReadQN"
				echo "The first read in the read file and the quality file are different."
			fi
		else 
			echo "$iFileN and $qualFileN have differnet line number. $iFileLineNo != $qFileLineNo"
			if [ $iFileLineNo -lt $qFileLineNo ]; then 
				checkExtraLines $qualFileN $linesPerRead 
			else
				checkExtraLines $qualFileN $linesPerRead  
				checkExtraLines $qualFileN $linesPerRead  
			fi
		fi
	else
		echo "No corresponding quality files are found"
	fi
fi