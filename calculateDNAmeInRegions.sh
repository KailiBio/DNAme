!/bin/bash

# calculateDNAmeInRegions.sh

#########################
echo "\033[32m  @o@  @o@  @o@  @o@  @o@  @o@  @o@  @o@  @o@  @o@  @o@  @o@  \033[0m"
start=$(date +%Y-%m-%d\ %H:%M:%S)
echo "\033[32m  BEGIN@ "$start" \033[0m"
echo "\033[32m  enjoy~ \033[0m"
echo "\033[32m  @o@  @o@  @o@  @o@  @o@  @o@  @o@  @o@  @o@  @o@  @o@  @o@  \033[0m"
echo ""

#########################
help_info(){
	echo "usage:"
	echo "sh calculateDNAmeInRegions.sh <option>* [-x input_1] [-y input_2] [-r region] [-s strand] [-f fasta] [-c chrome_size] [-o output_path] [-p prefix]"
	echo ""
	echo "This file is for calculating DNA methylation level for given region (based on Bismark output)."
	echo "Supporting one file or 2 replicates, strand specific or without strand."
	echo ""
	echo "Arguments:"
	echo "-x Input DNAme file (bedMethyl file)."
	echo "-y Input file for replicate 2 (if with replications)."
	echo "-r The interested regions with uniq name on the 4th column."
	echo "   Should be strandard 6 column bed file."
	echo "-s Logical parameter (1 or 0). "
	echo "   1 means region file with strand Info, this code will considering strand while doing intersectBed."
	echo "-f The genome fasta file."
	echo "-c The chromosome size file."
	echo "-o The path for output file."
	echo "-p Prefix for intermediate files."
	echo ""
	echo "Any questions, please contact me."
	echo "		--Kaili Fan (fankaili.bio@gmail.com)"
	echo ""
}

if [ $# -lt 6 ];then
	help_info
	exit 1
fi

input_2="NA"

while getopts "x:y:r:s:f:c:o:p:" Arg
do
	case $Arg in
		x)	input_1=$OPTARG;;
		y)	input_2=$OPTARG
			echo "with replicates";;
		r)	region=$OPTARG;;
		s)	strand=$OPTARG;;
		f)	fasta=$OPTARG;;
		c)	chrome_size=$OPTARG;;
		o)	output_path=$OPTARG;;
		p)	prefix=$OPTARG;;
		?)	echo "Wrong parameter!!!"
			exit 1;;
	esac
done

#########################
cd ${output_path}

echo ${prefix}

# 1. count number of CpG sites
echo "counting CpG sites"
bedtools slop -i ${region} -g ${chrome_size} -b 1 > temp_${prefix}_region_add.bed ;
bedtools getfasta -fi ${fasta} -bed temp_${prefix}_region_add.bed > temp_${prefix}_fa.txt ;
cat temp_${prefix}_fa.txt | tr a-z A-Z | awk '{if(NR%2!=1){split($0,a,"CG");print length(a)-1}}' > ${prefix}_cg_count.txt ;
paste ${region} ${prefix}_cg_count.txt > ${prefix}_cg_count_merge.txt ;

if [ ${strand} == "1" ];then
	echo "strand-specific" ;

	if [ "$input_2" = "NA" ]; then

		# 2. remove low coverage CpG sites
		echo "remove CpG sites with less than 5 reads coverage" ;
		aveMethy_1=`awk '{FS=OFS="\t"}{if($10>5){sum+=$11}}END{print sum/(100*NR)}' ${input_1}` ;
		awk 'BEGIN{FS=OFS="\t"} {if($10>5){print $1,$2,$3,$4,$11/100,$6}}' ${input_1} > temp_${prefix}_raw_methy_1.txt ;

		# 3. calculate DNA methylation level of given region in each file
		echo "DNAme in files" ;
		intersectBed -a ${prefix}_cg_count_merge.txt -b temp_${prefix}_raw_methy_1.txt -s -wa -wb | \
		awk -v aveMethy="$aveMethy_1" 'BEGIN{FS=OFS="\t"}{sm[$4]+=$12;sn[$4]+=1}END{for(i in sm){if(sn[i]){print i,sn[i],sm[i]/sn[i]}else{print i,0,aveMethy}}}' \
		> ${prefix}_CG_methy_count_1.txt ;

		# 4. fix
		echo "fix..." ;
		awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=$3}else{if($7==0){print $0,"NA"}else if(($7!=0)&&(a[$4]==0)){print $0,-0.1}else{print $0,b[$4]}}}' \
		${prefix}_CG_methy_count_1.txt ${prefix}_cg_count_merge.txt > ${prefix}_cg_DNAme.txt ;

	else

		# 2. clean methyl data
		echo "remove CpG sites with less than 5 reads coverage" ;
		aveMethy_1=`awk '{FS=OFS="\t"}{if($10>5){sum+=$11}}END{print sum/(100*NR)}' ${input_1}` ;
		awk 'BEGIN{FS=OFS="\t"} {if($10>5){print $1,$2,$3,$4,$11/100,$6}}' ${input_1} > temp_${prefix}_raw_methy_1.txt ;
		aveMethy_2=`awk '{FS=OFS="\t"}{if($10>5){sum+=$11}}END{print sum/(100*NR)}' ${input_2}` ;
		awk 'BEGIN{FS=OFS="\t"} {if($10>5){print $1,$2,$3,$4,$11/100,$6}}' ${input_2} > temp_${prefix}_raw_methy_2.txt ;

		# 3. calculate DNA methylation level of given region in each file
		echo "DNAme in files" ;
		intersectBed -a ${prefix}_cg_count_merge.txt -b temp_${prefix}_raw_methy_1.txt -s -wa -wb | \
		awk -v aveMethy="$aveMethy_1" 'BEGIN{FS=OFS="\t"}{sm[$4]+=$12;sn[$4]+=1}END{for(i in sm){if(sn[i]){print i,sn[i],sm[i]/sn[i]}else{print i,0,aveMethy}}}' \
		> ${prefix}_CG_methy_count_1.txt ;
		intersectBed -a ${prefix}_cg_count_merge.txt -b temp_${prefix}_raw_methy_2.txt -s -wa -wb | \
		awk -v aveMethy="$aveMethy_2" 'BEGIN{FS=OFS="\t"}{sm[$4]+=$12;sn[$4]+=1}END{for(i in sm){if(sn[i]){print i,sn[i],sm[i]/sn[i]}else{print i,0,aveMethy}}}' \
		> ${prefix}_CG_methy_count_2.txt ;
		awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=$3}else{m=($3+b[$1])/2; print $1,a[$1],b[$1],$2,$3,m}}' ${prefix}_CG_methy_count_1.txt \
		${prefix}_CG_methy_count_2.txt > ${prefix}_CG_methy_count.txt ;

		# 4. fix
		echo "fix..." ;
		awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$7}else{if((a[$1]!=0)&&($2==0)&&($4==0)){print $1,-0.1}else if((a[$1]!=0)&&($2==0)&&($4!=0)){print $1,$5}else if((a[$1]!=0)&&($2!=0)&&($4==0)){print $1,$3}else{print $1,$6}}}' \
		${prefix}_cg_count_merge.txt ${prefix}_CG_methy_count.txt > ${prefix}_CG_methy_count_fix.txt ;
		awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{if($4 in a){print $0,a[$4]}else{if($5==0){print $0,"NA"}else{print $0,-0.1}}}}' \
		${prefix}_CG_methy_count_fix.txt ${prefix}_cg_count_merge.txt > ${prefix}_cg_DNAme.txt ;

	fi
else
	echo "without strand Info";

	if [ "$input_2" = "NA" ]; then

		# 2. remove low coverage CpG sites
		echo "remove CpG sites with less than 5 reads coverage" ;
		aveMethy_1=`awk '{FS=OFS="\t"}{if($10>5){sum+=$11}}END{print sum/(100*NR)}' ${input_1}` ;
		awk 'BEGIN{FS=OFS="\t"} {if($10>5){print $1,$2,$3,$4,$11/100,$6}}' ${input_1} > temp_${prefix}_raw_methy_1.txt ;

		# 3. calculate DNA methylation level of given region in each file
		echo "DNAme in files" ;
		intersectBed -a ${prefix}_cg_count_merge.txt -b temp_${prefix}_raw_methy_1.txt -wa -wb | \
		awk -v aveMethy="$aveMethy_1" 'BEGIN{FS=OFS="\t"}{sm[$4]+=$12;sn[$4]+=1}END{for(i in sm){if(sn[i]){print i,sn[i],sm[i]/sn[i]}else{print i,0,aveMethy}}}' \
		> ${prefix}_CG_methy_count_1.txt ;

		# 4. fix
		echo "fix..." ;
		awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=$3}else{if($7==0){print $0,"NA"}else if(($7!=0)&&(a[$4]==0)){print $0,-0.1}else{print $0,b[$4]}}}' \
		${prefix}_CG_methy_count_1.txt ${prefix}_cg_count_merge.txt > ${prefix}_cg_DNAme.txt ;

	else

		# 2. clean methyl data
		echo "remove CpG sites with less than 5 reads coverage" ;
		aveMethy_1=`awk '{FS=OFS="\t"}{if($10>5){sum+=$11}}END{print sum/(100*NR)}' ${input_1}` ;
		awk 'BEGIN{FS=OFS="\t"} {if($10>5){print $1,$2,$3,$4,$11/100,$6}}' ${input_1} > temp_${prefix}_raw_methy_1.txt ;
		aveMethy_2=`awk '{FS=OFS="\t"}{if($10>5){sum+=$11}}END{print sum/(100*NR)}' ${input_2}` ;
		awk 'BEGIN{FS=OFS="\t"} {if($10>5){print $1,$2,$3,$4,$11/100,$6}}' ${input_2} > temp_${prefix}_raw_methy_2.txt ;

		# 3. calculate DNA methylation level of given region in each file
		echo "DNAme in files" ;
		intersectBed -a ${prefix}_cg_count_merge.txt -b temp_${prefix}_raw_methy_1.txt -wa -wb | \
		awk -v aveMethy="$aveMethy_1" 'BEGIN{FS=OFS="\t"}{sm[$4]+=$12;sn[$4]+=1}END{for(i in sm){if(sn[i]){print i,sn[i],sm[i]/sn[i]}else{print i,0,aveMethy}}}' \
		> ${prefix}_CG_methy_count_1.txt ;
		intersectBed -a ${prefix}_cg_count_merge.txt -b temp_${prefix}_raw_methy_2.txt -wa -wb | \
		awk -v aveMethy="$aveMethy_2" 'BEGIN{FS=OFS="\t"}{sm[$4]+=$12;sn[$4]+=1}END{for(i in sm){if(sn[i]){print i,sn[i],sm[i]/sn[i]}else{print i,0,aveMethy}}}' \
		> ${prefix}_CG_methy_count_2.txt ;
		awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2;b[$1]=$3}else{m=($3+b[$1])/2; print $1,a[$1],b[$1],$2,$3,m}}' ${prefix}_CG_methy_count_1.txt \
		${prefix}_CG_methy_count_2.txt > ${prefix}_CG_methy_count.txt ;

		# 4. fix
		echo "fix..." ;
		awk '{FS=OFS="\t"}{if(NR==FNR){a[$4]=$7}else{if((a[$1]!=0)&&($2==0)&&($4==0)){print $1,-0.1}else if((a[$1]!=0)&&($2==0)&&($4!=0)){print $1,$5}else if((a[$1]!=0)&&($2!=0)&&($4==0)){print $1,$3}else{print $1,$6}}}' \
		${prefix}_cg_count_merge.txt ${prefix}_CG_methy_count.txt > ${prefix}_CG_methy_count_fix.txt ;
		awk '{FS=OFS="\t"}{if(NR==FNR){a[$1]=$2}else{if($4 in a){print $0,a[$4]}else{if($5==0){print $0,"NA"}else{print $0,-0.1}}}}' \
		${prefix}_CG_methy_count_fix.txt ${prefix}_cg_count_merge.txt > ${prefix}_cg_DNAme.txt ;

	fi
fi

rm temp_${prefix}_*



echo ""
echo ""
echo "\033[32m  Cheers!!!✌️ ✌️ ✌️ \033[0m"
echo ""
echo ""

#########################
echo "\033[32m ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^\033[0m"
end=$(date +%Y-%m-%d\ %H:%M:%S)
echo "\033[32m END@ "$end" \033[0m"
echo "\033[32m Time used: $((${SECONDS} / 3600))h $((${SECONDS} / 60))m \033[0m"
echo "\033[32m ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^ ^m^\033[0m"
