#!/bin/sh

###DON'T FORGET: replace \r\n with \n###

module load r
module load plink/1.90 #have been having problems with just 'module load plink'

###The purpose of this script is to extract a list of individuals (from the database) from multiple different locations on the cluster
#It requires a list of IDs (usually generated from an extract that Ling Lin gives you, containing DBID [19052] and plate info [e.g. 55D02])
#-I started by making a list of all plate positions that I wanted to include, in the format that is in the plink .fam files (55_D02 format for some, 55D02 for others)
#-I called it IDlist.txt
#You must have the shell script in a directory above the directories that you will be searching through (could be the primary /srv/gsfs0/projects/mignot directory)

#grep the .fam (PLINK) file: (for example grep -f myplatepositions.txt plink.fam > samplestouse.txt)
grep -r --include=*.fam --file='IDlist.txt' ${PWD} > samplenames.txt

#can also grep with the date, may need it in order to pick the most updated file among duplicates
grep -r --include=*.fam --file='IDlist.txt' ${PWD} | awk -F: '{"stat -c %z "$1 | getline r; print r": "$0 }' > samplenamesWdate.txt

#awk file path and the (FID/IID), and use that as an extract list for the genotype file in plink: (for example by doing awk '{print $1}' samplenames.txt > extractlist)
awk '{print $1, $2}' samplenames.txt > extractlist.txt
awk '{print $4, $5 ";" $1, $2"\!"}' samplenamesWdate.txt > extractlistWdate.txt #added ; and ! for easier parsing, because I'm a sed noob
# sed -n 's/^.*:\(\)/\1/p' extractlist.txt gets what I want from the end of the string
# sed -n 's/\..*\(\)/\1/p' extractlist.txt gets what I want from the start of the string
sed -n 's/\(\/[^\.]*\)\..*:\([^:]*\)/\2 \1/p' extractlist.txt > PATHandID.txt
sed -n 's/\(\/[^\.]*\)\..*:\([^:]*\).*;\([^\;]*\)!/\2 \1 \3/p' extractlistWdate.txt > PATHandIDwDate.txt
	#What this sed is asking for
	#-n don't print every line automatically
	#s/.../ is the read process
	#\(\/[^\.]*\), everything between \( and \) will be remembered and returned for \1
	#	this is going to be the PLINK binary file names without .fam/.bed/.map
	#	essentially, it reads everything up to the first period in the file name (i.e. the start of the extension)
	#\..*: now read in the period and then read up to the colon (which is followed by the grepped sample ID)
	#\([^:]*\), everything between \( and \) will be remembered and returned for \2
	#	this is going to be the sample FID and IID pair
	#	remember everything that was just read before the colon based on .* reading
	#.*; now read up to the semicolon
	#\([^\;]*\), everything between \( and \) will be remembered and returned for \3
	#	this is going to be the time stamp for the file (date and time)
	#	remember everything that was just read before the semicolon based on .* reading
	#finish the line by read the !
	#return the desired data chunks sample/individual ID (\2) path to file (\1) date & time (\3)

cat <<EOF >removedups.R
samples <- read.table("PATHandIDwDate.txt", sep=" ")
#Make a POSIXct variable out of the date and time
samples[,"datetime"] <- as.POSIXct(paste(samples[,"V4"],samples[,"V5"]))
#select the most recent entry(ies) for a given sample/individual ID
mostrecent <- do.call(rbind, by(samples, samples[,"V1"], function(d) d[which.max(d[,"datetime"]),]))
#if duplicates entries still exist, eliminate them
if(anyDuplicated(mostrecent) != 0){
	mostrecent <- subset(mostrecent, !duplicated(mostrecent))
}
mrFAM <- as.character(mostrecent[,"V3"])
#ensure only extracting once from each unique PLINK binary set
if(anyDuplicated(mrFAM) != 0){
	mrFAM <- subset(mrFAM, !duplicated(mrFAM))
}
write.table(mrFAM, file="mostrecentFAMs.txt", quote=F, row.names=F, col.names=F)
write.table(mostrecent[,c("V1","V2")], file="SampleListFromPlates.txt", quote=F, row.names=F, col.names=F)
q()
EOF

R CMD BATCH removedups.R

#if repository directory doesn't exist, make it
mkdir -p PLINKextract

#Now, to go through all unique PLINK binaries and extract individuals from each location into PLINKextract

while read line
do
	#get PLINK binary name with path
	PLINKfiles=`echo $line`
	echo $PLINKfiles
	#get name of PLINK binaries for output
	source=`echo $line | awk -F "/" '{print $NF}'`
	echo $source
	
	plink --bfile $PLINKfiles --keep SampleListFromPlates.txt --allow-no-sex --make-bed --out PLINKextract/FROM_${source}
	
done<mostrecentFAMs.txt

