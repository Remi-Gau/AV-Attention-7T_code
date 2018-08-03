#!/bin/bash

pwd

clear

Subjects_List='02 03 04 06 07 08 09 11 12 13 14 15 16'

for Subject in $Subjects_List;
do


	echo "\nCompressing files for subject $Subject \n"
	
	cd /data/AV_Integration_2/Subjects_Data/Subject_$Subject/Transfer/



	FileList="$(ls S3rbeta*.nii)"

	echo $FileList

	tar -jcvf S3rBetaFiles.tar.bz2 $FileList

	rm $FileList



	FileList="$(ls S6rbeta*.nii)"

	echo $FileList

	tar -jcvf S6rBetaFiles.tar.bz2 $FileList

	rm $FileList



	mkdir /media/rxg243/BackUp2/AV_Integration_7T_2/Subjects_Data/Subject_$Subject/UpsampledBetas

	mv /data/AV_Integration_2/Subjects_Data/Subject_$Subject/Transfer/S*rBetaFiles.tar.bz2 /media/rxg243/BackUp2/AV_Integration_7T_2/Subjects_Data/Subject_$Subject/UpsampledBetas

	
done

