#!/bin/bash

cd ..

clear

Subjects_List='02 03 04 07 08 09 11 12 13 15 16'

for Subject in $Subjects_List;
do


	echo "\nCompressing files for subject $Subject \n"

	mkdir /media/rxg243/BackUp2/AV_Integration_7T_2/Subjects_Data/Subject_$Subject/UpsampledBetas
	
	cd /data/AV_Integration_2/Subjects_Data/Subject_$Subject/Transfer/

	mv -v /data/AV_Integration_2/Subjects_Data/Subject_$Subject/Transfer/rbeta*.nii /media/rxg243/BackUp2/AV_Integration_7T_2/Subjects_Data/Subject_$Subject/UpsampledBetas

	mv -v /data/AV_Integration_2/Subjects_Data/Subject_$Subject/Transfer/rcon*.nii /media/rxg243/BackUp2/AV_Integration_7T_2/Subjects_Data/Subject_$Subject/UpsampledBetas

	
done

