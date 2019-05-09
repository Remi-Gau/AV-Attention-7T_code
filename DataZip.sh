#!/bin/bash



data_folder=/media/remi/BackUp2/AV_Integration_7T_2/Subjects_Data/
dropbox_folder=/home/remi/Dropbox/PhD/Experiments/AV_Integration_7T/Subjects_Data/

clear

Subjects_List='02 03 04 06 07 08 09 11 12 13 14 15 16'
# Subjects_List='02'

for Subject in $Subjects_List;
do

	echo "\nMoving files files for subject $Subject \n"
	# cp $dropbox_folder/Subject_$Subject/Results/SVC/* $data_folder/Subject_$Subject/Results/SVC/
	# cp $dropbox_folder/Subject_$Subject/Results/Profiles/Surfaces/* $data_folder/Subject_$Subject/Results/Profiles/Surfaces/
	# cp $dropbox_folder/Subject_$Subject/ROI/* $data_folder/Subject_$Subject/ROI/

	# rm -r $data_folder/Subject_$Subject/Transfer/ROI


	echo "\nCompressing files for subject $Subject \n"

	# FileList="$(ls $data_folder/Subject_$Subject/BetaMapping/8Surf/*.vtk)"
	# echo $FileList
	# tar -zcvf $data_folder/Subject_$Subject/BetaMapping/8Surf/beta_vtk.tar.gz $FileList
	# rm $FileList

	# FileList="$(ls $data_folder/Subject_$Subject/BetaMapping/r4*.nii)"
	# for iFile in $FileList;
	# do
	# 	gzip -vf $iFile
	# done

	# tar -jcvf S3rBetaFiles.tar.bz2 $FileList
	# rm $FileList

	# FileList="$(ls S6rbeta*.nii)"
	# echo $FileList
	# tar -jcvf S6rBetaFiles.tar.bz2 $FileList
	# rm $FileList

	# mkdir /media/rxg243/BackUp2/AV_Integration_7T_2/Subjects_Data/Subject_$Subject/UpsampledBetas
	# mv /data/AV_Integration_2/Subjects_Data/Subject_$Subject/Transfer/S*rBetaFiles.tar.bz2 /media/rxg243/BackUp2/AV_Integration_7T_2/Subjects_Data/Subject_$Subject/UpsampledBetas

done
