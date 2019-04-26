#!/bin/bash



data_folder=/media/remi/BackUp2/AV_Integration_7T_2/Subjects_Data/

clear

Subjects_List='02 03 04 06 07 08 09 11 12 13 14 15 16'
# Subjects_List='02'

for Subject in $Subjects_List;
do

	echo "\nMoving files files for subject $Subject \n"
	mv $data_folder/Subject_$Subject/Transfer/Profiles/T1_06_Layers $data_folder/Subject_$Subject/Results/Profiles

	# rm -r $data_folder/Subject_$Subject/Transfer/Profiles

	mkdir $data_folder/Subject_$Subject/Results/SVM
	mv $data_folder/Subject_$Subject/Transfer/SVM/* $data_folder/Subject_$Subject/Results/SVM

	# rm -r $data_folder/Subject_$Subject/Transfer/SVM

	mkdir $data_folder/Subject_$Subject/ROI
	mkdir $data_folder/Subject_$Subject/ROI/MNI
	mkdir $data_folder/Subject_$Subject/ROI/MIPAV
	mv $data_folder/Subject_$Subject/Transfer/ROI/* $data_folder/Subject_$Subject/ROI
	mv $data_folder/Subject_$Subject/ROI_MIPAV/* $data_folder/Subject_$Subject/ROI/MIPAV
	mv $data_folder/Subject_$Subject/ROI_MNI/* $data_folder/Subject_$Subject/ROI/MNI

	# rm -r $data_folder/Subject_$Subject/Transfer/ROI

	mv $data_folder/Subject_$Subject/Results/Results/Profiles/Surfaces/* $data_folder/Subject_$Subject/Results/Profiles/Surfaces/

	mv $data_folder/Subject_$Subject/Results/Results/Profiles/Surfaces/Cdtions* $data_folder/Subject_$Subject/Results/Profiles/Surfaces/Cdtions

	mv $data_folder/Subject_$Subject/Transfer/* $data_folder/Subject_$Subject/BetaMapping

	echo "\nCompressing files for subject $Subject \n"

	# FileList="$(ls $data_folder/Subject_$Subject/BetaMapping/8Surf/*.vtk)"
	# echo $FileList
	# tar -zcvf $data_folder/Subject_$Subject/BetaMapping/8Surf/beta_vtk.tar.gz $FileList
	# rm $FileList

	FileList="$(ls $data_folder/Subject_$Subject/BetaMapping/r4*.nii)"
	for iFile in $FileList;
	do
		gzip -vf $iFile
	done

	# tar -jcvf S3rBetaFiles.tar.bz2 $FileList
	# rm $FileList
	#
	# FileList="$(ls S6rbeta*.nii)"
	# echo $FileList
	# tar -jcvf S6rBetaFiles.tar.bz2 $FileList
	# rm $FileList
	#
	# mkdir /media/rxg243/BackUp2/AV_Integration_7T_2/Subjects_Data/Subject_$Subject/UpsampledBetas
	# mv /data/AV_Integration_2/Subjects_Data/Subject_$Subject/Transfer/S*rBetaFiles.tar.bz2 /media/rxg243/BackUp2/AV_Integration_7T_2/Subjects_Data/Subject_$Subject/UpsampledBetas
	#
	#
done
