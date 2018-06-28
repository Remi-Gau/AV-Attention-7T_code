#!/bin/bash

# freesurfer-Linux-centos6_x86_64-stable-pub-v5.3.0

#export FREESURFER_HOME=/home/rxg243/Programs/freesurfer
#export SUBJECTS_DIR=/data/AV_Integration_7T_2/Subjects_Data
#source $FREESURFER_HOME/SetUpFreeSurfer.sh

clear

SUBJECTS_DIR=/data/AV_Integration_2/Subjects_Data/

Subjects_List='02'

echo "\n"
echo '########################################################################'
echo '#                                                                      #'
echo '#                            Creating ROIs                             #'
echo '#                                                                      #'
echo '########################################################################'


for SubjectInd in $Subjects_List
do

	echo "PROCESSING SUBJECT $SubjectInd\n\n"

	mkdir $SUBJECTS_DIR/Subject_$SubjectInd/Retinotopy/RegData 
	SourceImage=`ls $SUBJECTS_DIR/Subject_$SubjectInd/Retinotopy/MIPAV/copy_mean*.nii`
	
	RegFile=$SUBJECTS_DIR/Subject_$SubjectInd/Retinotopy/RegData/register.dat
	
	Overlay1=$SUBJECTS_DIR/Subject_$SubjectInd/Retinotopy/MIPAV/copy_Polar_MeanSineBetas.nii
	Overlay2=$SUBJECTS_DIR/Subject_$SubjectInd/Retinotopy/MIPAV/copy_Polar_MeanCosineBetas.nii

	echo "Source image: \n$SourceImage"
	
	echo "\n\n"
	echo '########################################################################'
	echo '#                            REGISTRATION                              #'
	echo '########################################################################'
	
	# Register source image to anatomical surface

	
	tkregister2 --mov $SourceImage --s Subject_$SubjectInd/Structural/FS --regheader --noedit --reg $RegFile
	
	#freeview -v $SUBJECTS_DIR/Subject_$SubjectInd/Structural/FS/mri/brainmask.mgz \
	#-v Subject_$SubjectInd/Structural/FS/mri/brainmask.mgz $SourceImage 
	#-reg $RegFile
	#-v Subject_$SubjectInd/Structural/FS/mri/brainmask.mgz $Overlay1 
	#-reg $RegFile

	tksurfer Subject_$SubjectInd/Structural/FS lh inflated -gray \
	-overlay $Overlay1 \
	-overlay-reg $RegFile \
	-overlay $Overlay2 \
	-overlay-reg $RegFile -tcl Freesurfer_Set_Thresh.tcl
	
	#tksurfer Subject_$SubjectInd/Structural/FS rh inflated -gray \
	#-overlay $Overlay1 \
	#-overlay-reg $RegFile \
	#-overlay $Overlay2 \
	#-overlay-reg $RegFile -tcl Freesurfer_Set_Thresh.tcl

	
	
	#mri_label2vol --label $SUBJECTS_DIR/Subject_$SubjectInd/Structural/FS_ROI/V1/Right/V1_rh_Polar.label \
	#--temp $SourceImage \
	#--reg $RegFile \
	#--fillthresh .5 \
	#--proj frac -1 1 .1 \
	#--subject $SubjectInd/Structural/FSL_LowDef \
	#--hemi rh \
	#--o $SUBJECTS_DIR/Subject_$SubjectInd/Structural/FS_ROI/V1/Right/V1/Right/V1_Ret_rh.nii
	
	#mri_label2vol --label $SUBJECTS_DIR/Subject_$SubjectInd/Structural/FS_ROI/V1/Left/V1_lh_Polar.label \
	#--temp $SourceImage \
	#--reg $RegFile \
	#--fillthresh .5 \
	#--proj frac -1 1 .1 \
	#--subject $SubjectInd/Structural/FSL_LowDef \
	#--hemi lh \
	#--o $SUBJECTS_DIR/Subject_$SubjectInd/Structural/FS_ROI/V1/Left/V1_Ret_lh.nii
	
done
