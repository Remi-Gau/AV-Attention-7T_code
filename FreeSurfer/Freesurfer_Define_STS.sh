#!/bin/bash

# freesurfer-Linux-centos6_x86_64-stable-pub-v5.3.0

#export FREESURFER_HOME=/home/SHARED/Program/freesurfer
#export SUBJECTS_DIR=/home/SHARED/Experiment/AV_Integration_7T/Subjects_Data
#source $FREESURFER_HOME/SetUpFreeSurfer.sh

clear

SUBJECTS_DIR=/home/SHARED/Experiment/AV_Integration_7T/Subjects_Data/

Subjects_List='01 02 03 06 08 09 10 11 12'
Subjects_List='01 11'

echo "\n"
echo '########################################################################'
echo '#                                                                      #'
echo '#                            Creating ROIs                             #'
echo '#                                                                      #'
echo '########################################################################'


for SubjectInd in $Subjects_List
do

  	# ?h.G_temp_sup-G_T_transv.label : HG
  	
  	# ?h.S_temporal_sup.label: STS
  	
  	# ?h.G_temp_sup-Plan_tempo.label : PT

	echo "PROCESSING SUBJECT $SubjectInd\n\n"
  	
  	# Visualize to draw ROI
	tksurfer $SubjectInd/Structural/FSL_LowDef lh pial	
  	
	tksurfer $SubjectInd/Structural/FSL_LowDef rh pial	
	
	
	# transform label into mask
	#UpSampImage=`ls $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/DownSamp_1mm_UNI_700MNI.nii`
	
	#mkdir $SUBJECTS_DIR$SubjectInd/Structural/400MNI/RegData
	
	#tkregister2 --mov $UpSampImage --s $SubjectInd/Structural/FSL_LowDef --regheader --noedit --reg $SUBJECTS_DIR$SubjectInd/Structural/400MNI/RegData/register.dat
	
	
	#mri_label2vol --label $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/NewLabels/STS_Post_lh.label \
	#--temp $UpSampImage \
	#--reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register_2.dat \
	#--proj frac -1 1 .1 \
	#--subject $SubjectInd/Structural/FSL_LowDef \
	#--hemi lh \
	#--o $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/ROI/STS_Post_lh.nii

	#mri_label2vol --label $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/NewLabels/STS_Post_Dorsal_lh.label \
	#--temp $UpSampImage \
	#--reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register_2.dat \
	#--proj frac -1 1 .1 \
	#--subject $SubjectInd/Structural/FSL_LowDef \
	#--hemi lh \
	#--o $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/ROI/STS_Post_Dorsal_lh.nii
	
	
	#mri_label2vol --label $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/NewLabels/STS_Post_rh.label \
	#--temp $UpSampImage \
	#--reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register_2.dat \
	#--proj frac -1 1 .1 \
	#--subject $SubjectInd/Structural/FSL_LowDef \
	#--hemi rh \
	#--o $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/ROI/STS_Post_rh.nii
	
	#mri_label2vol --label $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/NewLabels/STS_Post_Dorsal_rh.label \
	#--temp $UpSampImage \
	#--reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/qT1/RegData/register_2.dat \
	#--proj frac -1 1 .1 \
	#--subject $SubjectInd/Structural/FSL_LowDef \
	#--hemi rh \
	#--o $SUBJECTS_DIR/$SubjectInd/Structural/FSL_LowDef/ROI/STS_Post_Dorsal_rh.nii
	

done
