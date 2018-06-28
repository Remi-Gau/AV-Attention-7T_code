#!/bin/bash

# freesurfer-Linux-centos6_x86_64-stable-pub-v5.3.0

# export FREESURFER_HOME=/home/SHARED/Program/freesurfer
# export SUBJECTS_DIR=/home/SHARED/Experiment/AV_Integration_7T/Subjects_Data
#source $FREESURFER_HOME/SetUpFreeSurfer.sh

SUBJECTS_DIR=/home/SHARED/Experiment/AV_Integration_7T/Subjects_Data/

Subjects_List='01 02 03 06 08 09 10 11 12'

HS='lh rh'

#LabelsLists="V1 V2 G_temp_sup-G_T_transv G_temp_sup-Plan_tempo S_temporal_sup"
LabelsLists="V1 V2 G_temp_sup-G_T_transv"

echo "\n"
echo '########################################################################'
echo '#                                                                      #'
echo '#                            Creating ROIs                             #'
echo '#                                                                      #'
echo '########################################################################'


for SubjectInd in $Subjects_List
do

	echo "PROCESSING SUBJECT $SubjectInd\n\n"

	rm -r $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/ROI	      	      
	mkdir $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/ROI
	 
	
	UpSampImage=`ls $SUBJECTS_DIR$SubjectInd/Structural/400MNI/DownSamp_1mm_UNI_700MNI.nii`	
	
	# coregister 1 mm isometric image 
	tkregister2 --mov $UpSampImage --s $SubjectInd/Structural/FSL_LowDef --regheader --noedit --reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/ROI/register.dat

	for hs in $HS
	do	
		for label in $LabelsLists
		do
		
			# transform label into mask
			mri_label2vol --label $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/NewLabels/${hs}.${label}.label \
			--temp $UpSampImage \
			--reg $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/ROI/register.dat \
			--fillthresh .5 \
			--proj frac -1 1 .1 \
			--subject $SubjectInd/Structural/FSL_LowDef \
			--hemi ${hs} \
			--o $SUBJECTS_DIR$SubjectInd/Structural/FSL_LowDef/ROI/${label}_${hs}.nii

		done
	done	      
done
