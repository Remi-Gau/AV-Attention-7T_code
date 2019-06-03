!/bin/bash

# small script to compute the smoothness of the mean EPI after realignement/unwarp
# and after upsmabling at 0.4 mm  to compare it with that of the initial time series

#docker run -it --rm \
#-v '/home/remi/github/AV-Attention-7T_code':/code \
#-v /media/remi/BackUp2/AV_Integration_7T_2/derivatives/spm12/:/data:ro \
#-v /media/remi/BackUp2/AV_Integration_7T_2/derivatives/:/out \
#afni/afni

# docker run -it --rm \
# -v '/home/remi/github/AV-Attention-7T_code':/code \
# -v /media/remi/BIDS/AV_Att/rawdata:/data:ro \
# -v ~/:/out \
# afni/afni


Subjects_List='02 03 04 06 07 08 09 11 12 13 14 15 16'

for Subject in $Subjects_List;
do
  file_ls=`ls /data/sub-$Subject/*/func/*attention*bold.nii.gz`

  # echo $file_ls

  # mkdir /out/afni/
  # mkdir /out/afni/sub-$Subject

  for File in $file_ls;
  do

    FileName="$(echo $File | cut -d/ -f 6 | cut -d. -f 1)"

    # echo $FileName

    3dFWHMx -automask \
    -out /out/afni/sub-$Subject/FWHM_${FileName}.txt \
    -input $File
  done

  #3dFWHMx -automask \
  #-out /out/afni/sub-$Subject/ses-1/func/FWHM_meanUR.txt \
  #-input /data/sub-$Subject/ses-1/func/meanURsub-${Subject}_ses-1_task-audiovisualattention_run-01_bold.nii

  #3dFWHMx -automask \
  #-out /out/afni/sub-$Subject/ses-1/func/FWHM_r4meanUR.txt \
  #-input /data/sub-$Subject/ses-1/func/r4meanURsub-${Subject}_ses-1_task-audiovisualattention_run-01_bold.nii.gz

done
