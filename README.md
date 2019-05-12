# AV-Attention-7T_code

Code to run the analysis of the [AV attention fMRI](https://github.com/Remi-Gau/AV_Attention-Presentation_code-fMRI) experiment (7T)

This still needs more commenting and documenting (I am still learning) but do reach out if you need any help with it.

A lot of the code is similar to that of the better documented [AVT experiment analysis](https://github.com/Remi-Gau/AVT_analysis).

## Data
Some of the data is available on the [open science framework project associated to this repo](https://osf.io/7ka5j/) as `.mat` or `.csv` files. The raw data (in a BIDS compatible format) of this project are available upon request.

## Dependencies

You will need the following softwares to run part of the analysis.

| Softwares                                                           | Used version | Purpose                                                              |
|---------------------------------------------------------------------|--------------|----------------------------------------------------------------------|
| [FSL](https://fsl.fmrib.ox.ac.uk/fsl)                               | 5.0          | coregistration quality visualization                                 |
| [ANTs](http://stnava.github.io/ANTs/)                               | 2.1.0        | intersubject coregistration (MMSR)                                   |
| [JIST and the CBS tools](https://www.nitrc.org/projects/cbs-tools/) | 2 & 3.0.8    | segmentation, laminae definition, intersubject coregistration (MMSR) |
| [MIPAV](https://mipav.cit.nih.gov/)                                 | 7.0.1        | segmentation, laminae definition, intersubject coregistration (MMSR) |
| cosmetic                                                            |              | A1 ROI delineation                                                   |


Many extra matlab functions from github and the mathwork file exchange are needed and are added to the path by the function `code/subfun/Get_dependencies`. Yeah this weird, tiring and cumbersome but that's matlab weirdness for you (“_And this why we can’t have nice things. Have you heard of [python](http://python.org)?_”)

| Matlab, toolbox and other dependencies                                                                                                            | Used version | Purpose                    |
|---------------------------------------------------------------------------------------------------------------------------------------------------|--------------|----------------------------|
| [Matlab](https://www.mathworks.com/products/matlab.html)                                                                                          | 2016a        |                            |
| SPM12                                                                                                                                             | v6685        | preprocessing, GLM, ...    |
| [SPM-RG](https://github.com/Remi-Gau/SPM-RG)                                                                                                      | NA          | manual coregistration      |
| [nansuite](https://fr.mathworks.com/matlabcentral/fileexchange/6837-nan-suite)                                                                    | V1.0.0       |                            |
| [distributionPlot](https://fr.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m) | v1.15.0      | violin plots for matlab    |
| [plotSpread](https://fr.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot)                                          | v1.2.0       | plot datta spread          |
| [shadedErrorBar](https://fr.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar)                                            | v1.65.0      | shaded error bar           |
| [herrorbar](https://fr.mathworks.com/matlabcentral/fileexchange/3963-herrorbar)                                                                   | V1.0.0       | horizontal error bar       |
| [mtit](https://fr.mathworks.com/matlabcentral/fileexchange/3218-mtit-a-pedestrian-major-title-creator)                                            | v1.1.0       | main title for figures     |
| [matlab_for_CBS_tools](https://github.com/Remi-Gau/matlab_for_cbs_tools)                                                                          | NA           | import CBS-tools VTK files |
| [brain_colours](https://github.com/CPernet/brain_colours)                                                                                         | NA           | brain color maps           |

## Reproduce the figures from the paper
You should be able to reproduce the laminar profile figures from the paper by using the following scripts on the CSV files available on [OSF](https://osf.io/7ka5j/).
- `BOLDProfiles/make_figures_BOLD.m`
- `MVPA/make_figures_MVPA.m`

## Data analysis workflow

I indicate here the different folders where the code is kept. I try to indicate and in which order the scripts (or other manual interventions) have to be  run.

**Preprocessing of EPIs: code/preprocess/**
1. `Preprocess_01_CreateVDM.m` : creates the voxel displacement map using the fieldmap
2. `Preprocess_02_RealignAndUnwarp.m` : realign and unwarp the EPIs

**Running subject level GLM: code/ffx/**
1. `Analysis_FFX_Block.m` : runs the subject level GLM. It must first be run a first time on smoothed images to get an inclusive mask (GLM-mask) that will be used for a second pass.

**Preprocessing anatomical: code/cbs/ or sub-xx/code/cbs/**
`segment-layer.LayoutXML` : high-res segmention and layering using the CBS tools
