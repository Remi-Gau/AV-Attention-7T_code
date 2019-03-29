# AV-Attention-7T_code

Code to run the analysis of the [AV attention fMRI](https://github.com/Remi-Gau/AV_Attention-Presentation_code-fMRI) experiment (7T)

This still needs more commenting and documenting (I am still learning) but do reach out if you need any help with it.

A lot of the code is similar to that of the better documented [AVT experiment analysis](https://github.com/Remi-Gau/AVT_analysis).


## Dependencies for analysis:

You will need the following softwares to run part of the analysis.

| Softwares                                                           | Used version | Purpose                                                              |
|---------------------------------------------------------------------|--------------|----------------------------------------------------------------------|
| [FSL](https://fsl.fmrib.ox.ac.uk/fsl)                               | 5.0          | coregistration quality visualization                                 |
| [ANTs](http://stnava.github.io/ANTs/)                               | 2.1.0        | intersubject coregistration (MMSR)                                   |
| [JIST and the CBS tools](https://www.nitrc.org/projects/cbs-tools/) | 2 & 3.0.8    | segmentation, laminae definition, intersubject coregistration (MMSR) |
| [MIPAV](https://mipav.cit.nih.gov/)                                 | 7.0.1       | segmentation, laminae definition, intersubject coregistration (MMSR) |
| cosmetic                                                            |              | A1 ROI delineation                                                |


Many extra matlab functions are needed and are added to the path by the function code/subfun/Get_dependencies

| Matlab, toolbox and other dependencies                                                                                                            | Used version | Purpose                    |   |
|---------------------------------------------------------------------------------------------------------------------------------------------------|--------------|----------------------------|---|
| [Matlab](https://www.mathworks.com/products/matlab.html)                                                                                          | 2016a        |                            |   |
| SPM12                                                                                                                                             | v6685        | preprocessing, GLM, ...    |   |
| [SPM-RG](https://github.com/Remi-Gau/SPM-RG)                                                                                                      | N/A          | manual coregistration      |   |
| [nansuite](https://fr.mathworks.com/matlabcentral/fileexchange/6837-nan-suite)                                                                    | ???          |                            |   |
| [distributionPlot](https://fr.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m) | ???          | violin plots for matlab    |   |
| [plotSpread](https://fr.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot)                                          | ???          | plot datta spread          |   |
| [shadedErrorBar](https://fr.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar)                                            | ???          | shaded error bar           |   |
| [herrorbar](https://fr.mathworks.com/matlabcentral/fileexchange/3963-herrorbar)                                                                   | ???          | horizontal error bar       |   |
| [mtit](https://fr.mathworks.com/matlabcentral/fileexchange/3218-mtit-a-pedestrian-major-title-creator)                                            |              | main title for figures     |   |
| [matlab_for_CBS_tools](https://github.com/Remi-Gau/matlab_for_cbs_tools)                                                                          | NA           | import CBS-tools VTK files |   |
| [brain_colours](https://github.com/CPernet/brain_colours)                                                                                         | NA           | brain color maps           |   |
