
# RAFL: Relative Analysis of Fluorescence Localisation

`RAFL.py` is a command-line Python script that executes image processing, segmentation, and data analysis to extract individual cell metrics for average and relative fluorophore intensities across distinct cellular regions (nucleus, cytoplasm). 

## Table of contents 
[Overview](https://github.com/cngsc/RAFL#overview)
[System Requirements](https://github.com/cngsc/RAFL#system-requirements)
[Installation Guide](https://github.com/cngsc/RAFL#installation-guide)
[Demo](https://github.com/cngsc/RAFL#demo)
[Parameters that can be changed](https://github.com/cngsc/RAFL#parameters-that-can-be-changed)
[IMPORTANT: Naming your input](https://github.com/cngsc/RAFL#important-naming-your-input)
[Contact](https://github.com/cngsc/RAFL#contact)

## Overview
### RAFL.py
`RAFL.py`  supports image folders and enables cell gating based on nuclear stain intensity. Users can define their own gating thresholds or allow the script to calculate them automatically (see usage). 

The resulting single-cell metrics include relative fluorophore intensities per region and the percentage of nuclear fluorophore intensity, organized into Excel sheets by experimental conditions. A summary sheet consolidates single-cell metrics per condition. 

The script generates an image folder containing histograms depicting nuclear stain intensity distribution across all analyzed cells in the dataset, clearly marking the gating thresholds. Additionally, it creates overlays of nuclear and cytoplasmic masks on fluorophore images, alongside segmentation and gating process images for each field of view. 

### nd2totif.py
To aid analysis of raw image data obtained using a Nikon microscope with file extension .nd2, `nd2totif.py` code can be run from the command line to split each .nd2 file into the 5 flurophore channels and save each channel for each image as a .tif into a 'convertedtifs' folder. For this, the .nd2 images should be in a folder titled 'raw_images'. See demo for more details. Depending on the microscope used, the number of channels will differ. 

For the demo set, images were taken with the Nikon A1R cnfocal microscope with laser settings:  
405-nm violet laser, 488-nm blue laser, 561-nm green laser and 639-nm red laser  
w1: DAPI w2: GFP, w3: RFP, w4: CY5, w5: TD 

## System Requirements
### OS Requirements
These command line scripts require only a standard computer with enough RAM to support the size of the dataset to be analysed by the user. 
The command line scripts have been tested on the following systems : 

MacOSX: Ventura 13.4.1  
Windows: Windows Ver 10.0.22621.1848  

### Python Dependencies
Before utilising, users should have Python version 3.0 or higher, see https://packaging.python.org/en/latest/tutorials/installing-packages/

To setup the necessary dependencies, run `setup.py` from the command line. This will install any packages necessary, if not already present. 

`RAFL.py` depends on the following Python libraries:
```
numpy
scipy
skimage
imageio
pandas
csv
matplotlib
```
`nd2totif.py` depends on the following Python libraries:
```
numpy
skimage
matplotlib
nd2
```
## Installation Guide 
RAFL is implemented in Python ver 3.0 and can be installed by downloading the folder, and copying the .py files into your project folder file (see Usage)

## Demo 
### Instructions
Within the demo folder, there should be a dataset of .nd2 images in a folder 'raw_images', and 3 .py files.Make sure that your current working directory is the DEMO folder. Code will be run here and output folders also generated here. 

_if using windows, replace **python3** with **python** in the following code snippets_

First, check and install the required dependencies if not already present  
```python3 setup.py```  

Next, unpack the .nd2 images in the raw_images folder, and split them to their individual channels as .tif  
```python3 nd2totif.py```  
Now you should have a new folder generated "convertedtifs" that have each .nd2 image now split into their 5 individual channel .tifs

Next, run the segmentation and analysis code  
```python3 RAFL.py```
### Expected outputs
After the code has been run, the outputs are: 
1) A new folder labelled "Overlay_and_segmentation_images" will be generated, containing the gating histogram, nuclear and cytoplasmic masks overlaid onto the respective fluorophore images, and segmentation process images from the nuclear stain fluorophore images. 
2) An excel file with the label currentfoldername_GFP-RFP_percellcombinations.xlsx.
This splits each condition into a seperate excel sheet, and tabulates the respective metrics for each individual cell analysed for that condition (index is well-site-nuclearlabelnumber)
3) An excel file with the label currentfoldername_Results.xlsx
This excel file has summaries of the respective metrics for each condition (index = condition) for each fluoropohore, and the following excel sheets each individual metric for individual cells (tabulated in rows) grouped into conditions (columns). 
### Expected run time for demo set:
on MacOS Ventura 13.4.1: 82 seconds  
on windows x64 ver 10.0.22621.1848 : 166 seconds

## Usage

### Command line usage 
**Quick guide**
1) Download the folder from https://github.com/cngsc/RAFL
2) Copy the three .py scripts into your experiment folder. 
3) Run setup.py from command line to install dependencies (make sure experiment folder is your working directory)  
```python3 setup.py```
4) Ensure that images to be processed are in a folder "convertedtifs" in the current experiment folder where code is present. Else, change folder name to desired name on line 30 of RAFL.py
5) Ensure that your tif image naming system follows the required naming system (see Naming your input section)
6) Run RAFL.Py from the command line  
```"python3 RAFL.py"```

## Parameters that can be changed:
**Cell growth size in pixels [line 15]**  
This determines how much distance in pixels the nuclear mask is circumferentially grown from the nucleus to define "cell area". Users can visualise 

**Segmentation footprint and minimum distance [line 16, 17]**  
These affect the segmentation of nuclei using using the skimage package. (see peak_local_max under skimage documentation https://scikit-image.org/docs/stable/api/skimage.feature.html#skimage.feature.peak_local_max)
  
**Defining gates [line 24]**  
Here, we set the bounds to gate out dead, dying cells or cells by gating out nuclei that have intense nuclear stain intensities over a small area. Gating histogram is plotted in the ouput of `RAFL.py` and can be visualised in the generated 'Overlay_and_segmentation_images' folder. 

Users can choose to either use:
1. Automatic gating  
```gatebyfixedvalues = 'No' ```
2. Define their own gates as a tuple  
```gatebyfixedvalues = (leftboundvalue, rightboundvalue) ```

## IMPORTANT: Naming your input
`RAFL.py` takes in 3 fluorophore images (wavelengths) for each field of view :

w1 is the nuclear stain for segmentation  
w2 is fluorophore 1 (written as 'GFP')  
w3 is fluorophore 2 (written as 'RFP')  

and also takes into account multiple sites (field of views) for a single condition (well)

Name the files like so:  
**(condition or well)-(site)_(wavelength)**  
eg: 1-1_w1 means condition1, site1, wavelength1 (dapi)

if files are named a different way, please go to:  
    line 373, 234 to change site  
    line 374 to change well  
    line 231 to change nuclear stain wavelength  
    line 488, 491, 503 to change fluorophore wavelengths  


## Contact

If you have any problems, questions, ideas or suggestions, please email christinengsc97@gmail.com
