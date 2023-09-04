
# Relative Analysis of Fluorescence Localisation (RAFL) 

This command-line Python script executes image processing, segmentation, and data analysis to extract individual cell metrics for average fluorophore intensities across distinct cellular regions (nucleus, cytoplasm). 
The script supports image folders and enables cell gating based on nuclear stain intensity. Users can define their own gating thresholds or allow the script to calculate them automatically. 
The resulting single-cell metrics include relative fluorophore intensities per region and the percentage of nuclear fluorophore intensity, organized into Excel sheets by experimental conditions. A summary sheet consolidates single-cell metrics per condition

The script generates an image folder containing histograms depicting nuclear stain intensity distribution across all analyzed cells in the dataset, clearly marking the gating thresholds. Additionally, it creates overlays of nuclear and cytoplasmic masks on fluorophore images, alongside segmentation and gating process images for each field of view.

## USAGE

To run the python script, copy the script into your experiment folder. 
Images to be processed should be in a folder titled 'convertedtifs', which is within in the same folder as where this script is copied into. Else, change folder name to desired name on line 30

type python, then the name of this python script. 
Eg "python3 RAFL.py"

if using on a windows instead of macs, change the backslashes to slashes on line 30, 583, 583

## Parameters that can be changed: 
    line 15: cell growth size in pixels 
    line 16, 17: segmentation footprint and minimum distance for watershed segmentation using scipy package (see scipy documentation)
    line 24: gatebyfixed values : a tuple of bounds, or 'No' for automatic gating

if files are raw .nd2 files taken from the Nikon microscope, nd2totif.py can be run (command line script) to split each .nd2 into their 5 channels and save each channel as a .tif into a 'convertedtifs' folder. For this, the .nd2 images should be in a folder titled 'raw_images'

Script takes in 3 fluorophore images for each field of view : 
    w1 is the nuclear stain for segmentation
    w2 is fluorophore 1 (written as 'GFP')
    w3 is fluorophore 2 (written as 'RFP')

Name the files like so : (condition or well)-(site)_(wavelength)
eg: 1-1_w1 means condition1, site1, wavelegth1 (dapi)

if files are named a different way, please go to 
    line 373, 234 to change site, 
    line 374 to change well 
    line 231 to change nuclear stain wavelength
    line 488, 491, 503 to change fluorophore wavelengths


## CONTACT: 
    If you have any problems, questions, ideas or suggestions, please email cngsc@stanford.edu
