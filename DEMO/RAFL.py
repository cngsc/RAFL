import csv
import time
import pandas as pd
import scipy
import skimage
from skimage import io
import imageio as imageio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import sys,os
np.set_printoptions(threshold=np.inf)

# here we define main conditions that can be written to an excel sheet later (reports the analysis)
cell_growth_size_pixels = 6
segmentation_footprint = 8
segmentation_minimum_distance = 18
startTime = time.time()

# set your average nuclear dapi intensity fixed limits here.
# this is the average dapi intensity per pixel within the nuclues. 
# if you have set limits that you know you want to use, gatebyfixed values = (leftbound value, rightbound value)
# else, set it to gatebyfixedvalues = 'No', and the automated gating will be run - see README for details
gatebyfixedvalues = (0, 3000)


#Create an output folder for which images and data will be written to
path1 = os.getcwd()
path = os.path.join(path1, 'convertedtifs')

# make a directory to house upcoming generated overlay and segmentation images, overwrite if it already exists
outputimages_dir = os.path.join(path1, 'Overlay_and_segmentation_images') 
os.makedirs(outputimages_dir, exist_ok=True) 


def segmentation(img):
    # Takes in an image matrix, returns binary image, distance transform, segmented nuclei labels and mask
    bImg = blur(img)
    T = skimage.filters.threshold_otsu(bImg)     
    fImg = binary(bImg, T)
    fImg = fill_hole(fImg)
    dist_visual, labels, mask = watershed(fImg, segmentation_footprint,segmentation_minimum_distance)
    return fImg, dist_visual, labels, mask

def blur(img):
    # takes in an img, returns a gaussian blurred image
    return skimage.filters.gaussian(img, sigma = 1.0)
    
def fill_hole(fImg):
    # fills any hole within an image
    return scipy.ndimage.binary_fill_holes(fImg)

def binary(bImg, T):
    # Given an input of an image and threshold T, 
    # return a binary mask where anypixel above T is 1 and anything below is 0(background)
    return (bImg > T).astype(int)

def watershed(fImg, segmentation_footprint,segmentation_minimum_distance):
    # takes in an image, returns dist_visual, labels arrays. 
    # identifies where nuclei peaks are, then utilise watershed to determine nuclei edges. 
    # returns a distance transform array(dist_visual) and an an array where all pixels belonging to the same nuclei are labelled with the same label
    # pixels of different nuclei will have different labels 
    dist = scipy.ndimage.distance_transform_edt(fImg)
    dist_visual = dist.copy()
    coords = skimage.feature.peak_local_max(dist,footprint=np.ones((segmentation_footprint, segmentation_footprint)), min_distance=segmentation_minimum_distance, labels=fImg) #returns an array of coord tuples (y, x) of the maxima of dist 
    coordsT = coords.T #returns 2d array: first row is Y coords of dist max, second row is X coords of dist max 
    # initiate mask of TF in same array shape of pixels. 
    mask = np.zeros(dist.shape, dtype=bool)
    mask[tuple(coordsT)] = True # basically now you rewrite array of mask: for every tuple(y,x), mask[y][x] = true
    markers, nlabels = scipy.ndimage.label(mask, structure=np.ones((3, 3)))

    # watershed assumes markers represent local minima 
    # need to invert distance transform image (-dist) so light pixels = high elevation, dark pixel = low elevation
    # basically each dapi image, labels start from 0 (background), and go up to the number that is the total number of nuclei in that image 
    # all pixels belonging to one nuclei will have the same label number, eg each pixel belonging to first nucleus is labelled '1' in the labels array. 
    labels = skimage.segmentation.watershed(-dist, markers, mask=fImg, watershed_line=True)  
    nuclei_count=(nlabels)

    return dist_visual, labels, mask

def calculate_areas(labels):
    #takes in labels array
    #returns a dictionary with label number: area in micrometers
    dareas = {} # initiate dictionary to store label : area count values 
    nuc, counts = np.unique(labels, return_counts=True)
    realcounts= 0.4*counts # converting pixel numbers for each nucleus to micrometeres  (1 pixel = 0.4 micrometers)
    dareas = {nuc[i]: realcounts[i] for i in range(1, len(nuc))} # skip 0 because that is background pixel
    return dareas

def calculate_intensities(img, labels):
    # returns a dictionary of nuc label: intensity sums 
    # img contains pixel intensities at each x, y location in that fluorophore image (DAPI nuclear stain)
    # labels contains labels of each nucleus 
    rows = np.shape(img)[0]
    cols = np.shape(img)[1]

    di  = {} # initiate a dictionary to hold intensity sums per labeled nucleus 

    # loop through rows, columns of labels
    for y in range(rows):
        for x in range(cols):
            if labels[y][x] != 0: #avoid label0 because that is background pixel 
                nuc = labels[y][x]
                if nuc not in di.keys():
                    di[nuc] = 0           # initiate count for intensities sum for label if it doesnt already exist
                di[nuc] += img[y][x]   #now that you have initiated count, add the pixel value from img at that x y position for that label of nuc
    return di

def average_intensities(dareas, di):
    ai = {k: di[k] / float(dareas[k]) for k in dareas if k in di}
    return ai

def get_kde_norm_bounds(l_adi):
    # given a list of average intensity values(l_adi), returns a KDE (x, y)
    # plots a normal distribution based on the density histogram
    # and bounds calulated off curve 
    kde = scipy.stats.gaussian_kde(l_adi, bw_method="scott")

    x = np.linspace(np.min(l_adi), np.max(l_adi))
    y = kde.pdf(x)

    # calculate the mean and variance of the peak in the KDE plot 
    # obtain the x,y values of the max peak using index at max y value
    peak_mean = x[np.argmax(y)]  
    peak_y = y[np.argmax(y)]

    # get the half full peak x values 
    half_peak_y = peak_y/2 
    halffullheight_indices = np.where(y >= half_peak_y) # these are the indices in the y array where y is above half_peak_y
    halffullheight_xvals = x[halffullheight_indices]
    x_diff = max(halffullheight_xvals) - min(halffullheight_xvals) # this is FHWM

    # now convert this x_diff(FHWM) to standard dev
    standard_dev = x_diff/2.355

    # generate the values for the normal distribution
    x_norm = np.linspace(peak_mean - 3*standard_dev, peak_mean + 3*standard_dev, 100)
    y_norm = scipy.stats.norm.pdf(x_norm, peak_mean, standard_dev)
    leftboundx = x_norm[0]
    rightboundx = x_norm[-1]
    return x, y, x_norm, y_norm, leftboundx, rightboundx

def plot_histograms(l_adi, x, y, x_norm, y_norm, leftboundx, rightboundx):
    # takes in dictionary of values, x, y from kde function, x_norm,y_norm values for the normal distribution,
    # and the leftbounds rightbounds obtained from the normal distribution function limits 
    
    #returns a plot of histogram with kde, norm distribution plotted, and bounds labelled
    binlist = np.linspace(0, 10000, 200)
    fig = plt.figure()
    plt.hist(l_adi, bins=binlist, label='Histogram', density=True, zorder=0)
    #plot the KDE plot, local minima
    plt.plot(x, y, '--', c='r',label ='KDE plot', zorder=2)
    # plot the norm plot 
    plt.plot(x_norm, y_norm, '-.', c='k', label='Normal Distribution', zorder=2)
    nnuc = len(l_adi)
    ax_labs = ['Cumulative avg int vs area across images (%s nuclei) ' %nnuc, 'Average intensity per area', 'Density']
    plt.annotate("Leftbound x: \n %s" %leftboundx, xy=(leftboundx, 0), xytext =(0.60, 0.5), textcoords = "figure fraction", arrowprops=dict(facecolor='black', width=1, headlength=5),horizontalalignment='center', verticalalignment='top')
    plt.annotate("Rightbound x: \n %s" %rightboundx,  xy=(rightboundx, 0), xytext =(0.60, 0.4), textcoords = "figure fraction", arrowprops=dict(facecolor='black', width=1, headlength=5), horizontalalignment='center', verticalalignment='top')
    plt.title(ax_labs[0])
    plt.xlabel(ax_labs[1])
    plt.ylabel(ax_labs[2])
    plt.legend()
    
    plt.savefig(os.path.join(outputimages_dir, 'gating histogram'), bbox_inches="tight", dpi=720)
    plt.close(fig)

def plot_overlay(background, masks, titles, filename, channel):
    fig = plt.figure()
    fig.set_size_inches(16, 12)
    nimages = len(masks)

    for i in range(len(masks)):
        masked = np.ma.masked_where(masks[i] == 0, masks[i])
        plt.subplot(1,2,i+1) #row, column, position
        plt.imshow(background, cmap="gray", interpolation='none') # first plot the image
        plt.imshow(masked, cmap="nipy_spectral", interpolation='none', alpha=0.18)
        plt.title(titles[i], fontsize=30)
        plt.xticks([]), plt.yticks([])
        prefix_fn = filename + '_' + channel
        plt.savefig(os.path.join(outputimages_dir, '%s_overlay.png' %prefix_fn))
    plt.close(fig)

def plotimages(images, titles, filename):
# 
    fig = plt.figure()
    fig.set_size_inches(16, 12)
    nimages = len(images)
    if (nimages % 2) != 0: # if it is an odd number 
        col_plot = int((nimages +1)/2)
    else:
        col_plot = int(nimages/2)
    for i in range(nimages):
        plt.subplot(2,col_plot,i+1)
        if "Original" in titles[i]:
            cmap="gray"
        else:
            cmap="turbo" 
        plt.imshow(images[i], cmap)
        plt.title(titles[i], fontsize=20)
        plt.xticks([]), plt.yticks([])
        plt.savefig(os.path.join(outputimages_dir,'%s.png' %filename ))
    plt.close(fig)

def shrink_nucleus(nuc):
    # this shrinks the prelabelled nuclei by one pixel circumference everytime it is run
    rows = np.shape(nuc)[0]
    cols = np.shape(nuc)[1]
    newnuc = nuc.copy()
    for y in range(1, rows-1): # avoid boundary conditions
        for x in range(1, cols-1):
            currlabel = nuc[y][x]
            if nuc[y][x+1] != currlabel or nuc[y][x-1] != currlabel or nuc[y+1][x] !=currlabel or nuc[y-1][x] != currlabel:
                newnuc[y][x] = 0
    return newnuc

def initialize_all_labels(path, cell_growth_size_pixels):
    # takes in path of folder of images 
    # returns a list of average intensities for all nuclei across images (l_adi) for later calculations of bounds 
    # returns adi_dict: a dictionary of adi dictionaries(average intensity per label) where key = site, value = dictionary of average dapi intesnities for each label within that site
    # returns nlabels_dict : a dictionary of nuclei labels (matrix) 
    # returns clabels_dict : a dictionary of cell labels (matrix)

    l_adi = []   #iniate a list that will store the avg int per area of each nuclei as you loop through images
    dapi_dict = {} 
    nlabels_dict = {} # dictionary that will hold the labels of nuclei for each image
    clabels_dict = {} # dictionary that will hold the labels of cells (grown from nuclei via expandlabel) for each image. 
    adi_dict = {} # dictionary that holds adi dicts for each label
    da_dict = {}
    for file in os.listdir(path):
        if not file.endswith('w1.TIF'): #only if it is a dapi file
            continue

        site = str(file[-13:-7]) # obtain site title for each well image eg "1-1"
        filepath = os.path.join(path, file)
        img = io.imread(filepath, as_gray=True)
        print(file)

        # do the segmentation, build the dictionaries of dareas and di and average intensities dictionary (adi) 
        fImg, dist_visual, labels, mask = segmentation(img)
        da = calculate_areas(labels)
        di = calculate_intensities(img, labels)
        adi = average_intensities(da, di)

        # start building a list that combines all the avgintperarea of all nuclei across all images
        l_adi += list(adi.values()) 

        # store these in their respective large overarching libraries where key = site
        # dictionary nlabels holds the labels of each nuclei from nuclei segmentation (used to feed into expandlabels later)
        nlabels_dict[site] = labels
        adi_dict[site] = adi
        da_dict[site] = da
        
        # now call expand labels to draw a cell area around the nucleus
        cell = skimage.segmentation.expand_labels(labels, distance = cell_growth_size_pixels)
        # add it into the dictionary
        clabels_dict[site] = cell 
        dapi_dict[site] = img

    return dapi_dict, l_adi, adi_dict, nlabels_dict, clabels_dict, da_dict

def gate_by_bounds(leftboundx, rightboundx, nlabels, clabels, d, da):
    # intakes the leftbound, rightbound gates, nuclei labels matrix, celllabels matrix and dictionary(ie what you are gating on: avg dapi intensity)
    # returns a newlabels array containing gated nuclei labels (removing labels belonging to nuclei that are not within the specified bounds)
    # and returns a newlabels array containing gated cyto labels (removing cytoplasms grown from nuclei that are not within the specified bounds)
    gated_nlabels = nlabels.copy() #make a copy of label matrix based on old labels matrix 
    gated_clabels = clabels.copy()
    delabels = [] # initiate a list to keep the labels that have features that fall out of desired range

    # here d will be each adi (average intensity dictionary)
    for i in d:
        if d[i] <= leftboundx or d[i] >= rightboundx:
            delabels.append(i) # add the label to the list of labels to drop from labels matrix 
        if da[i] < 20 and i not in delabels: #gating out ultra tiny small speckles
            delabels.append(i)
    # now loop through labels, if label is in list of delabels, replace label with 0 (ie make it background)
    rows = np.shape(gated_nlabels)[0] 
    cols = np.shape(gated_nlabels)[1]
    for y in range(rows):
        for x in range(cols):        
            if nlabels[y][x] in delabels:
                gated_nlabels[y][x] = 0  #rewrite it to background rather than a label 
            if clabels[y][x] in delabels: #do the same for the cell labels 
                gated_clabels[y][x] = 0

    return gated_nlabels, gated_clabels


def get_all_labels(path, gatebyfixedvalues):
    get_all_labels_start_time = time.time()
    #takes in a folder of images, 
    # outputs a dictionary of gated nuclei labels, gated cell labels, cyto labels and ungated cell labels. 
    dapi_dict, l_adi, adi_dict, nlabels_dict, clabels_dict, da_dict = initialize_all_labels(path, cell_growth_size_pixels)

    getalllabel_initialize_all_labels_time = time.time()
    print("Time to initialize all labels", getalllabel_initialize_all_labels_time - get_all_labels_start_time)
    x, y, x_norm, y_norm, leftboundx, rightboundx =get_kde_norm_bounds(l_adi)

    #if you choose to gate by fixed values, leftboundx, rightboundx will be rewritten to those fixed values
    if gatebyfixedvalues != 'No':
        leftboundx = gatebyfixedvalues[0]
        rightboundx = gatebyfixedvalues[1]

    # plot histogram of avg dapi int per nuclei, showing you the gates and the probability density curve.
    plot_histograms(l_adi, x, y, x_norm, y_norm, leftboundx, rightboundx)
    get_all_labels_gate_bounds_time = time.time()
    print("Time to get gate bounds", get_all_labels_gate_bounds_time - getalllabel_initialize_all_labels_time)

    #now loop through dictionary nlabels_dict
    gated_nlabels_dict = {} #initiate a new dictionary to hold the gated nuclear label matrices at each site
    gated_cytolabels_dict = {} # initiate a new dictionary to hold gated cytoplasmic image label matrices at each site 
    for site in list(nlabels_dict.keys()):
        # pull out the nuclei labels and cell labels for each image. 
        nlabels = nlabels_dict[site]
        
        # grow nuclei by 1 micron (= ~ 3 pixels)
        # so that when you subtract from cell to give cyto, it avoids taking the nuclear boundary pixel count. 
        nlabels_exp = skimage.segmentation.expand_labels(nlabels, distance = 3)
       
        # extract the relevant labels and dictionaries for that site
        clabels = clabels_dict[site]
        adi = adi_dict[site]
        da = da_dict[site]

        # make cyto mask first, and subtract out expanded nuclei to avoid nuclear envelope boundary
        cytolabs = clabels - nlabels_exp
        #now call the gating function on that site image

        gated_nlabels, gated_cytolabels = gate_by_bounds(leftboundx, rightboundx, nlabels, cytolabs, adi, da)
        
        #shrink the nuclei labels circumferentially by 1um (~3 pixels)
        nn1 = shrink_nucleus(gated_nlabels)
        nn2 = shrink_nucleus(nn1)
        nn3 = shrink_nucleus(nn2)
            
        gated_nlabels = nn3
    
        # now check for exceptions - somehow maybe a label is in gated_nlabels but not in gated_cytolabels or vice versa
        # if so, then drop that corresponding label
        common_labels = np.intersect1d(gated_nlabels, gated_cytolabels)
        # now rewrite the labels, any that is not present in both is set to 0 (background). Basically dropping that label
        gated_nlabels[~ np.isin(gated_nlabels, common_labels)] = 0
        gated_cytolabels[~ np.isin(gated_cytolabels, common_labels)] = 0  


        # assign new checked and gated label matrices to the dictionary of gated labels for each site image. 
        gated_nlabels_dict[site] = gated_nlabels    
        gated_cytolabels_dict[site] = gated_cytolabels

        images = [dapi_dict[site], nlabels, gated_nlabels, clabels, cytolabs, gated_cytolabels]
        titles = ['Original DAPI image', 'Nuclei mask', 'Gated shrunken nuclei mask', 'Cells', 'Boundary excluded cyto mask', 'Gated cyto mask']
        name = site + '_segmentation'
        print(name) # so that you see the progress of the code running 

        plotimages(images, titles, name)

    return gated_nlabels_dict, gated_cytolabels_dict, clabels_dict

def summing_image_intensities(
        path,
        file,
        gated_nlabels_dict,
        gated_cytolabels_dict,
        ALL_percentage_intperarea_nuc_perlabel,
        ALL_intperarea_nuc_perlabel,
        ALL_intperarea_cyto_perlabel,
        channel
    ):
    # now that you have the gated labels for nuclear and cyto, 
    # in this function we extract the pixel values corresponding to those labels in the fluorophore images
    # and tabulate per nuclei values for each fluorophore in each space (nuc, cyto)
    filepath = os.path.join(path, file)
    img = io.imread(filepath, as_gray=True)
    site = str(file[-13:-7]) #obtain site title for each well image eg "1-4_w1"
    well = site[:-2]


    #now grab the labels array for the nucleus and cytoplasmic gated masks for this site 
    site_nuc_labels = gated_nlabels_dict[site]
    site_cyto_labels = gated_cytolabels_dict[site]
    
    rows = np.shape(site_nuc_labels)[0]
    cols = np.shape(site_nuc_labels)[1]

    # now we will sum the pixels for each nucleus in this site, and store in a dictionary
    # initiate temporary dictionaries to hold counts of pixels  - these are reinitailised for each image. 
    nuc_pixint_perlabel_dict = {} # note these only hold the sums for each site : eg 1-1
    cyto_pixint_perlabel_dict = {}
    total_pixint_perlabel_dict = {}

    # initiate dictionaries to hold per label nuclear and cyto areas
    area_nuc_perlabel_dict = {}
    area_cyto_perlabel_dict = {}
    prefix = site + '-'
    for y in range(rows):
        for x in range(cols):
            if site_nuc_labels[y][x] != 0: # if its labelled as a coord belonging to a nucleus
                currlabel = site_nuc_labels[y][x] #extract the current label being passed, eg "1"
                newlabel = prefix + str(currlabel) #new label has site infront of the nuclear label, (for future concatenation into large overarching dictionaries)
             
                if newlabel not in nuc_pixint_perlabel_dict.keys():
                    nuc_pixint_perlabel_dict[newlabel] = 0 # initiate sum if label not already in dictionaries 

                if newlabel not in area_nuc_perlabel_dict.keys():
                    area_nuc_perlabel_dict[newlabel] = 0

                # add that pixel intensity at that point to the labels value in the dictionary 
                nuc_pixint_perlabel_dict[newlabel]+= img[y][x] 
                # add that point the area of that label 
                area_nuc_perlabel_dict[newlabel] += 1  
                # note this is different from dareas tabulated earlier, because those are on unshrunken, ungated nuclei masks

            if site_cyto_labels[y][x] != 0: # if its labelled as a coord belonging to a cytoplasm
                currlabel = site_cyto_labels[y][x] #extract the current label being passed, eg "1"
                newlabel = prefix + str(currlabel)     
                # initiate sum if label not already in dictionaries 
                if newlabel not in cyto_pixint_perlabel_dict.keys():
                    cyto_pixint_perlabel_dict[newlabel] = 0 
                if newlabel not in area_cyto_perlabel_dict.keys():
                    area_cyto_perlabel_dict[newlabel] = 0
                
                cyto_pixint_perlabel_dict[newlabel]+= img[y][x] # add that pixel intensity at that point to the labels value in the dictionary
                area_cyto_perlabel_dict[newlabel] += 1 

    # After looping through the whole image, you now have four dictionaries: one that points to the sum of nuclei color intensity per nucleus label eg: nuc_pixint_perlabel_dict[1]=36
    #  one that points to the sum of cytoplasm color intensity per cyto label (which should correspond to its respective nuclei label if they belong in the same cell
    # one that points to the area of the nucleus for each label (eg area_nuccell_dict[1] = 12)
    # and one that points to the area of the cytoplasm for each label (eg area_cyto_perlabel_dict[1] = 20)

    # now we first calculate average int per area for nuclear and cyto for each label
    # generate new dictionaries to store these new numbers
    intperarea_nuc_perlabel = {key: nuc_pixint_perlabel_dict[key]/area_nuc_perlabel_dict[key] for key in nuc_pixint_perlabel_dict}
    intperarea_cyto_perlabel = {key: cyto_pixint_perlabel_dict[key]/area_cyto_perlabel_dict[key] for key in cyto_pixint_perlabel_dict}
    # these newly generated dictionaries now give eg:  intperarea_nuc_perlabel[1] = 3 (because took 36/12)

    # next, create a dictionary to tabulate "total" avg concentration wihtin each mask in the cell. Key is site_label
    intperarea_cytoplusnuc_perlabel = {}
    percentage_intperarea_nuc_perlabel = {}
    for key in intperarea_nuc_perlabel:
        intperarea_cytoplusnuc_perlabel[key] = intperarea_cyto_perlabel[key] + intperarea_nuc_perlabel[key] # sum the average intensities for nuc and cyto PER LABEL
        # next, we want to create a dictionary of % nuclear for each cell. Key is label 
        percentage_intperarea_nuc_perlabel[key] = intperarea_nuc_perlabel[key] / intperarea_cytoplusnuc_perlabel[key]

    # the dictionaries with prefix "ALL_" are initiated and defined outside of this function (in main)
    # so they are not rewritten to 0 everytime this function runs
    # and hence they hold a running cumulative count as the main iterates through the folder of images. 

    #previously, we renamed the label with their site prefix (eg 1-1-1: site 1-1, label1 in format [well]-[view]-[label]. site = [well]-[view])
    # so that we can concatenate dictionaries individual of different sites in the same well 
    # and append them to the same well which will be the key in the large overarching dict with prefix ALL_
    # make a massive dictionary that will hold % nuc per label per well across sites. The key will be well 
    # so eg ALL_percentage_intperarea_nuc_perlabel[1] : 1-1-1  : 0.5 , 1-2-1: 0.7 etc 
    if well not in ALL_percentage_intperarea_nuc_perlabel: 
        ALL_percentage_intperarea_nuc_perlabel[well] = {} # initiate a dictionary to hold the % nuc counts of that site. Each value is for one cell (or label within each site)
    ALL_percentage_intperarea_nuc_perlabel[well].update(percentage_intperarea_nuc_perlabel) # concatenate the dictionaries from this particular site to other sites in that same well
    
    # also make a massive dictionary that holds the raw avg nuc and avg cyto counts per label for each well. The key will be well 
    if well not in ALL_intperarea_nuc_perlabel:
        ALL_intperarea_nuc_perlabel[well] = {} # initiate dictionary
    if well not in  ALL_intperarea_cyto_perlabel:
        ALL_intperarea_cyto_perlabel[well] = {} #initiate dictionary
    ALL_intperarea_nuc_perlabel[well].update(intperarea_nuc_perlabel) # concatenate the label dicitionaries
    ALL_intperarea_cyto_perlabel[well].update(intperarea_cyto_perlabel)

    # plot to check overlay of masks and fluorophore images. these are output in the overlay_and_segmentation folder
    masks = [site_nuc_labels, site_cyto_labels]
    titles = ['Gated nuc mask overlay_'+ channel, 'Gated cyto mask overlay_'+ channel]
    plot_overlay(img, masks, titles, site, channel)

def main():

    start_time = time.time()

    gated_nlabels_dict, gated_cytolabels_dict, clabels_dict = get_all_labels(path, gatebyfixedvalues)
    get_all_labels_time = time.time()
    print("get_all_labels:", get_all_labels_time - start_time)
    print('Dictionaries generated, now summing intensities')
    
    # here we define and initiate the large overarching dictionaries to pass into summing_image_intensities
    GFP_ALL_percentage_intperarea_nuc_perlabel = {}
    GFP_ALL_intperarea_nuc_perlabel = {}
    GFP_ALL_intperarea_cyto_perlabel = {}
    RFP_ALL_percentage_intperarea_nuc_perlabel = {}
    RFP_ALL_intperarea_nuc_perlabel = {}
    RFP_ALL_intperarea_cyto_perlabel = {}

    summing_intensities_start_time = time.time()
    for file in os.listdir(path):
        if not file.endswith('w2.TIF') and not file.endswith('w3.TIF'):
            continue

        if file.endswith('w2.TIF'): 
            # for gfp images
            summing_image_intensities(
                path,
                file,
                gated_nlabels_dict,
                gated_cytolabels_dict,
                GFP_ALL_percentage_intperarea_nuc_perlabel,
                GFP_ALL_intperarea_nuc_perlabel,
                GFP_ALL_intperarea_cyto_perlabel,
                'GFP'
            )
        elif file.endswith('w3.TIF'):
            # for rfp images 
            summing_image_intensities(
                path,
                file,
                gated_nlabels_dict,
                gated_cytolabels_dict,
                RFP_ALL_percentage_intperarea_nuc_perlabel,
                RFP_ALL_intperarea_nuc_perlabel,
                RFP_ALL_intperarea_cyto_perlabel,
                'RFP'
            )

    summing_intensities_end_time = time.time()
    print("Summing intensities time:", summing_intensities_end_time - summing_intensities_start_time)

    # now write results to dataframes
    dataframewriting_starttime = time.time()
    print('Writing to dataframe')
    
    # extracts just the values from the large overarching dictionary, so that you just have gfp_percent_nuc[1] = [0.5, 0.7, ...]
    # ie you drop the assignment to specific labels, and just tabulate their values in a column. 
    gfp_percentnuc = {key: list(value.values()) for key, value in sorted(GFP_ALL_percentage_intperarea_nuc_perlabel.items())}
    gfp_intperarea_nuc = {key: list(value.values()) for key, value in sorted(GFP_ALL_intperarea_nuc_perlabel.items())}
    gfp_intperarea_cyto = {key: list(value.values()) for key, value in sorted(GFP_ALL_intperarea_cyto_perlabel.items())}

    rfp_percentnuc = {key: list(value.values()) for key, value in sorted(RFP_ALL_percentage_intperarea_nuc_perlabel.items())}
    rfp_intperarea_nuc = {key: list(value.values()) for key, value in sorted(RFP_ALL_intperarea_nuc_perlabel.items())}
    rfp_intperarea_cyto = {key: list(value.values()) for key, value in sorted(RFP_ALL_intperarea_cyto_perlabel.items())}
   
    # convert these new dictionaries to dataframes, with well number as column name, values in the column
    gfp_percentnuc_df = pd.DataFrame.from_dict(gfp_percentnuc, orient= 'index').transpose()
    gfp_intperarea_nuc_df = pd.DataFrame.from_dict(gfp_intperarea_nuc, orient= 'index').transpose()
    gfp_intperarea_cyto_df = pd.DataFrame.from_dict(gfp_intperarea_cyto, orient= 'index').transpose()

    rfp_percentnuc_df = pd.DataFrame.from_dict(rfp_percentnuc, orient = 'index').transpose()
    rfp_intperarea_nuc_df = pd.DataFrame.from_dict(rfp_intperarea_nuc, orient= 'index').transpose()
    rfp_intperarea_cyto_df = pd.DataFrame.from_dict(rfp_intperarea_cyto, orient= 'index').transpose()
    
    # make a new dictionary where you take the average of all the individual label values for each well
    # so dict[well] = 20 which is in an average of the values in the list. 
    averaged_gfp_percentnuc_dict = {key: value.mean() for key, value in sorted(gfp_percentnuc_df.items())}
    averaged_gfp_intperarea_nuc_dict = {key: value.mean() for key, value in sorted(gfp_intperarea_nuc_df.items())}
    averaged_gfp_intperarea_cyto_dict = {key: value.mean() for key, value in sorted(gfp_intperarea_cyto_df.items())}

    averaged_rfp_percentnuc_dict = {key: value.mean() for key, value in sorted(rfp_percentnuc_df.items())}
    averaged_rfp_intperarea_nuc_dict = {key: value.mean() for key, value in sorted(rfp_intperarea_nuc_df.items())}
    averaged_rfp_intperarea_cyto_dict = {key: value.mean() for key, value in sorted(rfp_intperarea_cyto_df.items())}

    # make a dictionary for number of cells analysed in each well. Key is well, value is ncells analysed for that well
    ncellsanalysed_dict = {key: value.count() for key, value in sorted(rfp_percentnuc_df.items())}

    # create dataframe for summary data. Row is well, column is data type averaged for each well(%, avgnuc, avgcyto etc)
    gfp_summary_df = pd.DataFrame.from_dict(ncellsanalysed_dict, orient = 'index', columns = ['no of cells'])
    rfp_summary_df = pd.DataFrame.from_dict(ncellsanalysed_dict, orient = 'index', columns = ['no of cells'])
    
    gfp_summary_df.reset_index(inplace=True)
    rfp_summary_df.reset_index(inplace=True)

    # initialise the columns
    # code takes in every value in the rfp_summary_df["index"] column which corresponds to the well number
    # and will look for that value in your dictionary
    # if value is found in the dict, corresponding value will populate the column
    gfp_summary_df['GFP avg percent nuclear'] = gfp_summary_df["index"].apply(lambda x: averaged_gfp_percentnuc_dict.get(x))
    gfp_summary_df['GFP avg nuclear intensity'] = gfp_summary_df["index"].apply(lambda x: averaged_gfp_intperarea_nuc_dict.get(x))
    gfp_summary_df['GFP avg cyto intensity'] = gfp_summary_df["index"].apply(lambda x: averaged_gfp_intperarea_cyto_dict.get(x))

    rfp_summary_df['RFP avg percent nuclear'] = rfp_summary_df["index"].apply(lambda x: averaged_rfp_percentnuc_dict.get(x))
    rfp_summary_df['RFP avg nuclear intensity'] = rfp_summary_df["index"].apply(lambda x: averaged_rfp_intperarea_nuc_dict.get(x))
    rfp_summary_df['RFP avg cyto intensity'] = rfp_summary_df["index"].apply(lambda x: averaged_rfp_intperarea_cyto_dict.get(x))
  
    executionTime = (time.time() - startTime)
    dataframewriting_endtime = time.time()
    # finally, write settings to excel sheet 
    settings_df = pd.DataFrame({'Settings': ['segmentation_footprint', 'segmentation_minimum_distance', 'Cell growth pixels', 'gatebyfixedvalues'], 'Pixels': [segmentation_footprint, cell_growth_size_pixels, segmentation_minimum_distance, 'avg DAPI nucint bounds' + str(gatebyfixedvalues)]})
    efficiency_df = pd.DataFrame({'Times': ['total execution time', 'time for generating all labels', 'Time for summing intensities', 'Time to write to dataframe'], 'Values': [executionTime, get_all_labels_time - start_time, summing_intensities_end_time - summing_intensities_start_time, dataframewriting_endtime - dataframewriting_starttime]})
    
    print('Execution time in seconds: ' + str(executionTime))

    # now, we write the dataframes generated to excel sheets 

    expt = os.path.basename(path1)
    with pd.ExcelWriter(path1 + "/%s_Results.xlsx" %expt) as writer:
        gfp_summary_df.to_excel(writer, sheet_name='Summary GFP', index=False)
        rfp_summary_df.to_excel(writer, sheet_name='Summary RFP', index=False)

        rfp_percentnuc_df.to_excel(writer, sheet_name='RFP_percentnuc_perlabel', index=False)
        rfp_intperarea_nuc_df.to_excel(writer, sheet_name='RFP_avgnucint_perlabel', index=False)
        rfp_intperarea_cyto_df.to_excel(writer, sheet_name='RFP_avgcytoint_perlabel', index=False)
        
        gfp_percentnuc_df.to_excel(writer, sheet_name='GFP_percentnuc_perlabel', index=False)
        gfp_intperarea_nuc_df.to_excel(writer, sheet_name='GFP_averagenucint_perlabel', index=False)
        gfp_intperarea_cyto_df.to_excel(writer, sheet_name='GFP_averagecytoint_perlabel', index=False)

        settings_df.to_excel(writer, sheet_name='Settings', index=False)
        efficiency_df.to_excel(writer, sheet_name='Settings', index=False, startrow=6)
    
    # here, we create another excel file, whereby we tabulate the gfp and rfp individual values assigned to each individual label eg 1-1-5
    # for each well, create a new sheet. 
    with pd.ExcelWriter(path1 + "/%s_GFP-RFP_percellcombinations.xlsx" %expt) as writer:
        # all your dictionaries should have similar keys (which are wells), and similar inner keys (site-labels)
        common_keys = sorted(GFP_ALL_percentage_intperarea_nuc_perlabel.keys()) # pick one dictionary to extract the keys 
        for key in common_keys:  # for each well 
            # following code sorts the inner keys for the current well in the dictionary GFP_ALL_percentage_intperarea_nuc_perlabel
            # then the list comprehension extracts the corresponding values for those inner keys in a different dictionary, and populates them in a list
            sorted_gfp_percentnuc_values = [GFP_ALL_percentage_intperarea_nuc_perlabel[key].get(inner_key, '') for inner_key in sorted(GFP_ALL_percentage_intperarea_nuc_perlabel[key])]
            sorted_rfp_percentnuc_values = [RFP_ALL_percentage_intperarea_nuc_perlabel[key].get(inner_key, '') for inner_key in sorted(GFP_ALL_percentage_intperarea_nuc_perlabel[key])]
            sorted_gfp_avgnucint_values = [GFP_ALL_intperarea_nuc_perlabel[key].get(inner_key, '') for inner_key in sorted(GFP_ALL_percentage_intperarea_nuc_perlabel[key])]
            sorted_rfp_avgnucint_values = [RFP_ALL_intperarea_nuc_perlabel[key].get(inner_key, '') for inner_key in sorted(GFP_ALL_percentage_intperarea_nuc_perlabel[key])]
            sorted_gfp_avgcytoint_values = [GFP_ALL_intperarea_cyto_perlabel[key].get(inner_key, '') for inner_key in sorted(GFP_ALL_percentage_intperarea_nuc_perlabel[key])]
            sorted_rfp_avgcytoint_values = [RFP_ALL_intperarea_cyto_perlabel[key].get(inner_key, '') for inner_key in sorted(GFP_ALL_percentage_intperarea_nuc_perlabel[key])]
            df_combined = pd.DataFrame({
                # these lists then populate the columns 
                'GFP_avgcytoint' : sorted_gfp_avgcytoint_values,
                'RFP_avgcytoint' : sorted_rfp_avgcytoint_values,
                'GFP_avgnucint' : sorted_gfp_avgnucint_values,
                'RFP_avgnucint' : sorted_rfp_avgnucint_values,
                'GFP_percentnuc': sorted_gfp_percentnuc_values,
                'RFP_percentnuc': sorted_rfp_percentnuc_values
                
            }, index=sorted(GFP_ALL_percentage_intperarea_nuc_perlabel[key])) # set the index to the inner keys that have been used to extract the values 
            
            sheet_name = "well_" + key  # Prepend "well" to the sheet name
            df_combined.to_excel(writer, sheet_name=sheet_name)

# This boiler plate invokes the main() function when the script is run in
# python.sn
if __name__ == '__main__': 
  main()
