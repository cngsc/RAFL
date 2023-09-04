import sys, os
import numpy as np
import nd2
from skimage import io
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore", lineno=44)

# this command gets the folder that you have copied this code file into
path1 = os.getcwd()

#this is the path to your folder of input images, which should be within the folder you copied this code file into
input_folder_path= path1 + '/raw_images/'

#this is an output folder that is generated in this code to store your converted tif files
output_folder_path = os.path.join(path1, 'convertedtifs')

os.makedirs(output_folder_path, exist_ok=True) # make output folder if it doesnt exist

def openfile(nd2file, output_folder_path):
    allchannels_array = nd2.imread(nd2file)


    num_channels = allchannels_array.shape[0]

    # if you run the following code, channel 0 is dapi, then gfp, rfp, cy5, td
    # with nd2.ND2File(nd2file) as ndfile:
    #     print(ndfile.metadata)
    
    # now loop through the indexes of the arraays 
    for i in range(num_channels):

        #obtain the image for that particular channel
        channel_data = allchannels_array[i]
        suffix = f"_w{i + 1}"

        # split string so that you obtain the image name (ie the last str after the path backslash, minus the extension of .nd2)
        filename=nd2file.rsplit('/', 1)[-1][:-4]
        
        outputfilename = f"{filename}{suffix}.TIF"
        output_path = os.path.join(output_folder_path, outputfilename)
        
        #Save the channel data as a .TIF image
        io.imsave(output_path, channel_data)


# This boiler plate invokes the main() function when the script is run in
# python.sn

def main():

    for file in os.listdir(input_folder_path):
        if file.endswith(".nd2"):
            nd2_file = os.path.join(input_folder_path, file)
            openfile(nd2_file, output_folder_path)
if __name__ == '__main__': 
  main()
