#------------------------------------
# Data Settings
#
# Consolidated set of variables that will
# be used throughout the code for processing.
# These may need to be adjusted according to your data
#------------------------------------

#------------------------------------
# Directory and filename variables
#------------------------------------
# Raw File Directories
dir_Rraw = "./RED/RAW/"
dir_Graw = "./GREEN/RAW/"
dir_Braw = "./BLUE/RAW/"
dir_Wraw = "./WHITE/RAW/"

# Processed File Directories
dir_Rprocessed = "./RED/PROCESSED/"
dir_Gprocessed = "./GREEN/PROCESSED/"
dir_Bprocessed = "./BLUE/PROCESSED/"
dir_Wprocessed = "./WHITE/PROCESSED/"
dir_data = "./DATA/"

# Pulls list of files in Raw File Directories
R_images = list.files(dir_Rraw)
G_images = list.files(dir_Graw)
B_images = list.files(dir_Braw)
W_images = list.files(dir_Wraw)

#------------------------------------
# Mask processing variables
#------------------------------------

#Toggle Processing of Specific Colors
R_toggle = TRUE # R = Red
G_toggle = TRUE # G = Green
B_toggle = TRUE # B = Blue

# X and Y Coordinates for image cropping
x_crop = 50:1550
y_crop = 50:1550

#Red Image to Mask Processing Settings
R1 = 7 # Amount of Guassian Blur
R2 = 40 #  Threshold Size (both heigh and width) in Pixels
R3 = 0.019 # Threshold Offset (higher is more selective)

#Green Image to Mask Processing Settings
G1 = 6 # Amount of Guassian Blur
G2 = 40 #  Threshold Size (both heigh and width) in Pixels
G3 = 0.006 # Threshold Offset 

#Blue Image to Mask Processing Settings
B1 = 5 # Amount of Guassian Blur
B2 = 40 #  Threshold Size (both heigh and width) in Pixels
B3 = 0.015 # Threshold Offset  

#------------------------------------
# Data Analysis variables
#------------------------------------

#Select Method for labeling nuclei
label_select = TRUE # FALSE for simple bwlabel() method
                    # TRUE for watershed() method

#Nuclei Size Exclusion Settings
Nuclei_Lower_Size = 200
Nuclei_Upper_Size = 100000

# Select base color as reference for nuclei identification
C = "B" # R = Red, G = Green, B = Blue

# AZ-Timelapse information
T = 10 # Minute interval between frames set 0 for single time analysis


