#------------------------------------
# Folder and File Structure Setup
#
# Creates appropriate folder structure
# that is anticipated by the subsequent code
#------------------------------------

# Create folders for each of the filter colors
dir.create(file.path("RED"))
dir.create(file.path("GREEN"))
dir.create(file.path("BLUE"))
dir.create(file.path("WHITE"))

# Create subfolders intended to contain the original
# raw images that will be processed
dir.create(file.path("RED/RAW"))
dir.create(file.path("GREEN/RAW"))
dir.create(file.path("BLUE/RAW"))
dir.create(file.path("WHITE/RAW"))

# Create subfolders intended to contain the mask images
# after the originals have been processed, and a 
# data folder for storing all quantitative data
dir.create(file.path("RED/PROCESSED"))
dir.create(file.path("GREEN/PROCESSED"))
dir.create(file.path("BLUE/PROCESSED"))
dir.create(file.path("WHITE/PROCESSED"))
dir.create(file.path("DATA"))


#***  PAUSE   ***#
# At this stage, place microscope images in the appropriate
# RAW folders based on their fluorescence color

