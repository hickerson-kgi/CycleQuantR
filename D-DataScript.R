#<><><><><><><><><><><><><><><><><><><><><><>
#------------------------------------
# CREATE A DATAFRAME OF NUCLEI INFORMATION
# extracted from images and masks.
#
# Data extracted currently includes color, size in pixels,
# relative intensity, and position.
#
# At this point, no information linking nuclei across colors / frames
# is assessed or stored. Further use of the masks will no longer be required.
#------------------------------------
#<><><><><><><><><><><><><><><><><><><><><><>

# declare EBImage library for processing images
library(EBImage)

# run B-DataSettings.R script to pull-in variable settings
source("./B-DataSettings.R") 

#----------------------------------------------
# initialize a data frame that will contain the extracted nuclei information
#----------------------------------------------
nuclei = data.frame(
  frame = character(),      # filename of the frame
  color = character(),      # character identifying the color of the nuclei (B=Blue, G=Green, R=Red)
  size = integer(),         # size of nuclei in number of pixels
  intensity = double(),     # average intensity (brightness)
  x_pos = double(),         # center of x-position
  y_pos = double(),         # center of y-position
  stringsAsFactors = FALSE)

#----------------------------------------------
# step through all RED images and masks
#----------------------------------------------
if(R_toggle) {
  print("Analyzing Nuclei of Red Images")
  
  for (i in R_images) {
    
    # Load original and mask images - convert to grayscale
    #----------------------------------------------
    R_orig = readImage(paste(dir_Rraw, i, sep=""))[x_crop, y_crop,] # Load and crop raw image
    img = channel(R_orig, "red")                                    # Convert red channel to grayscale
    bg = gblur(img, 100)                                            # Broad Gaussian Blur
    img = img-bg
    img = img-min(img)
    img = img/max(img)
    
    R_mask = readImage(paste(dir_Rprocessed, i, sep=""))            # Load paired mask image
    mask = channel(R_mask, "red")                                   # Convert red channel to grayscale
    
    
    # identify all unique nuclei in image
    #----------------------------------------------
    dist = distmap(mask)
    if(label_select) {label = watershed(dist, 2)}
    else {label = bwlabel(dist)}
    
    # step through all nuclei in image
    #----------------------------------------------
    for (j in 1:max(label)) {
      
      pixels = which(label==j, arr.ind=TRUE) # identify coordinates of all pixels that match this label
      
      size = length(pixels)                  # size of nuclei
      intensity = sum(img[pixels])/size      # average of all intensity values within the nuclei
      x_pos = mean(pixels[,'row'])           # mean x coordinate of nuclei
      y_pos = mean(pixels[,'col'])           # mean y coordinate of nuclei
      
      nuclei[nrow(nuclei)+1,] = c(i, "R", size, intensity, x_pos, y_pos)
      print(paste("cell", j, "of", max(label), "in red image", i, sep=" "))
    }
  }
}
#----------------------------------------------
# step through all GREEN images and masks
#----------------------------------------------
if(G_toggle) {
  print("Analyzing Nuclei of Green Images")
  
  for (i in G_images) {
    
    # Load original and mask images - convert to grayscale
    #----------------------------------------------
    G_orig = readImage(paste(dir_Graw, i, sep=""))[x_crop, y_crop,] # Load and crop raw image
    img = channel(G_orig, "green")                                  # Convert green channel to grayscale
    bg = gblur(img, 100)                                            # Broad Gaussian Blur
    img = img-bg
    img = img-min(img)
    img = img/max(img)
    
    G_mask = readImage(paste(dir_Gprocessed, i, sep=""))            # Load paired mask image
    mask = channel(G_mask, "green")                                 # Convert green channel to grayscale
  
  
    # identify all unique nuclei in image
    #----------------------------------------------
    dist = distmap(mask)
    if(label_select) {label = watershed(dist, 2)}
    else {label = bwlabel(dist)}
    
    # step through all nuclei in image
    #----------------------------------------------
    for (j in 1:max(label)) {
   
      pixels = which(label==j, arr.ind=TRUE) # identify coordinates of all pixels that match this label
  
      size = length(pixels)                  # size of nuclei
      intensity = sum(img[pixels])/size      # average of all intensity values within the nuclei
      x_pos = mean(pixels[,'row'])           # mean x coordinate of nuclei
      y_pos = mean(pixels[,'col'])           # mean y coordinate of nuclei
      
      nuclei[nrow(nuclei)+1,] = c(i, "G", size, intensity, x_pos, y_pos)
      print(paste("cell", j, "of", max(label), "in green image", i, sep=" "))
    }
  }
}

#------------------------------------
# step through all BLUE images and masks
#------------------------------------
if(B_toggle) {
  print("Analyzing Nuclei of Blue Images")
  
  for (i in B_images) {
    
    # Load original and mask images - convert to grayscale
    #----------------------------------------------
    B_orig = readImage(paste(dir_Braw, i, sep=""))[x_crop, y_crop,] # Load and crop raw image
    img = channel(B_orig, "blue")                                   # Convert blue channel to grayscale
    bg = gblur(img, 100)                                            # Broad Gaussian Blur
    img = img-bg
    img = img-min(img)
    img = img/max(img)
    
    B_mask = readImage(paste(dir_Bprocessed, i, sep=""))            # Load paired mask image
    mask = channel(B_mask, "blue")                                  # Convert blue channel to grayscale
    
    
    # identify all unique nuclei in image
    #----------------------------------------------
    dist = distmap(mask)
    if(label_select) {label = watershed(dist, 2)}
    else {label = bwlabel(dist)}
    
    # step through all nuclei in image
    #----------------------------------------------
    for (j in 1:max(label)) {
      
      pixels = which(label==j, arr.ind=TRUE) # identify coordinates of all pixels that match this label
      
      size = length(pixels)                  # size of nuclei
      intensity = sum(img[pixels])/size      # average of all intensity values within the nuclei
      x_pos = mean(pixels[,'row'])           # mean x coordinate of nuclei
      y_pos = mean(pixels[,'col'])           # mean y coordinate of nuclei
      
      nuclei[nrow(nuclei)+1,] = c(i, "B", size, intensity, x_pos, y_pos)
      print(paste("cell", j, "of", max(label), "in blue image", i, sep=" "))
    }
  }
}

# SAVE results
save(nuclei, file=paste(dir_data, "nuclei_data.data", sep=""))
write.table(nuclei, file=paste(dir_data, "nuclei_data.csv", sep=""),row.names = FALSE, col.names = TRUE, sep = ", ")
