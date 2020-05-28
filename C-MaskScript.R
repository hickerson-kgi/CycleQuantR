#------------------------------------
# Script for Batch Processing Images
# 
# Raw images located in each of the RAW folders
# is processed to generate a mask identifying
# individual nuclei using color specific data settings
#------------------------------------

# declare EBImage library for processing images
library(EBImage)

# run B-DataSettings.R script to pull-in variable settings
source("./B-DataSettings.R") 

#------------------------------------------------------
# FUNCTION FOR GENERATING COLOR MASKS
#------------------------------------------------------

Mask <- function(img, blur, thresh, offset) {
  img = gblur(img, blur)                      # Gaussian Blur
  img = thresh(img, thresh, thresh, offset) # Adaptive Threshold
}

#------------------------------------------------------
# SCRIPT TO BATCH PROCESS EACH COLOR AND GENERATE MASKS
#------------------------------------------------------

if (R_toggle)  {
  # open, crop, and process Red images
  #----------------------------------------------
  print("Generating Masks of Red Images")
  
  for (i in R_images) {
    R_orig = readImage(paste(dir_Rraw, i, sep=""))        # open file
    R_orig = R_orig[x_crop, y_crop, ]                     # crop image according to data settings
    R_mask = Mask(R_orig, R1, R2, R3)                     # generate mask
    writeImage(R_mask, paste(dir_Rprocessed, i, sep=""))  # save mask image 
    print(paste("Red Mask", i))                           # update console on status
  }
}

if (G_toggle)  {
  # open, crop, and process Green images
  #----------------------------------------------
  print("Generating Masks of Green Images")
  
  for (i in G_images) {
    G_orig = readImage(paste(dir_Graw, i, sep=""))        # open file
    G_orig = G_orig[x_crop, y_crop, ]                     # crop image according to data settings
    G_mask = Mask(G_orig, G1, G2, G3)                     # generate mask
    writeImage(G_mask, paste(dir_Gprocessed, i, sep=""))  # save mask image 
    print(paste("Green Mask", i))                         # update console on status
  }
}

if (B_toggle)  {
  # open, crop, and process Blue images
  #----------------------------------------------
  print("Generating Masks of Blue Images")
  
  for (i in B_images) {
    B_orig = readImage(paste(dir_Braw, i, sep=""))        # open file
    B_orig = B_orig[x_crop, y_crop, ]                     # crop image according to data settings
    B_mask = Mask(B_orig, B1, B2, B3)                     # generate mask
    writeImage(B_mask, paste(dir_Bprocessed, i, sep=""))  # save mask image 
    print(paste("Blue Mask", i))                          # update console on status
  }
}

#------------------------------------------------------
# SCRIPT TO GENERATE OVERLAY IMAGES BASED ON MASKED AND PHASE IMAGE
#------------------------------------------------------

print("Generating Overlay of Masks and Phase Images")

for (i in W_images) {
  
  # open phase image and crop to mask dimensions  
  overlay = readImage(paste(dir_Wraw, i, sep=""))[x_crop, y_crop,]

  # load red mask and add to overlay
  if (R_toggle)
  {
    R_mask = readImage(paste(dir_Rprocessed, i, sep=""))
    R_mask = channel(R_mask, 'red')
    overlay = paintObjects(R_mask, overlay, col=c(NA, "red"), opac=c(1,0.3))

    dist = distmap(R_mask)
    if(label_select) {label = watershed(dist, 2)}
    else {label = bwlabel(dist)}
    overlay = paintObjects(label, overlay, col="red")
  }
  
  # load green mask and add to overlay
  if (G_toggle)
  {
    G_mask = readImage(paste(dir_Gprocessed, i, sep=""))
    G_mask = channel(G_mask, 'green')
    overlay = paintObjects(G_mask, overlay, col=c(NA, "green"), opac=c(1,0.3))

    dist = distmap(G_mask)
    if(label_select) {label = watershed(dist, 2)}
    else {label = bwlabel(dist)}
    overlay = paintObjects(label, overlay, col="green")
  }
  
  # load blue mask and add to overlay
  if (B_toggle)
  {
    B_mask = readImage(paste(dir_Bprocessed, i, sep=""))
    B_mask = channel(B_mask, 'blue')
    overlay = paintObjects(B_mask, overlay, col=c(NA, "blue"), opac=c(1,0.3))

    dist = distmap(B_mask)
    if(label_select) {label = watershed(dist, 2)}
    else {label = bwlabel(dist)}
    overlay = paintObjects(label, overlay, col="blue")
  }
  
  # save overlay image
  writeImage(overlay, paste(dir_Wprocessed, i,sep=""), type="tif")
  print(paste("Overlay", i))
}