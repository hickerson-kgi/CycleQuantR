#----------------------------------------------
# Script for Analyzing and Plotting Nuclei Info
# from Alex Zambon Data
# 
# Written by: Anna Hickerson, 03-12-2018
# updated: 05-04-2020 to accomodate revised data file format
#----------------------------------------------

# Load data file
load("./DATA/nuclei_master.data", envir = parent.frame(), verbose = TRUE)
master$n_frame = as.numeric(substr(master$frame, 0, 4)) # create numerical frame numbers

source("./B-DataSettings.R")

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#------------------------------------------------------------------------------------------------------------------------------------------------
# Calling Functions for Analyzing and Graphing Nuclei Data
#------------------------------------------------------------------------------------------------------------------------------------------------
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# manually designate two cells with IDs c1 and c2 as a single cell with whichever ID is lower
ID_fix(c1, c2)

# identify nuclei ID present in the most number of frames
ID = Duration_max(nuclei)

# Plot the Intensity, Size, and Position of a specific nuclei ID as a function of time 
Plot_ID(nuclei, ID)

# identify nuclei IDs with a duration above a specified number of frames (Ex: 12 frames = 60 min)
IDs = Duration_above(nuclei, 12)

# Plot Sizes of multiple nuclei identified by IDS as a function of time (in color)
Plot_Areas(nuclei, IDs)

# Plot Intensities of multiple nuclei identified by IDS as a function of time (in color)
Plot_Intensities(nuclei, IDs)

# Plot Positions of multiple nuclei identified by IDS as a function of time (in color)
Plot_Positions(nuclei, IDs)

# Horizontal bar graph of Ages of multiple nuclei identified by IDS as a function of time (in color)
Plot_Ages(nuclei, IDs)

# Plot Histogram of all the Nuclei Sizes
Plot_Histogram(nuclei)

# 3D Scatter plot of a single cell with ID in color
Plot_Scatterplot(nuclei, ID)

# Plot Nuclei positions and colors per frame and save as images to externally create movie
Create_Position_Frames()
  
# Plot Nuclei positions and search area on top of blue mask per frame and save as images to externally create movie
Overlay_Search_Area()


#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#------------------------------------------------------------------------------------------------------------------------------------------------
# Functions for Analyzing and Graphing Nuclei Data
#------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------
# manually designate two cells with IDs c1 and c2
# as a single cell with whichever ID is lower
# includes replacing parent.
ID_fix <- function(nuclei, c1, c2)
{
  nuclei$id[nuclei$id == max(c1, c2)] = min(c1, c2)
  nuclei$parent[nuclei$parent == max(c1, c2)] = min(c1, c2)
}

#------------------------------------------------------------------------------------------------------------------------------------------------
# identify nuclei ID present in the most number of frames
Duration_max <- function(nuclei) {
  ID = which.max(table(nuclei$id))
  ID = as.numeric(ID)
}
#------------------------------------------------------------------------------------------------------------------------------------------------
# identify nuclei IDs with a duration above a specified number of frames (Ex: 12 frames = 60 min)
Duration_above <- function(nuclei, n_frames) {
  IDs = table(nuclei$id[-1])
  IDs = IDs[IDs > n_frames]
  IDs = as.numeric(names(IDs))
}
#------------------------------------------------------------------------------------------------------------------------------------------------
# Plot the Size and Position of a specific nuclei ID as a function of time 
Plot_ID <- function(nuclei, ID) {

  x = nuclei$x_pos[nuclei$id==ID]
  y = nuclei$y_pos[nuclei$id==ID]
  n_frame = nuclei$n_frame[nuclei$id==ID]
  
  if(C == "R") {
    colors=rgb(0, nuclei$green[nuclei$id==ID], nuclei$blue[nuclei$id==ID])
    size = nuclei$size_red[nuclei$id==ID]
  }
  else if (C=="G") {
    colors=rgb(nuclei$red[nuclei$id==ID], 0, nuclei$blue[nuclei$id==ID])
    size = nuclei$size_green[nuclei$id==ID]
  }
  else {
    colors=rgb(nuclei$red[nuclei$id==ID], nuclei$green[nuclei$id==ID], 0)
    size = nuclei$size_blue[nuclei$id==ID]
  }
  
  par(mfrow=c(2,1))
  
  # SIZE
  plot(n_frame[1], size[1], col=colors[1], type="l", xlab="time (frames)", ylab="area (pixels)", pch = 2,
       xlim=c(min(n_frame), max(n_frame)), ylim=c(min(size), max(size)), main=paste("Area of Nuclei ID =", ID))
  for (i in 2:length(n_frame)) {lines(n_frame[(i-1):i], size[(i-1):i], col=colors[i])}
  
  # POSITION
  plot(x[1], y[1], col=colors[1], type="l", xlab="x-position", ylab="y-position", pch = 3, asp=1,
       xlim=c(min(x), max(x)), ylim=c(max(y), min(y)), main=paste("Position of Nuclei ID =", ID))
  for (i in 2:length(x)) {lines(x[(i-1):i], y[(i-1):i], col=colors[i])}
}
#------------------------------------------------------------------------------------------------------------------------------------------------
# Plot Sizes of multiple nuclei identified by IDs as a function of time (in color)
Plot_Areas <- function(nuclei, IDs) {
  par(mfrow=c(1,1))

  x = nuclei$x_pos[nuclei$id==IDs[1]]
  y = nuclei$y_pos[nuclei$id==IDs[1]]
  n_frame = nuclei$n_frame[nuclei$id==IDs[1]]
  
  if(C == "R") {
    colors=rgb(0, nuclei$green[nuclei$id==IDs[1]], nuclei$blue[nuclei$id==IDs[1]])
    size = nuclei$size_red[nuclei$id==IDs[1]]
    m_size = max(nuclei$size_red)
  }
  else if (C=="G") {
    colors=rgb(nuclei$red[nuclei$id==IDs[1]], 0, nuclei$blue[nuclei$id==IDs[1]])
    size = nuclei$size_green[nuclei$id==IDs[1]]
    m_size = max(nuclei$size_green)
  }
  else {
    colors=rgb(nuclei$red[nuclei$id==IDs[1]], nuclei$green[nuclei$id==IDs[1]], 0)
    size = nuclei$size_blue[nuclei$id==IDs[1]]
    m_size = max(nuclei$size_blue)
  }
  
  plot(n_frame[1], size[1], col=colors[1], type="l", xlab="time (frames)", ylab="area (pixels)",
       xlim=c(min(nuclei$n_frame), max(nuclei$n_frame)), ylim=c(0, m_size), main="Area of Nuclei")

  for(n in IDs) {

    x = nuclei$x_pos[nuclei$id==n]
    y = nuclei$y_pos[nuclei$id==n]
    n_frame = nuclei$n_frame[nuclei$id==n]
    
    if(C == "R") {
      colors=rgb(0, nuclei$green[nuclei$id==n], nuclei$blue[nuclei$id==n])
      size = nuclei$size_red[nuclei$id==n]
    }
    else if (C=="G") {
      colors=rgb(nuclei$red[nuclei$id==n], 0, nuclei$blue[nuclei$id==n])
      size = nuclei$size_green[nuclei$id==n]
    }
    else {
      colors=rgb(nuclei$red[nuclei$id==n], nuclei$green[nuclei$id==n], 0)
      size = nuclei$size_blue[nuclei$id==n]
    }
    
    for (i in 2:length(n_frame)) {lines(n_frame[(i-1):i], size[(i-1):i], col=colors[i])}
  }
}
#------------------------------------------------------------------------------------------------------------------------------------------------
# Plot Intensities of multiple nuclei identified by IDs as a function of time (in color)
Plot_Intensities <- function(nuclei, color, IDs) {
  par(mfrow=c(1,1))
  
  x = nuclei$x_pos[nuclei$id==IDs[1]]
  y = nuclei$y_pos[nuclei$id==IDs[1]]
  n_frame = nuclei$n_frame[nuclei$id==IDs[1]]
  
  if(C == "R") {
    colors=rgb(0, nuclei$green[nuclei$id==IDs[1]], nuclei$blue[nuclei$id==IDs[1]])
    size = nuclei$size_red[nuclei$id==IDs[1]]
  }
  else if (C=="G") {
    colors=rgb(nuclei$red[nuclei$id==IDs[1]], 0, nuclei$blue[nuclei$id==IDs[1]])
    size = nuclei$size_green[nuclei$id==IDs[1]]
  }
  else {
    colors=rgb(nuclei$red[nuclei$id==IDs[1]], nuclei$green[nuclei$id==IDs[1]], 0)
    size = nuclei$size_blue[nuclei$id==IDs[1]]
  }
  
  if(color == "R") {
    intensity = nuclei$intensity_red[nuclei$id==IDs[1]]
    m_intensity = max(nuclei$intensity_red)
  }
  else if(color == "G") {
    intensity = nuclei$intensity_green[nuclei$id==IDs[1]]
    m_intensity = max(nuclei$intensity_green)
  }
  else {
    intensity = nuclei$intensity_blue[nuclei$id==IDs[1]]
    m_intensity = max(nuclei$intensity_blue)
  }
  
  plot(n_frame[1], intensity[1], col=colors[1], type="l", xlab="time (hr)", ylab="intensity (0-1 scale)",
       xlim=c(min(nuclei$n_frame), max(nuclei$n_frame)), ylim=c(0,m_intensity), main="Intensity of Nuclei")
  
  for(n in IDs) {
    x = nuclei$x_pos[nuclei$id==n]
    y = nuclei$y_pos[nuclei$id==n]
    n_frame = nuclei$n_frame[nuclei$id==n]
    
    if(C == "R") {
      colors=rgb(0, nuclei$green[nuclei$id==n], nuclei$blue[nuclei$id==n])
      size = nuclei$size_red[nuclei$id==n]
    }
    else if (C=="G") {
      colors=rgb(nuclei$red[nuclei$id==n], 0, nuclei$blue[nuclei$id==n])
      size = nuclei$size_green[nuclei$id==n]
    }
    else {
      colors=rgb(nuclei$red[nuclei$id==n], nuclei$green[nuclei$id==n], 0)
      size = nuclei$size_blue[nuclei$id==n]
    }
    
    if(color == "R") {
      intensity = nuclei$intensity_red[nuclei$id==n]
      m_intensity = max(nuclei$intensity_red)
    }
    else if(color == "G") {
      intensity = nuclei$intensity_green[nuclei$id==n]
      m_intensity = max(nuclei$intensity_green)
    }
    else {
      intensity = nuclei$intensity_blue[nuclei$id==n]
      m_intensity = max(nuclei$intensity_blue)
    }
    
    for (i in 2:length(n_frame)) {lines(n_frame[(i-1):i], intensity[(i-1):i], col=colors[i])}
  }
}
#------------------------------------------------------------------------------------------------------------------------------------------------
# Plot Positions of multiple nuclei identified by IDs as a function of time (in color)
Plot_Positions <- function(nuclei, IDs) {
  par(mfrow=c(1,1))
  
  x = nuclei$x_pos[nuclei$id==IDs[1]]
  y = nuclei$y_pos[nuclei$id==IDs[1]]
  n_frame = nuclei$n_frame[nuclei$id==IDs[1]]
  
  if(C == "R") {
    colors=rgb(0, nuclei$green[nuclei$id==IDs[1]], nuclei$blue[nuclei$id==IDs[1]])
    size = nuclei$size_red[nuclei$id==IDs[1]]
  }
  else if (C=="G") {
    colors=rgb(nuclei$red[nuclei$id==IDs[1]], 0, nuclei$blue[nuclei$id==IDs[1]])
    size = nuclei$size_green[nuclei$id==IDs[1]]
  }
  else {
    colors=rgb(nuclei$red[nuclei$id==IDs[1]], nuclei$green[nuclei$id==IDs[1]], 0)
    size = nuclei$size_blue[nuclei$id==IDs[1]]
  }
  
  plot(x[1], y[1], col=colors[1], type="o", xlab="x-position", ylab="y-position",
       xlim=c(min(nuclei$x_pos), max(nuclei$x_pos)), ylim=c(max(nuclei$y_pos), min(nuclei$y_pos)), asp=1, main="Motion of Nuclei")
  
  for(n in IDs) {
    x = nuclei$x_pos[nuclei$id==n]
    y = nuclei$y_pos[nuclei$id==n]
    n_frame = nuclei$n_frame[nuclei$id==n]
    
    if(C == "R") {
      colors=rgb(0, nuclei$green[nuclei$id==n], nuclei$blue[nuclei$id==n])
      size = nuclei$size_red[nuclei$id==n]
    }
    else if (C=="G") {
      colors=rgb(nuclei$red[nuclei$id==n], 0, nuclei$blue[nuclei$id==n])
      size = nuclei$size_green[nuclei$id==n]
    }
    else {
      colors=rgb(nuclei$red[nuclei$id==n], nuclei$green[nuclei$id==n], 0)
      size = nuclei$size_blue[nuclei$id==n]
    }
    
    points(x[1], y[1], col=colors[1], pch=16)
    for (i in 2:length(x)) {lines(x[(i-1):i], y[(i-1):i], col=colors[i])}
    points(x[length(x)], y[length(y)], col=colors[length(x)], pch=1)
  }
}

#------------------------------------------------------------------------------------------------------------------------------------------------
# Horizontal bar graph of Ages of multiple nuclei identified by IDS as a function of time (in color)
Plot_Ages <- function(nuclei, IDs) {
  par(mfrow=c(1,1))
  
  x = nuclei$x_pos[nuclei$id==IDs[1]]
  y = nuclei$y_pos[nuclei$id==IDs[1]]
  n_frame = nuclei$n_frame[nuclei$id==IDs[1]]
  
  if(C == "R") {
    colors=rgb(0, nuclei$green[nuclei$id==IDs[1]], nuclei$blue[nuclei$id==IDs[1]])
    size = nuclei$size_red[nuclei$id==IDs[1]]
  }
  else if (C=="G") {
    colors=rgb(nuclei$red[nuclei$id==IDs[1]], 0, nuclei$blue[nuclei$id==IDs[1]])
    size = nuclei$size_green[nuclei$id==IDs[1]]
  }
  else {
    colors=rgb(nuclei$red[nuclei$id==IDs[1]], nuclei$green[nuclei$id==IDs[1]], 0)
    size = nuclei$size_blue[nuclei$id==IDs[1]]
  }
  
  plot(n_frame[1], IDs[1], col=colors[1], type="l", xlab="Age (frames)", ylab="ID",
       xlim=c(min(master$n_frame),max(master$n_frame)), ylim=c(0, max(IDs)), main="Age and Color of Nuclei")
  
  for(n in IDs) {
    x = nuclei$x_pos[nuclei$id==n]
    y = nuclei$y_pos[nuclei$id==n]
    n_frame = nuclei$n_frame[nuclei$id==n]
    
    if(C == "R") {
      colors=rgb(0, nuclei$green[nuclei$id==n], nuclei$blue[nuclei$id==n])
      size = nuclei$size_red[nuclei$id==n]
    }
    else if (C=="G") {
      colors=rgb(nuclei$red[nuclei$id==n], 0, nuclei$blue[nuclei$id==n])
      size = nuclei$size_green[nuclei$id==n]
    }
    else {
      colors=rgb(nuclei$red[nuclei$id==n], nuclei$green[nuclei$id==n], 0)
      size = nuclei$size_blue[nuclei$id==n]
    }
    
    for (i in 2:length(n_frame)) {lines(n_frame[(i-1):i], c(n, n), col=colors[i], lwd = 5)}
  }
}


#------------------------------------------------------------------------------------------------------------------------------------------------
# Plot Histogram of all the Nuclei Sizes
Plot_Histogram <- function(nuclei)
{
  if(C == "R") {
    size = nuclei$size_red[nuclei$id==n]
  }
  else if (C=="G") {
    size = nuclei$size_green[nuclei$id==n]
  }
  else {
    size = nuclei$size_blue[nuclei$id==n]
  }

    hist(size, nclass = 20, main="Nuclei Size Histogram", 
       xlab="Nuclei Size", ylab="Count", xaxt="n")
  axis(side=1,at=seq(0,5000,100))
}

#------------------------------------------------------------------------------------------------------------------------------------------------
# 3D Scatter plot of a single cell with ID in color
Plot_Scatterplot <- function(nuclei, ID)
{
  library(scatterplot3d)
  scatterplot3d(nuclei$x_pos[nuclei$id==ID],
                nuclei$y_pos[nuclei$id==ID],
                nuclei$time[nuclei$id==ID],
                color=rgb(nuclei$red[nuclei$id==ID], nuclei$green[nuclei$id==ID], nuclei$blue[nuclei$id==ID]))
}


#------------------------------------------------------------------------------------------------------------------------------------------------
# Plot Nuclei positions and colors per frame and save as images to externally create movie
Create_Position_Frames <- function(nuclei)
{
  for (f in nuclei$frame) {
    
    plot(nuclei$x_pos[nuclei$frame==f], nuclei$y_pos[nuclei$frame==f],
         col='black',
         main=paste("nuclei positions at frame ", f),
         xlab="x-position", ylab="y-position",
         xlim=c(min(x_crop), max(x_crop)), ylim=c(max(y_crop),min(y_crop)),
         pch=21,  bg=rgb(nuclei$red[nuclei$frame==f], nuclei$green[nuclei$frame==f], 0))
    dev.print(png, filename = paste(f, ".png", sep="") , width = (max(x_crop)-min(x_crop))/2, height = (max(y_crop)-min(y_crop))/2)
  }
}


#------------------------------------------------------------------------------------------------------------------------------------------------
# DRAFT - NEEDS UPDATING
#Plot Nuclei positions and search area on top of blue mask per frame and save as images to externally create movie
Overlay_Search_Area <- function(nuclei)
{
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  par(mfrow=c(1,1))
  
  for (i in Bimages) {
    print(i)
    # open image and get frame number
    img = readImage(paste(dir_Bprocessed, i, sep=""))
    f = substrRight(i, 10)
    f = as.numeric(substr(f, 1, 6))
    
    # draw circle
    for (n in which(nuclei$frame==f)) {
      img = drawCircle(img, nuclei$x_pos[n], nuclei$y_pos[n], sqrt(0.9*nuclei$size[n]/pi), col='white', fill=FALSE)
      img = drawCircle(img, nuclei$x_pos[n], nuclei$y_pos[n], sqrt(9*nuclei$size[n]/pi), col='gray', fill=FALSE)
    }
    
    # add labels
    display(img, method="raster")
    for (n in which(nuclei$frame==f)) {
      text(x = nuclei$x_pos[n], y = nuclei$y_pos[n], label = nuclei$id[n], adj = c(0.5, 0.5), col = "white", cex = 1.5)
    }
    
    # save image
    f_pad = sprintf("%03d", f)
    dev.print(jpeg, filename = paste("./data/annotated/annotated", f_pad, ".jpg", sep="") , width = dim(img)[1], height = dim(img)[2])
  }
}
#------------------------------------------------------------------------------------------------------------------------------------------------
# Use the overlay image, and successively add label and trace over time.
Overlay_Labels <- function(nuclei)
{
  # function to get frame number substring
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  par(mfrow=c(1,1))

    # iterate through all processed images
  Wimages = list.files(dir_Wprocessed)
  for (i in Wimages) {
    
    # print frame to console
    print(i)
    
    # open image and get frame number
    img = readImage(paste(dir_Wprocessed, i, sep=""))
    f = substrRight(i, 7)
    f = as.numeric(substr(f, 1, 3))
    print(f)
    
    # add labels
    display(img, method="raster")
    for (n in which(nuclei$n_frame==f-1)) {
      text(x = nuclei$x_pos[n], y = nuclei$y_pos[n], label = nuclei$id[n], adj = c(0.5, 0.5), col = "white", cex = 1.5)
    }
      
    # draw lines
    IDs = unique(nuclei$id)
    IDs = sort(IDs)
    
    for(n in IDs) {
      x = nuclei$x_pos[nuclei$id==n & nuclei$n_frame < f]
      y = nuclei$y_pos[nuclei$id==n & nuclei$n_frame < f]

      for (i in 2:length(x)) {lines(x[(i-1):i], y[(i-1):i], col="white")}
    }
    
    # save image
    f_pad = sprintf("%03d", f)
    dev.print(jpeg, filename = paste("./DATA/annotated/annotated", f_pad, ".jpg", sep="") , width = dim(img)[1], height = dim(img)[2])
  }
  
  
}
