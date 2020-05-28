#<><><><><><><><><><><><><><><><><><><><><><>
#------------------------------------
# FILTER AND ANALYZE NUCLEI DATAFRAME
#
# Nuclei co-located across different colors are identified
# and tagged with matching ID numbers
#------------------------------------
#<><><><><><><><><><><><><><><><><><><><><><>

# declare EBImage library for processing images
library(EBImage)

# run B-DataSettings.R script to pull-in variable settings
source("./B-DataSettings.R") 

# load data file back into nuclei dataframe
load(paste(dir_data, "nuclei_data.data", sep=""), envir = parent.frame(), verbose = TRUE)

# convert columns to numeric where applicable
nuclei$size= as.numeric(nuclei$size)
nuclei$intensity= as.numeric(nuclei$intensity)
nuclei$x_pos= as.numeric(nuclei$x_pos)
nuclei$y_pos= as.numeric(nuclei$y_pos)

# filter data by size of nuclei
nuclei = nuclei[which(nuclei$size>Nuclei_Lower_Size),] # remove any nuclei tagged smaller than Nuclei_Lower_Size
nuclei = nuclei[which(nuclei$size<Nuclei_Upper_Size),] # remove any nuclei tagged larger than Nuclei_Upper_Size

#----------------------------------------------
# Split into Data Frames Differentiated by Color / Channel
#----------------------------------------------

red    = nuclei[which(nuclei$color=="R"), ]
green  = nuclei[which(nuclei$color=="G"), ]
blue   = nuclei[which(nuclei$color=="B"), ]
master = nuclei[which(nuclei$color==C), ]

#----------------------------------------------
# MASTER: create a master dataframe to include nuclei data
# of base color, other colors present, and assign unique IDs
#----------------------------------------------

frames = unique(master$frame)     # names of all the unique frames #ANNA 5/1/2020#
frames = sort(frames)             # sort frames by alpha-numeric order

master$id = 0         # numerical nuclei (ID)entifier to remain consistent between frames
master$parent = 0     # numerical identifier of parent nucleus
c_id = 1              # current global numerical nuclei identifier

for(f in 1:(length(frames)-1))      # f is the index for the list of (f)rames
{                                   # skipping last frame becuase we reference f+1 in loop
  print(paste("Base Frame", frames[f])) # Update progress to console
  
  # Iterate through each nuclei (i)nstance within specified frame
  for(i in which(master$frame==frames[f]))
  {
    if(master$id[i] == 0) {
      master$id[i] = c_id   # give the nucleus an ID
      c_id = c_id+1         # increment current (next) nucleus ID
    }
    
    x_i = master$x_pos[i]         # x position of current (i)nstance nuclei in frame f
    y_i = master$y_pos[i]         # y position of current (i)nstance nuclei in frame f
    min_r2 = 10000000000          # minimum radial distance squared intentionally initialized very high
    min_n = 0                     # search for the index of nearest neighbor
    search_r2 = master$size[i]/pi # effective radius
    
    # AND iterate through nuclei in (n)ext frame
    for(n in which(master$frame==frames[f+1])) {
      
      x_n = master$x_pos[n]  # x position of nuclei in (n)ext frame f+1
      y_n = master$y_pos[n]  # y position of nuclei in (n)ext frame f+1
      
      r2 = (x_i-x_n)^2 + (y_i-y_n)^2 # radial distance between the two nuclei squared
      min_r2 = min(min_r2, r2)
      if(r2 == min_r2) {min_n = n}
      
      # if within small search radius,
      # and size has not changed significantly
      # assign matching ID
      if( r2 < 0.8*search_r2 &
          abs(master$size[n]-master$size[i]) < 0.2*master$size[n] ) {
        master$id[n] = master$id[i]
        master$parent[n] = master$id[i]
        min_n = 0
      }
    }

    # if no match was found,
    # use nearest neighbor in expanded search radius 
    if(min_r2 < 9*search_r2 & min_n > 0 ) {
      if(abs(master$size[min_n]-master$size[i]) < .5*master$size[min_n]) {
        master$id[min_n] = master$id[i]
        master$parent[min_n] = master$id[i]
      }
    }
  }

  # Create new cell IDs when duplicate IDs arise in the same frame
  #----------------------------------------------
  
  # create a table of all nuclei in next frame
  nuclei_in_frame = table(master$id[master$frame==frames[f+1]])
  
  # iterate through (t)able to find duplicate nuclei ID
  for(t in 1:length(nuclei_in_frame)) {
    
    # if duplicate ID is found
    if (nuclei_in_frame[[t]]> 1) {
      
      # set parent ID to current ID  ... names(nuclei_in_frame[t]) ...
      master$parent[master$frame==frames[f+1] & master$id==as.numeric(names(nuclei_in_frame[t]))] = as.numeric(names(nuclei_in_frame[t]))
      
      # provide new IDs to the "child" nuclei
      for(n in master$id[master$frame==frames[f+1] & master$id==as.numeric(names(nuclei_in_frame[t]))]) {
        master$id[n]= c_id
        c_id = c_id+1 # increment current (next) nucleus ID
        print(paste("Dup. Loop - Frame", frames[f], "c_id", c_id))
      }
    }
  }
}

# Complete the assignments for new nuclei in the last frame
f=max(frames)            # f is the last index in the list of (f)rames
print(paste("Base Frame", frames[f])) # Update progress to console

# Iterate through each nuclei (i)nstance within last frame
for(i in which(master$frame==frames[f])) {
  if(master$id[i] == 0) {
    master$id[i] = c_id   # give the nucleus an ID
    c_id = c_id+1         # increment current (next) nucleus ID
  }
}

#----------------------------------------------
# RED: search nuclei in red frame for nearest neighbor
# in master frame (base color) and assign common ID
#----------------------------------------------
master$red = FALSE
master$size_red = 0
master$intensity_red = 0

if(R_toggle) {
  red$id = 0
  
  for(f in frames) # Iterate through all (f)rames
  {
    # Iterate through each (r)ed nuclei instance within specified frame
    for(r in which(red$frame==f))
    {
      x_r = red$x_pos[r]    # x position of (r)ed nuclei in frame f
      y_r = red$y_pos[r]    # y position of (r)ed nuclei in frame f
      min_r2 = 10000000000  # minimum radial distance squared intentionally initialized very high
      
      # AND iterate through nuclei within matching (m)aster frame
      for(m in which(master$frame==f))
      { 
        x_m = master$x_pos[m]  # x position of (m)aster nuclei in frame f
        y_m = master$y_pos[m]  # y position of (m)aster nuclei in frame f
        
        r2 = (x_r-x_m)^2 + (y_r-y_m)^2 # radial distance between the two nuclei squared
        
        if(r2 < min_r2) {               # If the nuclei pair is closer than any other pair previosly found,
          min_r2 = r2                   # re-assign minimum radial distance found.
          if (r2 < master$size[m]/pi) { # And, if the pair is within 1 radius distance of the master nuclei,
            master$red[m] = TRUE        # assign RED values MASTER dataframe,
            master$size_red[m] = red$size[r]
            master$intensity_red[m] = red$intensity[r]
            red$id[r] = master$id[m]    # and assign the RED id to match the master ID.
          }
        }
      }
      print(paste("red nuclei", r, "in frame", f))
    }
  }
}

#----------------------------------------------
# GREEN: search nuclei in green frame for nearest neighbor
# in master frame (base color) and assign common ID
#----------------------------------------------
master$green = FALSE
master$size_green = 0
master$intensity_green = 0

if(G_toggle) {
  green$id = 0
  
  for(f in frames) # Iterate through all (f)rames
  {
    # Iterate through each (g)reen nuclei instance within specified frame
    for(g in which(green$frame==f))
    {
      x_g = green$x_pos[g]  # x position of (g)reen nuclei in frame f
      y_g = green$y_pos[g]  # y position of (g)reen nuclei in frame f
      min_r2 = 10000000000  # minimum radial distance squared intentionally initialized very high
      
      # AND iterate through nuclei within matching (m)aster frame
      for(m in which(master$frame==f))
      { 
        x_m = master$x_pos[m]  # x position of (m)aster nuclei in frame f
        y_m = master$y_pos[m]  # y position of (m)aster nuclei in frame f
        
        r2 = (x_g-x_m)^2 + (y_g-y_m)^2 # radial distance between the two nuclei squared
        
        if(r2 < min_r2) {               # If the nuclei pair is closer than any other pair previosly found,
          min_r2 = r2                   # re-assign minimum radial distance found.
          if (r2 < master$size[m]/pi) { # And, if the pair is within 1 radius distance of the master nuclei,
            master$green[m] = TRUE      # assign GREEN values MASTER dataframe,
            master$size_green[m] = green$size[g]
            master$intensity_green[m] = green$intensity[g]
            green$id[g] = master$id[m]   # and assign the GREEN id to match the master ID.
          }
        }
      }
      print(paste("green nuclei", g, "in frame", f))
    }
  }
}

#----------------------------------------------
# BLUE: search nuclei in blue frame for nearest neighbor
# in master frame (base color) and assign common ID
#----------------------------------------------
master$blue = FALSE
master$size_blue = 0
master$intensity_blue = 0

if(B_toggle) {
  blue$id = 0
  
  for(f in frames) # Iterate through all (f)rames
  {
    # Iterate through each (b)lue nuclei instance within specified frame
    for(b in which(blue$frame==f))
    {
      x_b = blue$x_pos[b]   # x position of (b)lue nuclei in frame f
      y_b = blue$y_pos[b]   # y position of (b)lue nuclei in frame f
      min_r2 = 10000000000  # minimum radial distance squared intentionally initialized very high
      
      # AND iterate through nuclei within matching (m)aster frame
      for(m in which(master$frame==f))
      { 
        x_m = master$x_pos[m]  # x position of (m)aster nuclei in frame f
        y_m = master$y_pos[m]  # y position of (m)aster nuclei in frame f
        
        r2 = (x_b-x_m)^2 + (y_b-y_m)^2 # radial distance between the two nuclei squared
        
        if(r2 < min_r2) {               # If the nuclei pair is closer than any other pair previosly found,
          min_r2 = r2                   # re-assign minimum radial distance found.
          if (r2 < master$size[m]/pi) { # And, if the pair is within 1 radius distance of the master nuclei,
            master$blue[m] = TRUE       # assign BLUE values MASTER dataframe,
            master$size_blue[m] = blue$size[b]
            master$intensity_blue[m] = blue$intensity[b]
            blue$id[b] = master$id[m]   # and assign the BLUE id to match the master ID.
          }
        }
      }
      print(paste("blue nuclei", b, "in frame", f))
    }
  }
}

colnames(master)[colnames(master)=="color"] <- "base_color"
master$size = NULL
master$intensity = NULL

#----------------------------------------------
# Count Number of Cells in a Frame by Colors
#----------------------------------------------

cellcount = data.frame(
  R = integer(),
  G = integer(),
  B = integer(),
  RG = integer(),
  RB = integer(),
  GB = integer(),
  RGB = integer(),
  stringsAsFactors = FALSE)


for(f in 1:length(frames)) { 
  cellcount[f, "R"] = length(master[master$frame==frames[f] &
                                      master$red == TRUE &
                                      master$green == FALSE &
                                      master$blue == FALSE,1])
  cellcount[f, "G"] = length(master[master$frame==frames[f] &
                                      master$red == FALSE &
                                      master$green == TRUE &
                                      master$blue == FALSE,1])
  cellcount[f, "B"] = length(master[master$frame==frames[f] &
                                      master$red == FALSE &
                                      master$green == FALSE &
                                      master$blue == TRUE,1])
  cellcount[f, "RG"] = length(master[master$frame==frames[f] &
                                      master$red == TRUE &
                                      master$green == TRUE &
                                      master$blue == FALSE,1])
  cellcount[f, "RB"] = length(master[master$frame==frames[f] &
                                      master$red == TRUE &
                                      master$green == FALSE &
                                      master$blue == TRUE,1])
  cellcount[f, "GB"] = length(master[master$frame==frames[f] &
                                      master$red == FALSE &
                                      master$green == TRUE &
                                      master$blue == TRUE,1])
  cellcount[f, "RGB"] = length(master[master$frame==frames[f] &
                                      master$red == TRUE &
                                      master$green == TRUE &
                                      master$blue == TRUE,1])
}

#------------------------------------------------------------------
# SAVE results
#------------------------------------------------------------------
# SAVE results into native R data files
save(cellcount, file=paste(dir_data, "cell_count.data",    sep=""))
save(red,       file=paste(dir_data, "nuclei_red.data",    sep=""))
save(green,     file=paste(dir_data, "nuclei_green.data",  sep=""))
save(blue,      file=paste(dir_data, "nuclei_blue.data",   sep=""))
save(master,    file=paste(dir_data, "nuclei_master.data", sep=""))

# SAVE results as comma seperated value files
write.table(cellcount, file=paste(dir_data, "cell_counts.csv", sep=""),   row.names = FALSE, col.names = TRUE, sep = ", ")
write.table(red,       file=paste(dir_data, "nuclei_red.csv", sep=""),    row.names = FALSE, col.names = TRUE, sep = ", ")
write.table(green,     file=paste(dir_data, "nuclei_green.csv", sep=""),  row.names = FALSE, col.names = TRUE, sep = ", ")
write.table(blue,      file=paste(dir_data, "nuclei_blue.csv", sep=""),   row.names = FALSE, col.names = TRUE, sep = ", ")
write.table(master,    file=paste(dir_data, "nuclei_master.csv", sep=""), row.names = FALSE, col.names = TRUE, sep = ", ")

#----------------------------------------------
# Tag Overlays with Nuclei IDs and Colors
#----------------------------------------------

for(f in frames) # Iterate through all (f)rames
{
  img = readImage(paste(dir_Wprocessed, f, sep="")) # open overlay image
  display(img, method="raster")                     # display image

  for (n in which(master$frame==f)) {               # add labels to each (n)uclei
    label_color = rgb(master$red[n], master$green[n], master$blue[n])
    text(x = master$x_pos[n], y = master$y_pos[n],
         label = master$id[n], adj = c(-0.2,.5),
         col = label_color, cex = 1)
  }
  # save image with labels
  tiff(filename = paste(dir_Wprocessed, "/labeled-", f, sep=""))
  dev.off()
  print(paste(dir_data, "labeled-", f, sep =""))
}

