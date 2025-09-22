# Read Data ####
library(astrochron)
library(extrafont)
library(WaverideR)
library(DescTools)
library(magick) 
library(png)
library(readxl)
library(parallel)
library(extrafont)
library(showtext)
font_add("Calibri", "C:/Windows/Fonts/calibri.ttf")  # adjust if file is elsewhere

readPNG_from_github <- function(url) {
  tmp <- tempfile(fileext = ".png")
  download.file(url, tmp, mode = "wb", quiet = TRUE)
  readPNG(tmp)
}


#images with climate and lithology###
img  <- readPNG_from_github("https://raw.githubusercontent.com/stratigraphy/cyclo_strat_gifs/main/limestone_1.png")
img3 <- readPNG_from_github("https://raw.githubusercontent.com/stratigraphy/cyclo_strat_gifs/main/marl_1.png")
img5 <- readPNG_from_github("https://raw.githubusercontent.com/stratigraphy/cyclo_strat_gifs/main/shale_1.png")

img_naut_1 <- readPNG_from_github("https://raw.githubusercontent.com/stratigraphy/cyclo_strat_gifs/main/Silurian_naut_1.png")
img_naut_2 <- readPNG_from_github("https://raw.githubusercontent.com/stratigraphy/cyclo_strat_gifs/main/Silurian_naut_2.png")
img_naut_3 <- readPNG_from_github("https://raw.githubusercontent.com/stratigraphy/cyclo_strat_gifs/main/Silurian_naut_3.png")


graphics.off()
ALT_GR <- getLaskar(sol="insolation")
ALT_GR <- ALT_GR[, c(1, 2)]
depth <- ALT_GR[, 1]
GR <- ALT_GR[, 2]
gr <- GR

#generate Litholog based on thresholds ####

lower_threshold <- mean(GR) - 0.5 * sd(GR)
middle_threshold <- mean(GR) + sd(GR)

get_lithology <- function(gr) {
  if (is.na(gr)) {
    return(NA)
  } else if (gr < lower_threshold) {
    return("Limestone")
  } else if (gr < middle_threshold) {
    return("Marl")
  } else {
    return("Shale")
  }
}

lithology <- sapply(GR, get_lithology)

lith_factor <- factor(lithology, levels = c("Limestone", "Marl", "Shale"))
bed_changes <- which(diff(as.numeric(lith_factor)) != 0) + 1
bed_edges <- c(1, bed_changes, length(depth) + 1)
bed_type <- lithology[bed_edges[-length(bed_edges)]]

# Lithology-specific widths and colors 
widths <- list(Limestone = 2.5, Marl = 2, Shale = 1.75)
colors <- list(Limestone = "#F4F4E4", Marl = "#E0E0DE", Shale = "#7C7C7C")

# set plot limits 
fixed_xlim <- c(0, 300)  


###plot litholog and insolation curve static ####
graphics.off()
layout(matrix(c(1, 2), nrow = 2), heights = c(2, 1.5))  # stack vertically

par(mar = c(4, 4, 2, 8), family = "Calibri")  # extra margin on right for legend
plot(NA, ylim = c(0, max(unlist(widths))), xlim = fixed_xlim,
     ylab = "Lithology", xlab = "Time kyr",
     main = "Litholog", yaxt = "n", bty = "n")
axis(1, las = 1)

for (i in 1:(length(bed_edges) - 1)) {
  lith <- bed_type[i]
  width <- widths[[lith]]
  
  top_depth <- depth[bed_edges[i]]
  base_depth <- depth[bed_edges[i + 1] - 1]
  
  if (base_depth >= 0 && top_depth <= 300) {
    poly_left <- max(top_depth, 0)
    poly_right <- min(base_depth, 300)
    
    x <- c(poly_left, poly_right, poly_right, poly_left)
    y <- c(0, 0, width, width)
    
    polygon(x, y, col = colors[[lith]], border = "black", lwd = 0.5)
  }
}

# Legend moved outside to the right
legend("right", inset = c(-0.25, 0), legend = names(colors),
       fill = unlist(colors), bty = "n", border = "black", xpd = TRUE)

#  GR Log (rotated) 
valid_gr <- !is.na(GR) & !is.na(depth)
plot(depth[valid_gr], GR[valid_gr], type = "l", col = "darkgreen", lwd = 1.5,
     xlab = "Time kyr", ylab = "insolation (W/m^2)",
     ylim = range(GR, na.rm = TRUE), xlim = fixed_xlim,
     main = "insolation")


polygon(x=c(fixed_xlim[1]-50,fixed_xlim[2]+250,fixed_xlim[2]+250,fixed_xlim[1]-50),
        y=c(middle_threshold,middle_threshold,650,650),
        col = as.character(colors[1]))

polygon(x=c(fixed_xlim[1]-50,fixed_xlim[2]+250,fixed_xlim[2]+250,fixed_xlim[1]-50),
        y=c(middle_threshold,middle_threshold,lower_threshold,lower_threshold),
        col = as.character(colors[2]))

polygon(x=c(fixed_xlim[1]-50,fixed_xlim[2]+250,fixed_xlim[2]+250,fixed_xlim[1]-50),
        y=c(lower_threshold,lower_threshold,400,400),
        col = as.character(colors[3]))


lines(depth[valid_gr], GR[valid_gr], type = "l", col = "black", lwd = 3,
     ylim = range(GR, na.rm = TRUE), xlim = fixed_xlim)

abline(h = c(lower_threshold, middle_threshold), col = "red", lty = 2)



#make it growing boxes
lith_widths <- cbind(depth,as.numeric( unlist(widths[as.character(lith_factor)]) ))
lith_widths_405 <- taner(hilbert(lith_widths,genplot=FALSE),flow=1/505,fhigh=1/305,roll=10^20,xmax=1/50)
lith_widths_100 <- taner(hilbert(lith_widths,genplot=FALSE),flow=1/150,fhigh=1/85,roll=10^20,xmax=1/50)

mins_405 <- min_detect(lith_widths_405)
mins_100 <- min_detect(lith_widths_100)

#top <- c(0,0.25,-1)
#names(top) <- colnames(mins_100)
#mins_100 <- rbind(top,mins_100)
#mins_405 <- rbind(top,mins_405)

#  Plot minima as horizontal boxes 

# Parameters for y placement
ymax <- max(unlist(widths))

# # Draw boxes for mins_405
# for (i in 1:(nrow(mins_405)-1)) {
#   x  <- mins_405[i, 1]
#   x2 <- mins_405[i+1, 1]
#   
#   rect(x, ymax*0.7, x2, ymax*0.8,
#        col = "pink", border = "black", lwd = 2)
#   
#   # Add centered text
#   text((x + x2)/2, (ymax*0.7 + ymax*0.8)/2,
#        labels = paste0("405-kyr"), cex = 0.8)
# }
# 
# # Draw boxes for mins_100
# for (i in 1:(nrow(mins_100)-1)) {
#   x  <- mins_100[i, 1]
#   x2 <- mins_100[i+1, 1]
#   
#   rect(x, ymax*0.9, x2, ymax,
#        col = "lightblue", border = "black", lwd = 2)
#   
#   # Add centered text
#   text((x + x2)/2, (ymax*0.9 + ymax)/2,
#        labels = paste0("~100-kyr"), cex = 0.8)
# }
# 
# 




##animate cyclostrat time ####
# animates the translation from the insolation curve to litholog
# 3 windows: 
#window 1: insolation curve
#window 2: litholog and the 100-kyr ecc curve
#window 3: lithoplog and the 405-kyr ecc curve
#setup 
n_beds <- length(bed_edges) - 1
fixed_xlim <- c(0, 2000)  
time_steps <- seq(0, fixed_xlim[2], by = 5)  # every 10 kyr
n_steps <- length(time_steps)

# set the working directory for the generated images and the gif
setwd("D:/Phd/R/R_results/gif_5")

# Loop to draw frame by frame 
n_beds <-sum(depth[bed_edges] < fixed_xlim[2], na.rm = TRUE)
getwd()
for (k in 1:n_steps) {
  #k  <- 20
  current_time <- time_steps[k]
  png(sprintf("frame_%03d.png", k), width = 1600, height = 800, res = 150)

  #graphics.off()
  layout(matrix(c(1, 2,3), nrow = 3,ncol=1), heights = c(1,1,1))
  par(mar = c(4, 4, 2, 8), family = "Calibri")

  
  valid_gr <- !is.na(GR) & !is.na(depth)
  plot(depth[valid_gr], GR[valid_gr], type = "l", col = "darkgreen", lwd = 1.5,
       xlab = "Time (kyr)", ylab = "Insolation (W/m^2)",
       ylim = range(GR, na.rm = TRUE), xlim = fixed_xlim,
       main = "insolation",
       cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.2)
  abline(h = c(lower_threshold, middle_threshold), col = "red", lty = 2)
  
  
  polygon(x=c(fixed_xlim[1]-100,fixed_xlim[2]+250,fixed_xlim[2]+250,fixed_xlim[1]-100),
          y=c(middle_threshold,middle_threshold,650,650),
          col = as.character(colors[1]))
  
  polygon(x=c(fixed_xlim[1]-100,fixed_xlim[2]+250,fixed_xlim[2]+250,fixed_xlim[1]-100),
          y=c(middle_threshold,middle_threshold,lower_threshold,lower_threshold),
          col = as.character(colors[2]))
  
  polygon(x=c(fixed_xlim[1]-100,fixed_xlim[2]+250,fixed_xlim[2]+250,fixed_xlim[1]-100),
          y=c(lower_threshold,lower_threshold,400,400),
          col = as.character(colors[3]))
  
  lines(depth[valid_gr], GR[valid_gr], col = "red", lwd = 2,
        ylim = range(GR, na.rm = TRUE), xlim = fixed_xlim)
  
  abline(h = c(lower_threshold, middle_threshold), col = "red", lty = 2)
  
  # Add a point for the current bed midpoint

  if (current_time <= max(depth)) {
    idx <- which.min(abs(depth - current_time))  
    points(depth[idx], GR[idx],
           col = "blue", pch = 19, cex = 2)
  }

  legend("right", inset = c(-0.12, 0), legend = names(colors),
         fill = unlist(colors), bty = "n", border = "black",
         xpd = TRUE)
  

  
  plot(NA, ylim = c(0, 1.5*max(unlist(widths))), xlim = fixed_xlim,
       ylab = "Lithology", 
       main = "Litholog ~100-kyr eccentricity cycle", 
       yaxt = "n", bty = "n",
       cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.2,xaxt="n",xlab="")
  
  for (j in 1:(length(bed_edges)-1)) {
    lith <- bed_type[j]
    width <- widths[[lith]]
    
    top_depth <- depth[bed_edges[j]]
    base_depth <- depth[bed_edges[j + 1] - 1]
    
    if (top_depth < current_time) {
      poly_left <- max(top_depth, 0)
      poly_right <- min(base_depth, current_time)
      
      if (poly_left < poly_right) {
        x <- c(poly_left, poly_right, poly_right, poly_left)
        y <- c(0, 0, width, width)
        polygon(x, y, col = colors[[lith]], border = "black", lwd = 0.5)
      }
    }
  }
  

  legend("right", inset = c(-0.12, 0), legend = names(colors),
         fill = unlist(colors), bty = "n", border = "black",
         xpd = TRUE)
  
  
  lines(lith_widths_100[,1],4*lith_widths_100[,2]+0.5)
  
  
  if (current_time <= max(lith_widths_100[,1])) {
    idx <- which.min(abs(lith_widths_100[,1] - current_time))  
    points(lith_widths_100[idx,1], 4*lith_widths_100[idx,2]+0.5,
           col = "blue", pch = 19, cex = 2)
  }
  
  mid_depth <- lith_widths_100[idx,1]
  mins_100_sel <- mins_100
  mins_100_sel[mins_100_sel[,1]<mid_depth,]
  mins_100_sel <- mins_100_sel[mins_100_sel[,1]<mid_depth,]
  
  if(nrow(mins_100_sel)>1){
    
  for (ijk in 1:(nrow(mins_100_sel)-1)) {
    x  <- mins_100_sel[ijk, 1]
    x2 <- mins_100_sel[ijk+1, 1]
    
    
    rect(x, ymax*1, x2, ymax*1.5,
         col = "lightblue", border = "black", lwd = 2)
    
    

    
    # Add centered text
    text((x + x2)/2, (ymax*1.5 + ymax)/2,
         labels = paste0("~100-kyr"), cex = 0.9)
  }
    
    rect(x2, ymax*1, mid_depth, ymax*1.5,
         col = "lightblue", border = "black", lwd = 2)
    
    
  }
  
  if(nrow(mins_100_sel)==1){
    x2 <-  mins_100_sel[,1]
    
    rect(x2, ymax*1, mid_depth, ymax*1.5,
         col = "lightblue", border = "black", lwd = 2)
  }
  
  
  
  plot(NA, ylim = c(0, 1.5*max(unlist(widths))), xlim = fixed_xlim,
       ylab = "Lithology",
       main = "Litholog 405-kyr eccentricity cycle", yaxt = "n", bty = "n",
       cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.2,xaxt="n",xlab="")


  for (j in 1:(length(bed_edges)-1)) {
    lith <- bed_type[j]
    width <- widths[[lith]]
    
    top_depth <- depth[bed_edges[j]]
    base_depth <- depth[bed_edges[j + 1] - 1]
    
    if (top_depth < current_time) {
      poly_left <- max(top_depth, 0)
      poly_right <- min(base_depth, current_time)
      
      if (poly_left < poly_right) {
        x <- c(poly_left, poly_right, poly_right, poly_left)
        y <- c(0, 0, width, width)
        polygon(x, y, col = colors[[lith]], border = "black", lwd = 0.5)
      }
    }
  }
  
  
  legend("right", inset = c(-0.12, 0), legend = names(colors),
         fill = unlist(colors), bty = "n", border = "black",
         xpd = TRUE)
  
  
  lines(lith_widths_405[,1],4*lith_widths_405[,2]+0.5)
  
  if (current_time <= max(lith_widths_405[,1])) {
    idx <- which.min(abs(lith_widths_405[,1] - current_time))  
    points(lith_widths_405[idx,1], 4*lith_widths_405[idx,2]+0.5,
           col = "blue", pch = 19, cex = 2)
  }
  
  mid_depth <- lith_widths_405[idx,1]
  
  mins_405_sel <- mins_405
  mins_405_sel[mins_405_sel[,1]<mid_depth,]
  mins_405_sel <- mins_405_sel[mins_405_sel[,1]<mid_depth,]
  
  
  
  if(nrow(mins_405_sel)>1){
  
  for (ijk in 1:(nrow(mins_405_sel)-1)) {
    x  <- mins_405_sel[ijk, 1]
    x2 <- mins_405_sel[ijk+1, 1]
    
    
    rect(x, ymax*1, x2, ymax*1.5,
         col = "lightblue", border = "black", lwd = 2)
    
    
        # Add centered text
    text((x + x2)/2, (ymax*1.5 + ymax)/2,
         labels = paste0("405-kyr"), cex = 2)
  }
    rect(x2, ymax*1, mid_depth, ymax*1.5,
         col = "lightblue", border = "black", lwd = 2)
    
  }
  
  if(nrow(mins_405_sel)==1){
    x2 <-  mins_405_sel[,1]
    
  rect(x2, ymax*1, mid_depth, ymax*1.5,
       col = "lightblue", border = "black", lwd = 2)
  }
  
  
  
  
  dev.off()
  
  }
#getwd()
gc()
graphics.off()
frames <- list.files(pattern = "frame_\\d+\\.png", full.names = TRUE)
img_list <- image_read(frames)
animation <- image_animate(img_list, fps = 20, loop = 10)  # fps=4 frames/sec
image_write(animation, "animation_5.gif")


##animate cyclostrat no time axis #####
# 2 windows: 
#window 1: litholog and the 100-kyr ecc curve
#window 2: lithoplog and the 405-kyr ecc curve
# set the working directory for the generated images and the gif
setwd("D:/Phd/R/R_results/gif_6")

n_beds <-sum(depth[bed_edges] < fixed_xlim[2], na.rm = TRUE)
getwd()
for (k in 1:n_steps) {
  #k  <- 150
  current_time <- time_steps[k]
  png(sprintf("frame_%03d.png", k), width = 1600, height = 800, res = 150)
  
  #graphics.off()
  layout(matrix(c(1, 2), nrow = 2,ncol=1), heights = c(1,1,1))
  par(mar = c(4, 4, 2, 8), family = "Calibri")
  
  plot(NA, ylim = c(0, 1.5*max(unlist(widths))), xlim = fixed_xlim,
       ylab = "Lithology", xlab = "",
       main = "Litholog ~100-kyr eccentricity cycle", 
       yaxt = "n", bty = "n",xaxt="n",
       cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.2)
  
  for (j in 1:(length(bed_edges)-1)) {
    lith <- bed_type[j]
    width <- widths[[lith]]
    
    top_depth <- depth[bed_edges[j]]
    base_depth <- depth[bed_edges[j + 1] - 1]
    
    if (top_depth < current_time) {
      poly_left <- max(top_depth, 0)
      poly_right <- min(base_depth, current_time)
      
      if (poly_left < poly_right) {
        x <- c(poly_left, poly_right, poly_right, poly_left)
        y <- c(0, 0, width, width)
        polygon(x, y, col = colors[[lith]], border = "black", lwd = 0.5)
      }
    }
  }
  
  
  legend("right", inset = c(-0.15, 0), legend = names(colors),
         fill = unlist(colors), bty = "n", border = "black",
         xpd = TRUE)
  
  
  lines(lith_widths_100[,1],4*lith_widths_100[,2]+0.5)
  
  
  if (current_time <= max(lith_widths_100[,1])) {
    idx <- which.min(abs(lith_widths_100[,1] - current_time))  
    points(lith_widths_100[idx,1], 4*lith_widths_100[idx,2]+0.5,
           col = "blue", pch = 19, cex = 2)
  }
  
  mid_depth <- lith_widths_100[idx,1]
  mins_100_sel <- mins_100
  mins_100_sel[mins_100_sel[,1]<mid_depth,]
  mins_100_sel <- mins_100_sel[mins_100_sel[,1]<mid_depth,]
  
  if(nrow(mins_100_sel)>1){
    
    for (ijk in 1:(nrow(mins_100_sel)-1)) {
      x  <- mins_100_sel[ijk, 1]
      x2 <- mins_100_sel[ijk+1, 1]
      
      
      rect(x, ymax*1, x2, ymax*1.5,
           col = "lightblue", border = "black", lwd = 2)
      
      
      
      
      # Add centered text
      text((x + x2)/2, (ymax*1.5 + ymax)/2,
           labels = paste0("~100-kyr"), cex = 0.6)
    }
    
    rect(x2, ymax*1, mid_depth, ymax*1.5,
         col = "lightblue", border = "black", lwd = 2)
    
    
  }
  
  if(nrow(mins_100_sel)==1){
    x2 <-  mins_100_sel[,1]
    
    rect(x2, ymax*1, mid_depth, ymax*1.5,
         col = "lightblue", border = "black", lwd = 2)
  }
  
  
  
  plot(NA, ylim = c(0, 1.5*max(unlist(widths))), xlim = fixed_xlim,
       ylab = "Lithology", xlab = "",xaxt="n",
       main = "Litholog 405-kyr eccentricity cycle", yaxt = "n", bty = "n",
       cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.2)
  
  
  for (j in 1:(length(bed_edges)-1)) {
    lith <- bed_type[j]
    width <- widths[[lith]]
    
    top_depth <- depth[bed_edges[j]]
    base_depth <- depth[bed_edges[j + 1] - 1]
    
    if (top_depth < current_time) {
      poly_left <- max(top_depth, 0)
      poly_right <- min(base_depth, current_time)
      
      if (poly_left < poly_right) {
        x <- c(poly_left, poly_right, poly_right, poly_left)
        y <- c(0, 0, width, width)
        polygon(x, y, col = colors[[lith]], border = "black", lwd = 0.5)
      }
    }
  }
  
  
  legend("right", inset = c(-0.15, 0), legend = names(colors),
         fill = unlist(colors), bty = "n", border = "black",
         xpd = TRUE)
  
  
  lines(lith_widths_405[,1],4*lith_widths_405[,2]+0.5)
  
  if (current_time <= max(lith_widths_405[,1])) {
    idx <- which.min(abs(lith_widths_405[,1] - current_time))  
    points(lith_widths_405[idx,1], 4*lith_widths_405[idx,2]+0.5,
           col = "blue", pch = 19, cex = 2)
  }
  
  mid_depth <- lith_widths_405[idx,1]
  
  mins_405_sel <- mins_405
  mins_405_sel[mins_405_sel[,1]<mid_depth,]
  mins_405_sel <- mins_405_sel[mins_405_sel[,1]<mid_depth,]
  
  
  
  if(nrow(mins_405_sel)>1){
    
    for (ijk in 1:(nrow(mins_405_sel)-1)) {
      x  <- mins_405_sel[ijk, 1]
      x2 <- mins_405_sel[ijk+1, 1]
      
      
      rect(x, ymax*1, x2, ymax*1.5,
           col = "lightblue", border = "black", lwd = 2)
      
      
      # Add centered text
      text((x + x2)/2, (ymax*1.5 + ymax)/2,
           labels = paste0("405-kyr"), cex = 2)
    }
    rect(x2, ymax*1, mid_depth, ymax*1.5,
         col = "lightblue", border = "black", lwd = 2)
    
  }
  
  if(nrow(mins_405_sel)==1){
    x2 <-  mins_405_sel[,1]
    
    rect(x2, ymax*1, mid_depth, ymax*1.5,
         col = "lightblue", border = "black", lwd = 2)
  }
  
  
  
  
  dev.off()
  
}

gc()
graphics.off()
frames <- list.files(pattern = "frame_\\d+\\.png", full.names = TRUE)
img_list <- image_read(frames)
animation <- image_animate(img_list, fps = 25, loop = 10)  # fps=4 frames/sec
image_write(animation, "animation_6.gif")

### Animate litholog on top insolation on the bottom#####
graphics.off()
ALT_GR <- getLaskar(sol="insolation")
ALT_GR <- ALT_GR[, c(1, 2)]
depth <- ALT_GR[, 1]
GR <- ALT_GR[, 2]
gr <- GR

lower_threshold <- mean(GR) - 0.5 * sd(GR)
middle_threshold <- mean(GR) + sd(GR)

get_lithology <- function(gr) {
  if (is.na(gr)) {
    return(NA)
  } else if (gr < lower_threshold) {
    return("Limestone")
  } else if (gr < middle_threshold) {
    return("Marl")
  } else {
    return("Shale")
  }
}

lithology <- sapply(GR, get_lithology)

# Identify bed boundaries
lith_factor <- factor(lithology, levels = c("Limestone", "Marl", "Shale"))
bed_changes <- which(diff(as.numeric(lith_factor)) != 0) + 1
bed_edges <- c(1, bed_changes, length(depth) + 1)
bed_type <- lithology[bed_edges[-length(bed_edges)]]

#  Lithology-specific widths and colors 
widths <- list(Limestone = 2.5, Marl = 2, Shale = 1.75)
colors <- list(Limestone = "#F4F4E4", Marl = "#E0E0DE", Shale = "#7C7C7C")

# Plotting setup 
fixed_xlim <- c(0, 300)  # time/ depth window
n_beds <- length(bed_edges) - 1

# set the working directory for the generated images and the gif
setwd("D:/Phd/R/R_results/gif")
getwd()


# Loop to draw frame by frame 
n_beds <-sum(depth[bed_edges] < 300, na.rm = TRUE)

for (i in 1:n_beds) {
  png(sprintf("frame_%03d.png", i), width = 1200, height = 800, res = 150)
  #graphics.off()
  layout(matrix(c(1, 2), nrow = 2), heights = c(2, 1.5))
  
  #  Top: Lithology column 
  par(mar = c(4, 4, 2, 8), family = "Calibri")
  plot(NA, ylim = c(0, max(unlist(widths))), xlim = fixed_xlim,
       ylab = "Lithology", xlab = "Time (kyr)",
       main = "Litholog (building step by step)", yaxt = "n", bty = "n",
       cex.main = 1.6,   # title size
       cex.lab = 1.4,    # axis labels size
       cex.axis = 1.2)   # axis tick labels size
  axis(1, las = 1, cex.axis = 1.2)
  
  for (j in 1:i) {
    lith <- bed_type[j]
    width <- widths[[lith]]
    
    top_depth <- depth[bed_edges[j]]
    base_depth <- depth[bed_edges[j + 1] - 1]
    
    if (base_depth >= 0 && top_depth <= 300) {
      poly_left <- max(top_depth, 0)
      poly_right <- min(base_depth, 300)
      
      x <- c(poly_left, poly_right, poly_right, poly_left)
      y <- c(0, 0, width, width)
      
      polygon(x, y, col = colors[[lith]], border = "black", lwd = 0.5)
    }
  }
  
  legend("right", inset = c(-0.25, 0), legend = names(colors),
         fill = unlist(colors), bty = "n", border = "black", xpd = TRUE,
         cex = 1.2)  # legend text size
  
  #  Bottom: GR log 
  par(mar = c(4, 4, 2, 8), family = "Calibri")
  valid_gr <- !is.na(GR) & !is.na(depth)
  
  plot(depth[valid_gr], GR[valid_gr], type = "l", col = "darkgreen", lwd = 1.5,
       xlab = "Time (kyr)", ylab = "insolation (W/m^2)",
       ylim = range(GR, na.rm = TRUE), xlim = fixed_xlim,
       main = "insolation",
       cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.2)
  
  abline(h = c(lower_threshold, middle_threshold), col = "red", lty = 2)
  
  polygon(x=c(fixed_xlim[1]-100,650,fixed_xlim[2]+100,fixed_xlim[1]-100), y=c(middle_threshold,middle_threshold,650,650),
          col = as.character(colors[1]))
  
  polygon(x=c(fixed_xlim[1]-100,650,fixed_xlim[2]+100,fixed_xlim[1]-100), y=c(middle_threshold,middle_threshold,lower_threshold,lower_threshold),
          col = as.character(colors[2]))
  
  polygon(x=c(fixed_xlim[1]-100,650,fixed_xlim[2]+100,fixed_xlim[1]-100), y=c(lower_threshold,lower_threshold,400,400),
          col = as.character(colors[3]))
  
  lines(depth[valid_gr], GR[valid_gr], col = "red", lwd = 2,
        ylim = range(GR, na.rm = TRUE), xlim = fixed_xlim)
  
  abline(h = c(lower_threshold, middle_threshold), col = "red", lty = 2)
  
  # Add current bed midpoint
  mid_depth <- mean(c(depth[bed_edges[i]], depth[bed_edges[i + 1] - 1]))
  mid_gr <- mean(GR[bed_edges[i]:(bed_edges[i + 1] - 1)], na.rm = TRUE)
  points(mid_depth, mid_gr, col = "blue", pch = 19, cex = 1.2)
  
  dev.off()
}

#dev.off()
graphics.off()
# Create GIF animation 
frames <- list.files(pattern = "frame_\\d+\\.png", full.names = TRUE)
img_list <- image_read(frames)
animation <- image_animate(img_list, fps = 20, loop = 10)  # fps=4 frames/sec
image_write(animation, "litholog_animation.gif")


graphics.off()

#animate litholog with climate ####
# 2 windows:
#window 1: litholog 
#window 2: insolation curve
#setup 

#  1. Read Data 

graphics.off()
# ALT_GR <- getLaskar(sol="insolation")
ALT_GR <- ALT_GR[, c(1, 2)]
depth <- ALT_GR[, 1]
GR <- ALT_GR[, 2]

# img <- readPNG("D:/Phd/documents/Defense/limestone_1.png")  
# img3 <- readPNG("D:/Phd/documents/Defense/marl_1.png")   
# img5 <- readPNG("D:/Phd/documents/Defense/shale_1.png")   
# 

#  2. Lithology Classification from GR 
lower_threshold <- mean(GR) - 0.5 * sd(GR)
middle_threshold <- mean(GR) + sd(GR)

get_lithology <- function(gr) {
  if (is.na(gr)) {
    return(NA)
  } else if (gr < lower_threshold) {
    return("Limestone")
  } else if (gr < middle_threshold) {
    return("Marl")
  } else {
    return("Shale")
  }
}

lithology <- sapply(GR, get_lithology)

#  3. Identify bed boundaries 
lith_factor <- factor(lithology, levels = c("Limestone", "Marl", "Shale"))
bed_changes <- which(diff(as.numeric(lith_factor)) != 0) + 1
bed_edges <- c(1, bed_changes, length(depth) + 1)
bed_type <- lithology[bed_edges[-length(bed_edges)]]

#  4. Lithology-specific widths and colors 
widths <- list(Limestone = 2.5, Marl = 2, Shale = 1.75)
colors <- list(Limestone = "#F4F4E4", Marl = "#E0E0DE", Shale = "#7C7C7C")

#  5. Time steps for animation 
fixed_xlim <- c(0, 300)  
time_steps <- seq(0, 300, by = 1)  # every 10 kyr
n_steps <- length(time_steps)

# set the working directory for the generated images and the gif
setwd("D:/Phd/R/R_results/gif")
getwd()
#graphics.off()


for (k in 1:n_steps) {
  #k <- 20
  current_time <- time_steps[k]
  
  png(sprintf("frame_%03d.png", k), width = 1400, height = 800, res = 150)
  
  #graphics.off()
  # Layout: 3 columns -> [image | top plot] on first row, [image | bottom plot] on second row
  layout(matrix(c(1,2,
                  1,3), nrow = 2, byrow = TRUE),
         widths = c(1, 3), heights = c(2, 1.5))
  
  #  Left side: Image 
  par(mar = c(0,0,0,0), family = "Calibri")
  plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")

  
  
  if (current_time <= max(depth)) {
    idx <- which.min(abs(depth - current_time))  
   val <- GR[idx]
  }
  
  if (val > middle_threshold ){
    img_sel <- img5}
  if ((val < middle_threshold) &(val>lower_threshold)){
    img_sel <- img3}
  if (val<lower_threshold){
    img_sel <- img}

  h <- dim(img_sel)[1]   # image height (pixels)
  w <- dim(img_sel)[2]   # image width  (pixels)
  img_aspect <- w / h
  
  # get panel physical size (in inches)
  panel_size <- par("pin")         # (width, height)
  panel_aspect <- panel_size[1] / panel_size[2]
  
  if (img_aspect > panel_aspect) {
    # image is relatively wider → fit width, pad vertical
    scaled_h <- panel_aspect / img_aspect
    rasterImage(img_sel, 0, (1 - scaled_h)/2, 1, (1 + scaled_h)/2)
  } else {
    # image is relatively taller → fit height, pad horizontal
    scaled_w <- img_aspect / panel_aspect
    rasterImage(img_sel, (1 - scaled_w)/2, 0, (1 + scaled_w)/2, 1)
  }
  
  
  #  Top: Lithology column (beds growing gradually) 
  par(mar = c(4, 4, 2, 8), family = "Calibri")
  plot(NA, ylim = c(0, max(unlist(widths))), xlim = fixed_xlim,
       ylab = "Lithology", xlab = "Time (kyr)",
       main = sprintf("Litholog up to %.0f kyr", current_time),
       yaxt = "n", bty = "n")
  axis(1, las = 1)
  
  for (j in 1:(length(bed_edges)-1)) {
    lith <- bed_type[j]
    width <- widths[[lith]]
    
    top_depth <- depth[bed_edges[j]]
    base_depth <- depth[bed_edges[j + 1] - 1]
    
    if (top_depth < current_time) {
      poly_left <- max(top_depth, 0)
      poly_right <- min(base_depth, current_time)
      
      if (poly_left < poly_right) {
        x <- c(poly_left, poly_right, poly_right, poly_left)
        y <- c(0, 0, width, width)
        polygon(x, y, col = colors[[lith]], border = "black", lwd = 0.5)
      }
    }
  }
  
  legend("right", inset = c(-0.25, 0), legend = names(colors),
         fill = unlist(colors), bty = "n", border = "black", xpd = TRUE)
  
  par(mar = c(4, 4, 2, 8), family = "Calibri")
  valid_gr <- !is.na(GR) & !is.na(depth)
  plot(depth[valid_gr], GR[valid_gr], type = "l", col = "darkgreen", lwd = 1.5,
       xlab = "Time kyr", ylab = "insolation (W/m^2)",
       ylim = range(GR, na.rm = TRUE), xlim = fixed_xlim,
       main = "Insolation")
  abline(h = c(lower_threshold, middle_threshold), col = "red", lty = 2)
  
  
  
  polygon(x=c(fixed_xlim[1]-100,fixed_xlim[2]+100,fixed_xlim[2]+100,fixed_xlim[1]-100), y=c(middle_threshold,middle_threshold,650,650),
          col = as.character(colors[["Shale"]]))
  polygon(x=c(fixed_xlim[1]-100,fixed_xlim[2]+100,fixed_xlim[2]+100,fixed_xlim[1]-100), y=c(middle_threshold,middle_threshold,lower_threshold,lower_threshold),
          col = as.character(colors[["Marl"]]))
  polygon(x=c(fixed_xlim[1]-100,fixed_xlim[2]+100,fixed_xlim[2]+100,fixed_xlim[1]-100), y=c(lower_threshold,lower_threshold,400,400),
          col = as.character(colors[["Limestone"]]))
  
  lines(depth[valid_gr], GR[valid_gr], col = "red", lwd = 2)
  abline(h = c(lower_threshold, middle_threshold), col = "red", lty = 2)
  legend("right", inset = c(-0.25, 0), legend = names(colors),
         fill = unlist(colors), bty = "n", border = "black", xpd = TRUE)
  
  if (current_time <= max(depth)) {
    idx <- which.min(abs(depth - current_time))  
    points(depth[idx], GR[idx], col = "blue", pch = 19, cex = 1.5)
  }
  
  dev.off()
}

# Create GIF animation 
frames <- list.files(pattern = "frame_\\d+\\.png", full.names = TRUE)
img_list <- image_read(frames)
animation <- image_animate(img_list, fps = 10, loop = 10)  # stable speed
image_write(animation, "litholog_animation.gif")


#generate threshold figure ####
#this code generates a pseudo temperature record with an extinction threshold and a 
#temperature perturbation to show the link between the 2.4-Myr ecc cycle 
# and extinction events 
# set the working directory for the generated images and the gif
setwd("D:/Phd/documents/Thesis/SUP_DATA")

# img <- readPNG("D:/Phd/documents/Defense/Silurian_naut_1.png")   
# img2 <- readPNG("D:/Phd/documents/Defense/Silurian_naut_2.png")   
# img3 <- readPNG("D:/Phd/documents/Defense/Silurian_naut_3.png")   

img <- img_naut_1 <- readPNG_from_github("https://raw.githubusercontent.com/stratigraphy/cyclo_strat_gifs/main/Silurian_naut_1.png")
img2 <- img_naut_2 <- readPNG_from_github("https://raw.githubusercontent.com/stratigraphy/cyclo_strat_gifs/main/Silurian_naut_2.png")
img3 <- img_naut_3 <- readPNG_from_github("https://raw.githubusercontent.com/stratigraphy/cyclo_strat_gifs/main/Silurian_naut_3.png")

laskar_11 <- getLaskar(sol="la11")
laskar_04 <- getLaskar(sol="la04")
insolation  <- getLaskar(sol="insolation",verbose=T)
insolation_hilb <- hilbert(insolation)
insolation_sel <- insolation[insolation[, 1] < 5000, ]


laskar_11_2400 <- taner(laskar_11,flow=1/2500,fhigh=1/2000,roll=10^20,xmax=1/500)

graphics.off()
layout.matrix <- matrix(c(1, 2), nrow = 2, ncol = 1)
graphics::layout(mat = layout.matrix,
                 heights = c(1),
                 widths = c(1, 1)) # Ensure equal widths for all 5 plots

par(mar = c(4, 4,2, 1), 
    family = "Calibri", cex.lab = 1.5, cex.axis = 1.3)

xlims <- c(8000,15000)
plot(insolation[,1],insolation[,2],type="l",xlim=xlims,
     xlab="time (kyr)",main=" insolation at 65 deg North (W/m^2)",
     ylab="(W/m^2)")
lines(insolation_hilb[,1],insolation_hilb[,2]+mean(insolation[,2]),col="red",lwd=1.5,)

ylims=c(min(laskar_11[,2])*0.5,max(laskar_11[,2]))

# Copy and transform insolation
insol_clip <- insolation
graphics.off()
#insol_clip[,2] <- abs(insol_clip[,2] - mean(insol_clip[,2]))

# Create thresholded copy (top 5%)
insol_clip_2 <- insol_clip
insol_clip_2[insol_clip_2[,2] < quantile(insol_clip_2[,2], 0.95), 2] <- NA


insol_clip_2<- na.omit(insol_clip_2)
insol_clip_2[,3] <- c(1,insol_clip_2[2:nrow(insol_clip_2),1]-insol_clip_2[1:(nrow(insol_clip_2)-1),1])

polygon_list <- list()
sublist_index <- 1
list_index <- 1

polygon_list[[list_index]] <- matrix(nrow = 0, ncol = ncol(insol_clip_2))  # initialize first sublist

for (i in 1:nrow(insol_clip_2)) {
  current_row <- insol_clip_2[i, ]
  
  if (current_row[3] == 1) {
    # Still in the same segment
    polygon_list[[list_index]] <- rbind(polygon_list[[list_index]], current_row)
  } else {
    # Start a new segment
    list_index <- list_index + 1
    polygon_list[[list_index]] <- matrix(nrow = 0, ncol = ncol(insol_clip_2))  # init new sublist
    polygon_list[[list_index]] <- rbind(polygon_list[[list_index]], current_row)
  }
}

for (i in seq_along(polygon_list)) {
  data <- polygon_list[[i]]
  x <- data[, 1]
  y <- data[, 2]
}



graphics.off()
layout.matrix <- matrix(c(1, 2,3), nrow = 3, ncol = 1)
graphics::layout(mat = layout.matrix,
                 heights = c(1),
                 widths = c(1)) # Ensure equal widths for all 5 plots

par(mar = c(4, 4,2, 1), 
    family = "Calibri", cex.lab = 1.5, cex.axis = 1.3)


xlims <- c(6250,14000)
plot(insolation[,1],insolation[,2],type="l",xlim=xlims,xlab="time (kyr)",main=" Insolation at 65 deg North (W/m^2)",
     ylab="(W/m^2)")
lines(insolation_hilb[,1],insolation_hilb[,2]+mean(insolation[,2]),col="red",lwd=1.5,)

ylims=c(min(laskar_11[,2])*0.5,max(laskar_11[,2]))
plot(laskar_11[,1],laskar_11[,2],type="l",xlim=xlims,col="black",ylim=ylims,lwd=2,xlab="time (kyr)",
     main="Eccentricity solution of Laskar et al. (2011)",ylab="Eccetricity")
lines(laskar_11_2400,col="red")


plot(insolation[,1],insolation[,2],type="l",xlim=xlims,xlab="time (kyr)",main=" insolation threshold",
     ylab="(W/m^2)")

for (i in seq_along(polygon_list)) {
  data <- polygon_list[[i]]
  x <- data[, 1]
  y <- data[, 2]
  polygon(c(x, rev(x)), c(rep(400, length(y)), rep(600, length(y))), col = "red", border = NA)
}
lines(insolation[,1],insolation[,2],col="black")

ylims=c(min(laskar_11[,2])-0.05,max(laskar_11[,2]))


insolation_sel <- insolation

insolation_sel[,2] <- 3000^(insolation[,2]/max(insolation[,2]))
insolation_sel[,2] <- (insolation_sel[,2]/max(insolation[,2]))+18.5
insolation_sel[insolation_sel[,2] <= 19.6, 2] <- 19.6

laskar_11_temp <- cbind(laskar_11[,1],20*laskar_11[,2]+mean(insolation_sel[,2])-0.1)
insolation_sel <- linterp(insolation_sel,1)
insolation_sel <- mwStats(insolation_sel,win=10,ends=TRUE)
insolation_sel <- linterp(insolation_sel,1)

for(i in 1:nrow(insolation_sel)){
  if(insolation_sel[i,2]>laskar_11_temp[i,2]){
    insolation_sel[i,2] <- laskar_11_temp[i,2]
  }
  
}
insolation_sel <- mwStats(insolation_sel,win=5,ends=TRUE)
insolation_sel <- linterp(insolation_sel,1)

graphics.off()
plot(insolation_sel[,1],insolation_sel[,2],type="l",xlim=xlims,
     ylab= "Average Temperature C")
par(new=TRUE)
plot(laskar_11_temp[,1],laskar_11_temp[,2],type="l",xlim=xlims,ylim=range(insolation_sel[,2]),yaxt="n",
     xlab="time (kyr)",lwd=1,col="red",
     main="Eccentricity solution of Laskar et al. (2011)",ylab="Eccetricity")

insolation_pertub <- insolation
insolation_pertub[,2] <- 0
insolation_pertub

insolation_pertub$Insolation <- 0

# Parameters
rise_start <- 9750
rise_peak  <- rise_start-1
fall_end   <- rise_peak-600
temp_increase <- 1.25

# Linear rise: 0 → 2 between 10000 and 9750
rise_mask <- insolation_pertub$Time_ka >= rise_peak &
  insolation_pertub$Time_ka <= rise_start

insolation_pertub$Insolation[rise_mask] <-
  (insolation_pertub$Time_ka[rise_mask] - rise_start) /
  (rise_peak - rise_start) * temp_increase

# Linear fall: 2 → 0 between 9750 and 8000
fall_mask <- insolation_pertub$Time_ka >= fall_end &
  insolation_pertub$Time_ka <= rise_peak

insolation_pertub$Insolation[fall_mask] <-
  (insolation_pertub$Time_ka[fall_mask] - fall_end) /
  (rise_peak - fall_end) * temp_increase
insolation_pertub <- noLow(insolation_pertub,smooth=0.01,genplot=TRUE,output=2)

insolation_pertub_2 <- insolation_pertub
insolation_pertub_2[,2] <- insolation_pertub[,2]+insolation_sel[,2]


insolation_pertub_3 <- cbind(insolation_pertub[,1],insolation_pertub[,2]+mean(insolation_sel[,2]))
insolation_pertub_4 <- cbind(insolation_pertub_3[,1],insolation_pertub_3[,2]+insolation_sel[,2])

insolation_pertub <- insolation
insolation_pertub[,2] <- 0
insolation_pertub

insolation_pertub$Insolation <- 0  


# Parameters
rise_start <- 9150
rise_peak  <- rise_start-1
fall_end   <- rise_peak-600
temp_increase <- 1.25


# Linear rise: 0 → 2 between 10000 and 9750
rise_mask <- insolation_pertub$Time_ka >= rise_peak & 
  insolation_pertub$Time_ka <= rise_start

insolation_pertub$Insolation[rise_mask] <- 
  (insolation_pertub$Time_ka[rise_mask] - rise_start) / 
  (rise_peak - rise_start) * temp_increase

# Linear fall: 2 → 0 between 9750 and 8000
fall_mask <- insolation_pertub$Time_ka >= fall_end & 
  insolation_pertub$Time_ka <= rise_peak

insolation_pertub$Insolation[fall_mask] <- 
  (insolation_pertub$Time_ka[fall_mask] - fall_end) / 
  (rise_peak - fall_end) * temp_increase
insolation_pertub <- noLow(insolation_pertub,smooth=0.01,genplot=TRUE,output=2)


insolation_pertub_2 <- insolation_pertub
insolation_pertub_2[,2] <- insolation_pertub[,2]+insolation_sel[,2]


insol_clip_3 <- insolation_pertub_2
insol_clip_3[insol_clip_3[,2] < 22, 2] <- NA
insol_clip_3<- na.omit(insol_clip_3)
insol_clip_3[,3] <- c(1,insol_clip_3[2:nrow(insol_clip_3),1]-insol_clip_3[1:(nrow(insol_clip_3)-1),1])

polygon_list <- list()
sublist_index <- 1
list_index <- 1

polygon_list[[list_index]] <- matrix(nrow = 0, ncol = ncol(insol_clip_3))  # initialize first sublist

for (i in 1:nrow(insol_clip_3)) {
  current_row <- insol_clip_3[i, ]
  
  if (current_row[3] == 1) {
    # Still in the same segment
    polygon_list[[list_index]] <- rbind(polygon_list[[list_index]], current_row)
  } else {
    # Start a new segment
    list_index <- list_index + 1
    polygon_list[[list_index]] <- matrix(nrow = 0, ncol = ncol(insol_clip_3))  # init new sublist
    polygon_list[[list_index]] <- rbind(polygon_list[[list_index]], current_row)
  }
}

graphics.off()
plot(insolation_pertub_2,type="l",xlim=xlims)
lines(insolation_pertub[,1],insolation_pertub[,2]+mean(insolation_sel[,2]),col="red")
for (i in seq_along(polygon_list)) {
  data <- polygon_list[[i]]
  x <- data[, 1]
  y <- data[, 2]
  polygon(c(x, rev(x)), c(rep(-10, length(y)), rep(100, length(y))), col = "red", border = "red")
}
abline(h=22,lty=3)


#full plot####

dev.off()
graphics.off()

# 4 equal rows
layout(matrix(1:4, nrow=4), heights = rep(1,4))

# Outer margins for global axis labels
par(oma = c(4, 4, 2, 2))  # outer margins (bottom, left, top, right)

xlims <- c(7000,11000)

##  Plot 1: Insolation 
par(mar = c(0, 4, 2, 0), family = "Calibri") # no bottom, keep top margin
plot(insolation[,1], insolation[,2], type="l", xlim=xlims,
     xaxt="n", ylab="(W/m^2)")

polygon(x=c(7300,7700,7700,7300), y=c(0,0,600,600), 
        col = adjustcolor("grey", alpha.f = 0.2))

polygon(x=c(9750,9300,9300,9750), y=c(0,0,600,600), 
        col = adjustcolor("grey", alpha.f = 0.2))
mtext("Insolation at 65°N (W/m^2)", side=3, line=0.5, cex=1,adj=0.35)



##  Plot 2: Temperature 
par(mar = c(0, 4, 0, 0), family = "Calibri") # no bottom/top margins
plot(insolation_sel[,1], insolation_sel[,2], type="l", xlim=xlims, ylim=c(19.5,23),
     xaxt="n", ylab="Temperature (°C)", xlab="")
abline(h=22, lty=3,lwd=2,col="red")
polygon(x=c(7300,7700,7700,7300), y=c(0,0,600,600), 
        col = adjustcolor("grey", alpha.f = 0.2))

polygon(x=c(9750,9300,9300,9750), y=c(0,0,600,600), 
        col = adjustcolor("grey", alpha.f = 0.2))
mtext("Temperature (°C)", side=3, line=-1.5, cex=1,adj=0.35)

insolation_pertub_4 <- cbind(insolation_pertub_3[,1],insolation_pertub_3[,2]+insolation_sel[,2]-mean(insolation_sel[,2]))

##  Plot 3: Perturbation 2 with polygons 
par(mar = c(0, 4, 0, 0), family = "Calibri") # no margins, rely on outer margins
plot(insolation_pertub_4, type="l", xlim=xlims, ylim=c(19.5,23), xlab="",xaxt="n",
     ylab="Temperature (°C)")
lines(insolation_pertub_3[,1], insolation_pertub_3[,2], col="red")

abline(h=22, lty=3,lwd=2,col="red")
lines(insolation_pertub_4)
mtext("Temperature (°C)", side=3, line=-1.5, cex=1,adj=0.35)

##  Shared labels in outer margin 
mtext("time (kyr ago )", side=1, line=2.5, outer=TRUE)


polygon(x=c(7300,7700,7700,7300), y=c(0,0,600,600), 
        col = adjustcolor("grey", alpha.f = 0.2))

polygon(x=c(9750,9300,9300,9750), y=c(0,0,600,600), 
        col = adjustcolor("grey", alpha.f = 0.2))

##  Plot 4: Perturbation with polygons
par(mar = c(0, 4, 0, 0), family = "Calibri") # no bottom/top margins
plot(insolation_pertub_2, type="l", xlim=xlims, ylim=c(19.5,23), xlab="",ylab="Temperature (°C)",)
lines(insolation_pertub[,1], insolation_pertub[,2] + mean(insolation_sel[,2]), col="red")
for (i in seq_along(polygon_list)) {
  data <- polygon_list[[i]]
  x <- data[, 1]; y <- data[, 2]
  polygon(c(x, rev(x)), c(rep(-10, length(y)), rep(100, length(y))),
          col=adjustcolor("red", alpha.f=0.4), border="darkred")
}
abline(h=22, lty=3,lwd=2,col="red")
lines(insolation_pertub_2)
polygon(x=c(7300,7700,7700,7300), y=c(0,0,600,600), 
        col = adjustcolor("grey", alpha.f = 0.2))

polygon(x=c(9750,9300,9300,9750), y=c(0,0,600,600), 
        col = adjustcolor("grey", alpha.f = 0.2))
mtext("Temperature (°C)", side=3, line=-1.5, cex=1,adj=0.35)


### animate it 1 ####
#this code generates a gif 
# 3 windows: 
#left window : cephalopod alive or extinct 
#top right window : insolation record 
#bottom right window : pseudo temperature record with an extinction threshold

#setup 
threshold <- 22
xlims <- c(7000,10500)
time_steps <- seq(from=xlims[2],to=xlims[1],by=-5)
# set the working directory for the generated images and the gif
setwd("D:/Phd/R/R_results/gif_2")
getwd()

#graphics.off()
#dev.off()
n_steps <- length(time_steps)
n_steps
n_cores <- detectCores() - 1   # use all but one core
n_cores <- 4
cl <- makeCluster(n_cores)

clusterExport(cl, c("time_steps", "n_steps", "insolation", 
                    "insolation_sel", "img", "img2", "xlims", "threshold"))  
# export objects needed inside the loop

parLapply(cl, 1:n_steps, function(k) {
#for (k in 1:length(time_steps)) {
#k <- 50 
#graphics.off()
current_time <- time_steps[k]
#img_plot <- image_graph(width = 700, height = 500, res = 96)
png(sprintf("frame_%03d.png", k), width = 2000, height = 1000, res = 200)

layout(matrix(c(1,2,1,3), nrow = 2, byrow = TRUE),
       widths = c(1, 3), heights = c(2, 2))

#  Left side: Image 
par(mar = c(0,0,0,0), family = "Calibri")
plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")

if (current_time <= max(insolation[,1])) {
  idx <- which.min(abs(insolation[,1] - current_time))  
  val <- insolation[idx,1]
}

val <- 1000
if (val > threshold ){
  img_sel <- img2}else{img_sel <- img}

h <- dim(img_sel)[1]   # image height (pixels)
w <- dim(img_sel)[2]   # image width  (pixels)
img_aspect <- w / h

# get panel physical size (in inches)
panel_size <- par("pin")         # (width, height)
panel_aspect <- panel_size[1] / panel_size[2]

if (img_aspect > panel_aspect) {
  # image is relatively wider → fit width, pad vertical
  scaled_h <- panel_aspect / img_aspect
  rasterImage(img_sel, 0, (1 - scaled_h)/2, 1, (1 + scaled_h)/2)
} else {
  # image is relatively taller → fit height, pad horizontal
  scaled_w <- img_aspect / panel_aspect
  rasterImage(img_sel, (1 - scaled_w)/2, 0, (1 + scaled_w)/2, 1)
}

##  Plot 1: Insolation 
par(mar = c(0, 4, 2,2), family = "Calibri") # no bottom, keep top margin
plot(insolation[,1], insolation[,2], type="l", xlim=xlims,
     xaxt="n", ylab="(W/m^2)")

polygon(x=c(7300,7700,7700,7300), y=c(0,0,600,600), 
        col = adjustcolor("grey", alpha.f = 0.2))

polygon(x=c(9750,9300,9300,9750), y=c(0,0,600,600), 
        col = adjustcolor("grey", alpha.f = 0.2))
mtext("Insolation at 65°N (W/m^2)", side=3, line=0.5, cex=1,adj=0.35)


if (current_time <= max(insolation[,1])) {
  idx <- which.min(abs(insolation[,1] - current_time))  
  points(insolation[idx,1], insolation[idx,2],
         col = "blue", pch = 19, cex = 2)
}


##  Plot 2: Temperature 
par(mar = c(4, 4, 0, 2), family = "Calibri") # no bottom/top margins
plot(insolation_sel[,1], insolation_sel[,2], type="l",
     xlim=xlims, ylim=c(19.5,23), ylab="Temperature (°C)", xlab="Time (kyr)")
abline(h=22, lty=3,lwd=2,col="red")
polygon(x=c(7300,7700,7700,7300), y=c(0,0,600,600), 
        col = adjustcolor("grey", alpha.f = 0.2))

polygon(x=c(9750,9300,9300,9750), y=c(0,0,600,600), 
        col = adjustcolor("grey", alpha.f = 0.2))
mtext("Temperature (°C)", side=3, line=-1.5, cex=1,adj=0.35)


if (current_time <= max(insolation_sel[,1])) {
  idx <- which.min(abs(insolation_sel[,1] - current_time))  
  points(insolation_sel[idx,1], insolation_sel[idx,2],
         col = "blue", pch = 19, cex = 2)
}
dev.off()
})
graphics.off()
stopCluster(cl)
gc()
frames <- list.files(pattern = "frame_\\d+\\.png", full.names = TRUE)
img_list <- image_read(frames)
animation <- image_animate(img_list, fps = 20, loop = 0)  # stable speed
image_write(animation, "animation_2.gif")
getwd

### animate it 2 ####
#this code generates a gif 
# 3 windows: 
#left window : cephalopod alive or extinct 
#top right window : insolation record 
#bottom right window : pseudo temperature record with an extinction threshold and 
# a temperature pertubation

#setup 
# set the working directory for the generated images and the gif
setwd("D:/Phd/R/R_results/gif_3")

#graphics.off()
#dev.off()
n_steps <- length(time_steps)

n_cores <- detectCores() - 1   # use all but one core
n_cores <- 4
cl <- makeCluster(n_cores)

clusterExport(cl, c("time_steps", "n_steps", "insolation",
                    "insolation_pertub_3",
                    "insolation_sel", "img", "img2", "xlims", "threshold"))  
# export objects needed inside the loop

parLapply(cl, 1:n_steps, function(k) {
  #for (k in 1:length(time_steps)) {
  current_time <- time_steps[k]
  #img_plot <- image_graph(width = 700, height = 500, res = 96)
  png(sprintf("frame_%03d.png", k), width = 1400, height = 1000, res = 200)
  
  layout(matrix(c(1,2,1,3), nrow = 2, byrow = TRUE),
         widths = c(1, 3), heights = c(2, 2))
  
  #  Left side: Image 
  par(mar = c(0,0,0,0), family = "Calibri")
  plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")
  
  if (current_time <= max(insolation[,1])) {
    idx <- which.min(abs(insolation[,1] - current_time))  
    val <- insolation[idx,1]
  }
  
  val <- 1000
  if (val > threshold ){
    img_sel <- img2}else{img_sel <- img}
  
  h <- dim(img_sel)[1]   # image height (pixels)
  w <- dim(img_sel)[2]   # image width  (pixels)
  img_aspect <- w / h
  
  # get panel physical size (in inches)
  panel_size <- par("pin")         # (width, height)
  panel_aspect <- panel_size[1] / panel_size[2]
  
  if (img_aspect > panel_aspect) {
    # image is relatively wider → fit width, pad vertical
    scaled_h <- panel_aspect / img_aspect
    rasterImage(img_sel, 0, (1 - scaled_h)/2, 1, (1 + scaled_h)/2)
  } else {
    # image is relatively taller → fit height, pad horizontal
    scaled_w <- img_aspect / panel_aspect
    rasterImage(img_sel, (1 - scaled_w)/2, 0, (1 + scaled_w)/2, 1)
  }
  
  ##  Plot 1: Insolation 
  par(mar = c(0, 4, 2,2), family = "Calibri") # no bottom, keep top margin
  plot(insolation[,1], insolation[,2], type="l", xlim=xlims,
       xaxt="n", ylab="(W/m^2)")
  
  polygon(x=c(7300,7700,7700,7300), y=c(0,0,600,600), 
          col = adjustcolor("grey", alpha.f = 0.2))
  
  polygon(x=c(9750,9300,9300,9750), y=c(0,0,600,600), 
          col = adjustcolor("grey", alpha.f = 0.2))
  mtext("Insolation at 65°N (W/m^2)", side=3, line=0.5, cex=1,adj=0.325)
  
  if (current_time <= max(insolation[,1])) {
    idx <- which.min(abs(insolation[,1] - current_time))  
    points(insolation[idx,1], insolation[idx,2],
           col = "blue", pch = 19, cex = 2)}

  insolation_pertub_4 <- cbind(insolation_pertub_3[,1],insolation_pertub_3[,2]+insolation_sel[,2]-mean(insolation_sel[,2]))
  par(mar = c(4, 4, 0, 2), family = "Calibri") # no bottom/top margins
  plot(insolation_pertub_4, type="l", xlim=xlims, ylim=c(19.5,23),
       xlab="time (kyr)",
       ylab="Temperature (°C)")
  lines(insolation_pertub_3[,1], insolation_pertub_3[,2], col="red")
  abline(h=22, lty=3,lwd=2,col="red")
  mtext("Temperature (°C)", side=3, line=-1.5, cex=1,adj=0.35)

  polygon(x=c(7300,7700,7700,7300), y=c(0,0,600,600), 
          col = adjustcolor("grey", alpha.f = 0.2))
  
  polygon(x=c(9750,9300,9300,9750), y=c(0,0,600,600), 
          col = adjustcolor("grey", alpha.f = 0.2))
  
  
  if (current_time <= max(insolation_pertub_4[,1])) {
    idx <- which.min(abs(insolation_pertub_4[,1] - current_time))  
    points(insolation_pertub_4[idx,1], insolation_pertub_4[idx,2],
           col = "blue", pch = 19, cex = 2)}
  dev.off()
})
stopCluster(cl)
frames <- list.files(pattern = "frame_\\d+\\.png", full.names = TRUE)
img_list <- image_read(frames)
animation <- image_animate(img_list, fps = 20, loop = 0)  # stable speed
image_write(animation, "animation_3.gif")



### animate it 3 ####
#this code generates a gif 
# 3 windows: 
#left window : cephalopod alive or extinct 
#top right window : insolation record 
#bottom right window : pseudo temperature record with an extinction threshold and 
# a temperature pertubation

#setup 
threshold <- 22
# set the working directory for the generated images and the gif
setwd("D:/Phd/R/R_results/gif_4")
getwd()

#graphics.off()
#dev.off()
n_cores <- detectCores() - 1   # use all but one core
n_cores <- 4
cl <- makeCluster(n_cores)
event_occured <- 0
clusterExport(cl, c("event_occured","time_steps", "n_steps", "insolation", "insolation_pertub_2","insolation_pertub",
                    "insolation_sel", "img","img2","img3", "xlims", "threshold"))  
# export objects needed inside the loop

parLapply(cl, 1:n_steps, function(k) {
  #for (k in 1:length(time_steps)) {
  current_time <- time_steps[k]
  #img_plot <- image_graph(width = 700, height = 500, res = 96)
  png(sprintf("frame_%03d.png", k), width = 1400, height = 1000, res = 200)
  
  layout(matrix(c(1,2,1,3), nrow = 2, byrow = TRUE),
         widths = c(1, 3), heights = c(2, 2))
  
  #  Left side: Image 
  par(mar = c(0,0,0,0), family = "Calibri")
  plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")
  
  if (current_time <= max(insolation_pertub_2[,1])) {
    idx <- which.min(abs(insolation_pertub_2[,1] - current_time))  
    val <- max(c(insolation_pertub_2[idx+2,2],insolation_pertub_2[idx-2,2]))
  }
  
  if(max(insolation_pertub_2[idx:xlims[2],2])>=22){
    event_occured <- event_occured + 1}
  
  if (event_occured == 0){
    val <- max(c(insolation_pertub_2[idx+2+25,2],insolation_pertub_2[idx-2,2]))}
  
  
  if (val > threshold) {
    img_sel <- img
  } else if (event_occured == 0) {
    img_sel <- img2
  } else {
    img_sel <- img3
  }
  
  
  
  h <- dim(img_sel)[1]   # image height (pixels)
  w <- dim(img_sel)[2]   # image width  (pixels)
  img_aspect <- w / h
  
  # get panel physical size (in inches)
  panel_size <- par("pin")         # (width, height)
  panel_aspect <- panel_size[1] / panel_size[2]
  
  if (img_aspect > panel_aspect) {
    # image is relatively wider → fit width, pad vertical
    scaled_h <- panel_aspect / img_aspect
    rasterImage(img_sel, 0, (1 - scaled_h)/2, 1, (1 + scaled_h)/2)
  } else {
    # image is relatively taller → fit height, pad horizontal
    scaled_w <- img_aspect / panel_aspect
    rasterImage(img_sel, (1 - scaled_w)/2, 0, (1 + scaled_w)/2, 1)
  }
  
  ##  Plot 1: Insolation 
  par(mar = c(0, 4, 2,2), family = "Calibri") # no bottom, keep top margin
  plot(insolation[,1], insolation[,2], type="l", xlim=xlims,
       xaxt="n", ylab="(W/m^2)")
  
  polygon(x=c(7300,7700,7700,7300), y=c(0,0,600,600), 
          col = adjustcolor("grey", alpha.f = 0.2))
  
  polygon(x=c(9750,9300,9300,9750), y=c(0,0,600,600), 
          col = adjustcolor("grey", alpha.f = 0.2))
  mtext("Insolation at 65°N (W/m^2)", side=3, line=0.5, cex=1,adj=0.35)
  
  
  if (current_time <= max(insolation[,1])) {
    idx <- which.min(abs(insolation[,1] - current_time))  
    points(insolation[idx,1], insolation[idx,2],
           col = "blue", pch = 19, cex = 2)
  }
  
  
  ##  Plot 3: Perturbation with polygons 
  par(mar = c(4, 4, 0, 2), family = "Calibri") # no bottom/top margins
  plot(insolation_pertub_2, type="l", xlim=xlims, ylim=c(19.5,23), xlab="Time (kyr)",ylab="Temperature (°C)",)
  lines(insolation_pertub[,1], insolation_pertub[,2] + mean(insolation_sel[,2]), col="red")
  abline(h=22, lty=3,lwd=2,col="red")
  lines(insolation_pertub_2)
  polygon(x=c(7300,7700,7700,7300), y=c(0,0,600,600),
          col = adjustcolor("grey", alpha.f = 0.2))

  polygon(x=c(9750,9300,9300,9750), y=c(0,0,600,600),
          col = adjustcolor("grey", alpha.f = 0.2))
  mtext("Temperature (°C)", side=3, line=-1.5, cex=1,adj=0.35)
  
  
  if (current_time <= max(insolation_pertub_2[,1])) {
    idx <- which.min(abs(insolation_pertub_2[,1] - current_time))  
    points(insolation_pertub_2[idx,1], insolation_pertub_2[idx,2],
           col = "blue", pch = 19, cex = 2)
    if(val > threshold){
      abline(v=insolation_pertub_2[idx,1],col="red",lwd=2)
    }
  }
  dev.off()
})
stopCluster(cl)
frames <- list.files(pattern = "frame_\\d+\\.png", full.names = TRUE)
img_list <- image_read(frames)
img_list <- img_list[1:400]

animation <- image_animate(img_list, fps = 20, loop = 0)  # stable speed
image_write(animation, "animation_4.gif")



#2.4 node plot ####
#this code generates a plot 
#this plot shows the location of 2.4Myr ecc nodes in the insolation curve

#setup 
#graphics.off()
insolation_hilb <- hilbert(insolation)
insolation_hilb_2400 <- taner(insolation_hilb,flow=1/2800,fhigh=1/2000,xmax=1/500)

graphics.off()

xlims <- c(7000,14000)
plot(insolation[,1],insolation[,2],type="l",xlim=xlims,xlab="time (kyr)",main=" insolation at 65 deg North (W/m^2)",
     ylab="(W/m^2)")
lines(insolation_hilb[,1],insolation_hilb[,2]+mean(insolation[,2]),col="darkgrey",lwd=2,)
lines(insolation_hilb_2400[,1],insolation_hilb_2400[,2]+mean(insolation[,2]),col="red")

graphics.off()
plot(insolation,type="l",xlim=c(1500,6500))
lines(insolation_hilb_2400[,1],insolation_hilb_2400[,2]+mean(insolation[,2]),col="red")
lines(insolation_hilb[,1],insolation_hilb[,2]+mean(insolation[,2]),col="darkgrey")

abline(h=quantile(mean(insolation[,2])+insolation_hilb_2400[,2],
                  probs = c(0.25)))

graphics.off()

# Base plot

plot(insolation[,1], insolation[,2], type="l",xlab="Time (kyr)", ylab="(W/m^2)",
     xlim=c(7000,14000))

lines(insolation_hilb_2400[,1],
      insolation_hilb_2400[,2] + mean(insolation[,2]), col="red",lwd=2)
lines(insolation_hilb[,1],insolation_hilb[,2]+mean(insolation[,2]),col="darkgrey",lwd=2
      )

# Threshold
thr <- quantile(mean(insolation[,2]) + insolation_hilb_2400[,2],
                probs = 0.25)

# Get x and y for red curve
x <- insolation_hilb_2400[,1]
y <- insolation_hilb_2400[,2] + mean(insolation[,2])

# Below-threshold mask
below <- y <= thr

# Identify contiguous runs
runs <- rle(below)
ends <- cumsum(runs$lengths)
starts <- c(1, head(ends, -1) + 1)

# y-limits of the plot
ylim <- par("usr")[3:4]

# Draw grey translucent bars
for(i in seq_along(runs$values)) {
  if(runs$values[i]) {
    xleft  <- x[starts[i]]
    xright <- x[ends[i]]
    
    polygon(
      x = c(xleft, xright, xright, xleft),
      y = c(ylim[1], ylim[1], ylim[2], ylim[2]),
      col = rgb(0.5,0.5,0.5,0.3), border = NA
    )
    text(
      x = (xleft + xright)/2, 
      y = ylim[2] + 0.05 * diff(ylim),  # 5% above the top
      labels = "2.4-Myr \necc. node", 
      cex = 0.8, xpd = NA
    )  }
}
