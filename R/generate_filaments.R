#####
#Code written by Andrew Buist, updated 16/05/22
#A SMLM DNA-PAINT-like dataset generator for use in conjunction with crescent_kf() to produce
#categorised training datasets to optimise hyperparameters, for subsequent real-data segmentation
#####

generate_filaments = function(loop_number = 1,
                              field_settings = c(0, 100, 1),
                              filament_settings = c(10, 0.3, 85, 17),
                              single_bundling = TRUE,
                              bundling_dist = c(3,10),
                              smoothing_settings = c(0.5, 10, 100),
                              cylinder_settings = c(1000, 1, 1),
                              optics_settings = c(1, 0.3, 0.3),
                              visualise_code = c(FALSE, FALSE, 0),
                              export_data = c(FALSE, "data_"),
                              verbose = TRUE){
#library(stats)
#library(spdep)
#library(pracma)
#library(scales)
#library(truncnorm)

oldpar = par(no.readonly = TRUE)
on.exit(par(oldpar))
start.time = Sys.time()

##Do you want to see the process? How long between graphs (in s)?
visualise.code = visualise_code[1]
graph.wait = visualise_code[3]

#Do you want debug graphs to be exported?
debug.graphs = visualise_code[2]

#Do you only want bundled pairs? (MTs may be unstable if FALSE!)
singlebundling = single_bundling

#Do you want to export your data? How many repeat exports are required?
export.data = export_data[1]
loopnum = loop_number

###WORM-LIKE CHAIN SETTINGS###
#number of mts
rep = filament_settings[1]
#number of divisions per field unit
resolution = field_settings[3]
#minimum field value
plotmin = field_settings[1]
#maximum field value
plotmax = field_settings[2]
#skew informs the rnum of average microtubule orientation - 0 = all  horizontal, 1 = all vertical, 0.5 = random
skew = filament_settings[2]
#trajectory informs the rtruncnorm function which selects the angle at which lines cross the plane
#avoid = 0, as code is dependent on finding crosses and none exist if = 0
trajectory = filament_settings[3]
#stochasticity informs the chance of random coordinates - 0 = (linear) no chance, 1 = random
stochasticity = smoothing_settings[1]
#sharpness informs the radius within which a point may randomly travel (sd)
sharpness = smoothing_settings[2]
#reduction.epochs = rounds of smoothing
smoothing = smoothing_settings[3]
#max follow angle
max.follow.angle = filament_settings[4]
#minimum bundling distance
min_dist = bundling_dist[1]
#back-propagation for intersect smoothing
back_prop = bundling_dist[2]

###CYLINDER SETTINGS###
#cylinder resolution
cylres = cylinder_settings[1]
#cylinder width
cylwidth = cylinder_settings[2]
#divisor for nphoton
falloff = cylinder_settings[3]
#max nphoton
intensity = optics_settings[1]
#randomisation in nphoton
rintensity = optics_settings[2]
#optical resampling error
opticalrand = optics_settings[3]

###WORM ENGINE###
divisions = ((plotmax-plotmin)*resolution)
reduction.epochs = smoothing*resolution

for(loops in 1:loopnum){
if(verbose == TRUE){
cat("\014")
cat(paste("Current Loop: ", loops, sep = ""))
cat(paste("\n", "Linears..."))
}

#storage
data.master = as.data.frame(matrix(data = NA, ncol = 3, nrow = 0))
data.master.smooth = as.data.frame(matrix(data = NA, ncol = 3, nrow = 0))
colnames(data.master) = c("x", "y", "ID")

#point distances
point.distances = as.data.frame(matrix(data = NA, nrow = reduction.epochs, ncol = rep))

#follow-me
follow.me.angles = as.data.frame(matrix(data = NA, nrow = rep, ncol = rep))
follow.me.intersect.x = as.data.frame(matrix(data = NA, nrow = rep, ncol = rep))
follow.me.intersect.y = as.data.frame(matrix(data = NA, nrow = rep, ncol = rep))

#linear equations for bundling
linears = as.data.frame(matrix(NA,ncol = 3, nrow = rep))
names(linears) = c("m", "c", "h/v")

#coerce diagonal pairwise-matrices to x and y for plotting
crosses = as.data.frame(matrix(data = NA, ncol = 3, nrow = 0))
colnames(crosses) = c("x", "y", "acute angle")

#for loop, generates each microtubule independently, then adds this to the .master matrix (ALL MT CODE MUST BE INSIDE HERE)
for(i in 1:rep){
  data = as.data.frame(matrix(data = NA, ncol = 3, nrow = divisions))
  colnames(data) = c("x", "y", "id")
  data[,3] = i
  #skew.chance = random number 0:1, test against skew
  skew.chance = runif(1, min = 0, max = 1)

  #is the random number (skew.chance) greater than the threshold imposed by skew?
  # skew = 0, always the case - horizontal
  #skew = 1, never the case - vertical
  if(skew.chance > skew){
    #horizontal positioning
    data[1,1] = plotmin
    data[divisions,1] = plotmax
    data[1,2] = runif(1, min = plotmin, max = plotmax)
    data[divisions,2] = rtruncnorm(1, a = plotmin, b = plotmax, mean = data[1,2], sd = trajectory)
    linears[i,3] = 0
  } else {
    #vertical positioning
    data[1,1] = runif(1, min = plotmin, max = plotmax)
    data[divisions,1] = rtruncnorm(1, a = plotmin, b = plotmax, mean = data[1,1], sd = trajectory)
    data[1,2] = plotmin
    data[divisions,2] = plotmax
    linears[i,3] = 1
  }
  if(visualise.code == TRUE){
    plot(c(plotmin,plotmax),c(plotmin,plotmax), type = "n", xlab = "Arbitrary Units", ylab = "Arbitrary Units")
    points(data[,1:2])
    Sys.sleep(graph.wait)
  }

  #populate space between anchor points with (n = divisions) equally spaced points
  for(j in 2:(divisions-1)){
      data[j,1] = data[j-1,1] - (data[1,1]-data[(divisions),1])/(divisions-1)
      data[j,2] = data[j-1,2] - (data[1,2]-data[(divisions),2])/(divisions-1)
  }
  if(visualise.code == TRUE){
    plot(c(plotmin,plotmax),c(plotmin,plotmax), type = "n", xlab = "Arbitrary Units", ylab = "Arbitrary Units")
    points(data[,1:2])
    Sys.sleep(graph.wait)
  }

  #follow-me function
  #to resolve the need to calculate the linear moment of a = b AND b = a,
  #using a triangular number indexing system
  #i.e. when rep = 1, nothing is eval'd, = 2 :1, = 3 :3, = 4 :6 ...
  if(linears[i,3] == 0){
    #horizontal
    #y = mx + c
    linears[i,1] = (data[divisions,2]-data[1,2])/(data[divisions,1]-data[1,1])
    linears[i,2] = data[1,2]-(linears[i,1]*data[1,1])
  } else {
    #vertical
    #x = my + c
    linears[i,1] = (data[divisions,1]-data[1,1])/(data[divisions,2]-data[1,2])
    linears[i,2] = data[1,1]-(linears[i,1]*data[1,2])
  }

  if(i != 1){
      for(j in 1:(i-1)){
        #solve linear algebra to find x and y intercept
        #HH - tan(a) = |m1-m2/1+m1m2|
        #HV - tan(a) = |m1-(1/m2)/1+m1(1/m2)|
        #VH - tan(a) = |(1/m1)-m2/1+(1/m1)m2|
        #VV - tan(a) = |m1-m2/1+m1m2| (both are x = my + c, so standard)
        if(linears[i,3] == 0 & linears[j,3] == 0){
          #HH
          x.cross = (linears[i,2]-linears[j,2])/(linears[j,1]-linears[i,1])
          y.cross = (linears[i,1]*x.cross) + linears[i,2]

          angles.cross = atan(abs((linears[i,1]-linears[j,1])/(1+(linears[i,1]*linears[j,1]))))*(180/pi)
        } else if(linears[i,3] == 0 & linears[j,3] == 1){
          #HV
          x.cross = (-(linears[j,1]*linears[i,2])-linears[j,2])/((linears[i,1]*linears[j,1])-1)
          y.cross = (linears[i,1]*x.cross) + linears[i,2]

          angles.cross = atan(abs((linears[i,1]-(1/linears[j,1]))/(1+(linears[i,1]*(1/linears[j,1])))))*(180/pi)
        } else if(linears[i,3] == 1 & linears[j,3] == 0){
          #VH
          x.cross = (-(linears[i,1]*linears[j,2])-linears[i,2])/((linears[j,1]*linears[i,1])-1)
          y.cross = (linears[j,1]*x.cross) + linears[j,2]

          angles.cross = atan(abs(((1/linears[i,1])-linears[j,1])/(1+((1/linears[i,1])*linears[j,1]))))*(180/pi)
        } else {
          #VV
          y.cross = (linears[i,2]-linears[j,2])/(linears[j,1]-linears[i,1])
          x.cross = (linears[i,1]*y.cross) + linears[i,2]

          angles.cross = atan(abs((linears[i,1]-linears[j,1])/(1+(linears[i,1]*linears[j,1]))))*(180/pi)
        }

        follow.me.intersect.x[i,j] = x.cross
        follow.me.intersect.y[i,j] = y.cross

        #derive angles at crosspoint
        follow.me.angles[i,j] = angles.cross
      }
    }
  data.master = rbind(data.master, data)
}

if(verbose == TRUE){
cat(" Done!")
cat(paste("\n", "Bundling..."))
}

if(debug.graphs == TRUE){
  png(filename = "crosses_debug.png", height = 1000, width = 1000)
  par(cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, mar = c(5,5,0,0))
  plot(c(plotmin,plotmax),c(plotmin,plotmax), type = "n", xlab = "Arbitrary Units", ylab = "Arbitrary Units")
  #establish colour spacing
  cols = rainbow(max(data.master[,3]), start = 0, end = 5/6)
  #plot points
  points(data.master[,1:2], col = cols[data.master[,3]], pch = 16, cex = 1)
    points(x = unlist(apply(follow.me.intersect.x, 2, na.omit)), y = unlist(apply(follow.me.intersect.y, 2, na.omit)), pch = 13, cex = 3)
  text(x = unlist(apply(follow.me.intersect.x, 2, na.omit)), y = unlist(apply(follow.me.intersect.y, 2, na.omit))+5, labels = signif(unlist(apply(follow.me.angles, 2, na.omit)), digits = 3))
  legend('topright', legend = c(paste("Skew = ", skew), paste("Divisions = ", divisions), paste("MT Number = ", rep)), cex = 1.5)
  dev.off()
}

pair.bundle = 0

#follow-me functions continued
for(i in 2:rep){
  for(j in 1:(i-1)){
    c.x = follow.me.intersect.x[i,j]
    c.y = follow.me.intersect.y[i,j]
    c.angle = follow.me.angles[i,j]

    #this section creates indexes for the member points of each line to be bundled
    #this 1st if() function dictates whether bundling occurs before or after the intersect
    #as in prior if-block, these represent HH, HV, VH, VV scenarios
    if(linears[i,3] == 0){
      #H?
      a.c = c.x
      a.priors = (which.min(abs(data.master[which(data.master[,3] == i),1] - a.c))+((i-1)*divisions))
    } else {
      #V?
      a.c = c.y
      a.priors = (which.min(abs(data.master[which(data.master[,3] == i),2] - a.c))+((i-1)*divisions))
    }
    if(linears[j,3] == 0){
      #?H
      b.c = c.x
      b.priors = (which.min(abs(data.master[which(data.master[,3] == j),1] - b.c))+((j-1)*divisions))
    } else {
      #?V
      b.c = c.y
      b.priors = (which.min(abs(data.master[which(data.master[,3] == j),2] - b.c))+((j-1)*divisions))
    }

    #informs the min_dist usage case (l:302 & 319)
    flip_result = 0
    #"toss a coin" in the 'if' to decide if bundling occurs before or after crossover
    if(runif(1,0,1) < 0.5){
      a.back = a.priors - ((i-1)*divisions)
      b.back = b.priors - ((j-1)*divisions)
      #if the length from the intersect to start/end is greater than the backprop value, it is overwritten by this value
      if(a.back > back_prop){a.back = back_prop}
      if(b.back > back_prop){b.back = back_prop}

      a.index = (a.priors - a.back):(i*divisions)
      b.index = (b.priors - b.back):(j*divisions)
      flip_result = 0
    } else {
      a.back = (i*divisions) - a.priors
      b.back = (j*divisions) - b.priors

      if(a.back > back_prop){a.back = back_prop}
      if(b.back > back_prop){b.back = back_prop}

      a.index = (((i-1)*divisions)+1):(a.priors + a.back)
      b.index = (((j-1)*divisions)+1):(b.priors + b.back)
      flip_result = 1
    }

    if(c.x > plotmin &
       c.x < plotmax &
       c.y > plotmin &
       c.y < plotmax){
      crosses.append = cbind(follow.me.intersect.x[i,j], follow.me.intersect.y[i,j], follow.me.angles[i,j])
      crosses = rbind(crosses, crosses.append)

      #if angle below threshold *AND* bundling in same "plane" (H/V,V/H discounted, important)
      if(c.angle < max.follow.angle & linears[i,3] == linears[j,3]){
        pair.bundle = pair.bundle+1

        m1 = linears[i,1]
        m2 = linears[j,1]

        c1 = linears[i,2]
        c2 = linears[j,2]

        mu = (m1+m2)/2
        cu = (c1+c2)/2

        if((linears[i,3] == 0) & (linears[j,3] == 0)){
          #the lower slope will always need displacement added & the higher slope, subtracted
          #this is inverted if the MTs bundle in the other direction (when flip_result == 1), so signs flipped accordingly
          if(flip_result == 0){
            if(m1 < m2){
              data.master[a.index,2] = mu*data.master[a.index,1] + cu + min_dist
              data.master[b.index,2] = mu*data.master[b.index,1] + cu - min_dist
            } else {
              data.master[a.index,2] = mu*data.master[a.index,1] + cu - min_dist
              data.master[b.index,2] = mu*data.master[b.index,1] + cu + min_dist
            }
          } else {
            if(m1 < m2){
              data.master[a.index,2] = mu*data.master[a.index,1] + cu - min_dist
              data.master[b.index,2] = mu*data.master[b.index,1] + cu + min_dist
            } else {
              data.master[a.index,2] = mu*data.master[a.index,1] + cu + min_dist
              data.master[b.index,2] = mu*data.master[b.index,1] + cu - min_dist
            }
          }
        } else {
          if(flip_result == 0){
            if(m1 < m2){
              data.master[a.index,1] = mu*data.master[a.index,2] + cu + min_dist
              data.master[b.index,1] = mu*data.master[b.index,2] + cu - min_dist
            } else {
              data.master[a.index,1] = mu*data.master[a.index,2] + cu - min_dist
              data.master[b.index,1] = mu*data.master[b.index,2] + cu + min_dist
            }
          } else {
            if(m1 < m2){
              data.master[a.index,1] = mu*data.master[a.index,2] + cu - min_dist
              data.master[b.index,1] = mu*data.master[b.index,2] + cu + min_dist
            } else {
              data.master[a.index,1] = mu*data.master[a.index,2] + cu + min_dist
              data.master[b.index,1] = mu*data.master[b.index,2] + cu - min_dist
            }
          }
        }

        #force bundled pairs outside the valid x/y range
        if(singlebundling == TRUE){
        follow.me.intersect.x[,i] = follow.me.intersect.y[,i] =
          follow.me.intersect.x[i,] = follow.me.intersect.y[i,] =
          follow.me.intersect.x[,j] = follow.me.intersect.y[,j] =
          follow.me.intersect.x[j,] = follow.me.intersect.y[j,] =
          #value always too low, numeric value = on which run bundling occurred
          plotmin-pair.bundle
        }

        if(visualise.code == TRUE){
          plot(c(plotmin,plotmax),c(plotmin,plotmax), type = "n", xlab = "Arbitrary Units", ylab = "Arbitrary Units")
          points(data.master[,1:2])
          Sys.sleep(graph.wait)
        }
      }
    }
  }
}

if(verbose == TRUE){
cat(" Done!")
cat(paste("\n", "Randomisation & Smoothing..."))
}

if(debug.graphs == TRUE){
  png(filename = "bundling_debug.png", height = 1000, width = 1000)
  par(cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, mar = c(5,5,0,0))
  plot(c(plotmin,plotmax),c(plotmin,plotmax), type = "n", xlab = "Arbitrary Units", ylab = "Arbitrary Units")
  #establish colour spacing
  cols = rainbow(max(data.master[,3]), start = 0, end = 5/6)
  #plot points
  for(i in 1:rep){
    lines(data.master[which(data.master[,3] == i),1:2], col = cols[i], lwd = 3)
  }
  legend('topright', legend = c(paste("Bundling Angle Upper Bound = ", max.follow.angle)), cex = 1.5)
  dev.off()
}

for(i in 1:rep){
  data = data.master[(((i-1)*divisions)+1):(i*divisions),]
  #randomise points
  for(j in 2:(divisions-1)){
    #stochasticity.chance = random number 0:1, test against stochasticity
    stochasticity.chance = runif(1, min = 0, max = 1)
    #is the random number (stochasticity.chance) greater than the threshold imposed by stochasticity?
    #stochasticity = 0, always the case
    #stochasticity = 1, never the case
    if(stochasticity.chance > stochasticity){
      data[j,1] = data[j,1]
      data[j,2] = data[j,2]
      } else {
        data[j,1] = rnorm(1, mean = data[j,1], sd = sharpness)
        data[j,2] = rnorm(1, mean = data[j,2], sd = sharpness)
      }
  }

  #epoch reduction smoothing
  for(f in 1:reduction.epochs){
    point.distances.hold = rep(0, divisions-1)
    data.hold1 = data.hold2 = data[(nrow(data)-divisions):(nrow(data)),]
   for(j in 2:(divisions-1)){
     x.av = (data.hold1[j-1,1]+ data.hold1[j,1]+ data.hold1[j+1,1])/3
     y.av = (data.hold1[j-1,2]+ data.hold1[j,2]+ data.hold1[j+1,2])/3
     data.hold2[j,1] = x.av
     data.hold2[j,2] = y.av
   }
    for(h in 2:divisions){
      point.distances.hold[h] = sqrt((data.hold2[h,1] - data.hold2[(h-1),1])^2 + (data.hold2[h,2] - data.hold2[(h-1),2])^2)
    }
    #change to max or mean to see different key stats
    point.distances[f,i] = max(point.distances.hold)
    data = data.hold2
  }
    data.master.smooth = rbind(data.master.smooth, data)
    if(visualise.code == TRUE){
      plot(c(plotmin,plotmax),c(plotmin,plotmax), type = "n", xlab = "Arbitrary Units", ylab = "Arbitrary Units")
      points(data.master.smooth[,1:2])
    }

    #reduces bunching
    if(linears[i,3] == 0){
      data.master.smooth[which(data.master.smooth[,3] == i),1] = linspace(plotmin, plotmax, divisions)
    } else {
      data.master.smooth[which(data.master.smooth[,3] == i),2] = linspace(plotmin, plotmax, divisions)
    }
}

if(verbose == TRUE){
cat(" Done!")
cat(paste("\n", "Cylinders & Optics..."))
}

if(debug.graphs == TRUE){
  #plots
  png(filename = "splines_smooth.png", width = 2000, height = 1000)
  par(mfrow = c(1,2), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, mar = c(5,5,0,0))
  #empty plot
  plot(c(plotmin,plotmax),c(plotmin,plotmax), type = "n", xlab = "Arbitrary Units", ylab = "Arbitrary Units")
  #establish colour spacing
  cols = rainbow(max(data.master.smooth[,3]), start = 0, end = 5/6)
  #plot points
  points(data.master.smooth[,1:2], col = cols[data.master.smooth[,3]], pch = 16, cex = 1)

  plot(c(0,reduction.epochs), c(min(point.distances),max(point.distances)), type = "n", xlab = "Smoothing Epochs", ylab = "Maximum Point-Point Distance (LOG)", log = "y")
  #repeat point plotting in common colours
  for(i in 1:rep){
    points(x = 1:reduction.epochs, y = point.distances[1:reduction.epochs,i], col = cols[i], pch = 16, cex = 1)
  }
  dev.off()
}

###CYLINDER ENGINE###
cylinders = as.data.frame(matrix(data = NA, nrow = rep*cylres, ncol = 4))
names(cylinders) = c("id", "x", "y", "nphoton")

for(i in 1:rep){
  k = ((i-1)*cylres)+1

  for(j in k:(cylres*i)){
    cylinders[j,1] = i
    cylinders[j,2] = ((j+1)-k)

    cylinders[j,3] = (sin(j))

    cylinders[j,4] = abs(((cylinders[j,3]*intensity))+rnorm(1, mean = 0, sd = 1*rintensity))/falloff

    cylinders[j,3] = cylinders[j,3]*cylwidth
  }
}


for(i in 1:rep){
  if(linears[i,3] == 0){
    #if path is horizontal, treat x = x, y = y
    input.array = data.master.smooth[which(data.master.smooth[,3] == i),1:2]
    fun = splinefun(input.array, method = "monoH.FC")

    x.region = input.array[,1]
    x.region = linspace(x.region[1], x.region[divisions], cylres)
    height.displace = fun(x.region)
    cylangles = atan(fun(x.region, deriv = 1))

    cylinders[which(cylinders[,1] == i), 2] = rep(0, cylres)

    #Rotation() can only handle 1 angle argument, so must loop through each successive
    for(j in 1:cylres){
      cylinders[which(cylinders[,1] == i)[j], 2:3] = Rotation(cylinders[which(cylinders[,1] == i)[j], 2:3], cylangles[j])
    }
    cylinders[which(cylinders[,1] == i), 2] = cylinders[which(cylinders[,1] == i), 2] + x.region
    cylinders[which(cylinders[,1] == i), 3] = cylinders[which(cylinders[,1] == i), 3] + height.displace
  } else {
    #else, path must be vertical, treat x = y, y = x
    input.array = data.master.smooth[which(data.master.smooth[,3] == i),2:1]
    fun = splinefun(input.array, method = "monoH.FC")

    x.region = input.array[,1]
    x.region = linspace(x.region[1], x.region[divisions], cylres)
    height.displace = fun(x.region)
    cylangles = atan(fun(x.region, deriv = 1))

    cylinders[which(cylinders[,1] == i), 2] = rep(0, cylres)

    #Rotation() can only handle 1 angle argument, so must loop through each successive
    for(j in 1:cylres){
      cylinders[which(cylinders[,1] == i)[j], 2:3] = Rotation(cylinders[which(cylinders[,1] == i)[j], 2:3], cylangles[j])
    }

    cylinders[which(cylinders[,1] == i), 2] = cylinders[which(cylinders[,1] == i), 2] + x.region
    cylinders[which(cylinders[,1] == i), 3] = cylinders[which(cylinders[,1] == i), 3] + height.displace
    #this is really the only change that makes this work: calculate everything relative to x, then flip at the end of the cycle
    cylinders[which(cylinders[,1] == i), 2:3] = cylinders[which(cylinders[,1] == i), 3:2]
  }

  if(visualise.code == TRUE){
    plot(fun(input.array[,1]), type = "l", main = paste("Spline: #", i))
    Sys.sleep(graph.wait)
  }
}

if(debug.graphs == TRUE){
  png(filename = "cylinders_transposed.png", height = 1000, width = 1000)
  plot(cylinders[,2:3], cex = 0.5, col = alpha(cols[cylinders[,1]], alpha = cylinders[,4]))
  dev.off()
}

#resampling - optics
for(i in 1:nrow(cylinders)){
  #as x and y are treated individually, one might be in a ratio to the other, i.e. an elliptical PSF
  cylinders[i,2] = rnorm(1, mean = cylinders[i,2], sd = opticalrand)
  cylinders[i,3] = rnorm(1, mean = cylinders[i,3], sd = opticalrand)
}

if(verbose == TRUE){
cat(" Done!")
}

if(debug.graphs == TRUE){
  png(filename = "cylinders_perturbed.png", height = 1000, width = 1000)
  par(bg = "black", fg = "white", col.axis = "white", col.lab = "white")
  plot(cylinders[,2:3], col = alpha("white", alpha = cylinders[,4]), cex = cylinders[,4]*.5, pch = 16)
  dev.off()
}

  if(export.data == TRUE){
    write.csv(x = cylinders, file = paste(export_data[2], loops, ".csv", sep = ""), row.names = FALSE)
  }
  par(bg = "black", fg = "white", col.axis = "white", col.lab = "white")
  plot(cylinders[,2:3], col = alpha("white", (cylinders[,4])/2), pch = 16, cex = (cylinders[,4])/2)
}
end.time = Sys.time()
time.taken = difftime(end.time, start.time, units = "secs")

if(verbose == TRUE){
cat(paste("\n", "Duration: ", signif(time.taken, digits = 6), " seconds!"))
}

if(export.data[1] == FALSE){
  cylinders
}
}
