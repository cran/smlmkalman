\name{generate_filaments}
\alias{generate_filaments}
\title{generate_filaments}
\description{
A function which generates simulated SMLM datasets of filamentous objects capable of bundling. Modelled from the dynamics of plant microtubules imaged using DNA-PAINT techniques, but may be applicable to other scenarios.
}
\usage{
generate_filaments(loop_number, field_settings, filament_settings,
single_bundling, bundling_dist, smoothing_settings,
cylinder_settings, optics_settings, visualise_code,
export_data, verbose)
}
\arguments{
  \item{loop_number}{
  Number of times the code should be repeated (to be used in conjunction with export_data to generate multiple .csv datasets)
}
  \item{field_settings}{
  A vector containing three elements: the plot minimum, the plot maximum, and the number of divisions between field values.
  }
  \item{filament_settings}{
  A vector containing four elements: the number of filaments, the ratio of horizontal (0) to vertical (1) filaments, the standard deviation from beginning x/y location to endpoint x/y location, and the angle below which crossing filaments will bundle.
  }
  \item{single_bundling}{
  A TRUE/FALSE statement that prevents bundled filaments from going through further rounds of bundling (if FALSE, might cause unexpected results!)
  }
  \item{bundling_dist}{
  A vector containing two elements: the distance between which two bundled filaments will be separated, and the number of points before crossover that should be made to be similarly parallel.
  }
  \item{smoothing_settings}{
  A vector containing three elements: the probability of a point in a filament being randomised away from initial position, the radius about the initial position that randomised points may move within, and the number of rounds of Laplacian smoothing.
  }
  \item{cylinder_settings}{
  A vector containing three elements: the number of points to be distributed in a cylinder about the generated splines (per filament), the width of the cylinder, and the "falloff" value (if 1 or greater, then the n-photon value is higher at the edges of the cylinder).
  }
  \item{optics_settings}{
  A vector containing three elements: a multiplier for whole-field n-photon (gain), the randomisation of the n-photon value of each point in the field, and the optical resampling radius of each point in the field (optical error).
  }
  \item{visualise_code}{
  A vector containing three elements: a TRUE/FALSE statement that dictates whether the simulation process should be plotted in real-time, a TRUE/FALSE statement that dictates whether a set of debug graphs should be exported to the active working directory, and the number of seconds that should be waited between blocks of code to allow the user to view the real-time graphs (set to 0 if visualise_code[1] is FALSE to avoid unnecessarily long wait times)
  }
  \item{export_data}{
  A vector containing two elements: a TRUE/FALSE statement that dictates whether, on each iteration of loop_number, the dataset should be saved to the working directory as a .csv file, and the file prefix if export_data[1] is TRUE.
  }
  \item{verbose}{
  A TRUE/FALSE statement that dictates whether the current place in the function loop should be printed to the terminal.
  }
}

\value{
A dataframe with 4 columns: the ID of each datapoint (which filament it belongs to), the x-position, the y-position, and the n-photon (simulated photonic emission), which is dependent upon the location, being boosted at the edges of the filament.
}

\author{
Andrew Buist
}

\examples{
  \dontshow{
  #Generate dataset
data = generate_filaments(loop_number = 1,
                          field_settings = c(0, 100, 1),
                          filament_settings = c(2, 0, 1, 1),
                          single_bundling = TRUE,
                          bundling_dist = c(3,10),
                          smoothing_settings = c(0.3, 0, 1),
                          cylinder_settings = c(500, 2, 1),
                          optics_settings = c(1, 0.3, 0.5),
                          visualise_code = c(FALSE, FALSE, 0),
                          export_data = c(FALSE, "data_"),
                          verbose = TRUE)

#Plot dataset coloured by ID, and opacity = n-photon
library(scales)
plot(data[,2:3], col = alpha(data[,1], alpha = data[,4]), pch = 16, cex = 0.3)
  }
  \donttest{
#Generate dataset
data = generate_filaments(loop_number = 1,
                          field_settings = c(0, 100, 1),
                          filament_settings = c(5, 0.3, 85, 17),
                          single_bundling = TRUE,
                          bundling_dist = c(3,10),
                          smoothing_settings = c(0.3, 5, 25),
                          cylinder_settings = c(500, 2, 1),
                          optics_settings = c(1, 0.3, 0.5),
                          visualise_code = c(FALSE, FALSE, 0),
                          export_data = c(FALSE, "data_"),
                          verbose = TRUE)

#Plot dataset coloured by ID, and opacity = n-photon
library(scales)
plot(data[,2:3], col = alpha(data[,1], alpha = data[,4]), pch = 16, cex = 0.3)
  }
}
