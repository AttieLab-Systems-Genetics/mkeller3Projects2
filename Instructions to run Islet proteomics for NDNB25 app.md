\# Run these lines in R console

\# Install required packages if they don\'t have them

required_pkgs \<- c(\"shiny\", \"tidyverse\", \"plotly\", \"DT\")

new_pkgs \<- required_pkgs\[!(required_pkgs %in%
installed.packages()\[,\"Package\"\])\]

if(length(new_pkgs)) install.packages(new_pkgs)

\# Set working directory to the app folder on mkeller3 research drive

\# You will have to adjust your path designation in Apple to what I'm
using in Windows

\# This folder contains the data files used by the app, as well as the
source code

setwd(\"W:/General/Projects2/Lloyd Smith/islet_app\")

\# Run the app

shiny::runApp()
