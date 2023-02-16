#Spliner for old bioxpert
#

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #Setsworking directory to script location on hard drive (only in RStudio)

###########################
#Install required Packages#
###########################

req.packages <- c("tidyverse","writexl","readxl", "zoo", "scales", "ggpubr")
new.packages <- req.packages[!(req.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) {
  install.packages(new.packages)
}

lapply(req.packages, require, character.only = TRUE)

rm(list = ls())

###############
#Read meta_data#
###############



meta_data <- read_xlsx((list.files(pattern = "*MetaData.xlsx"))[1], sheet = "Metadata", col_names = FALSE) %>% 
  filter (!is.na(...1))



project_id  <- meta_data[2][meta_data[1] == "Project ID"]
exp_nr      <- meta_data[2][meta_data[1] == "Bioreactor experiment nr"]
controller_ids <- as.character(meta_data[meta_data[1] == "Controller ID",][-1])
experiment_ids <- paste(exp_nr,controller_ids, sep ="_")
names(meta_data) <- c("parameter",experiment_ids)


alkali_use <- as.numeric(meta_data[meta_data[1] == "Alkali, bottle weight before exp.",][-1]) - 
  as.numeric(meta_data[meta_data[1] == "Alkali, bottle weight after exp.",][-1])

acid_use <- as.numeric(meta_data[meta_data[1] == "Acid, bottle weight before exp.",][-1]) - 
  as.numeric(meta_data[meta_data[1] == "Acid, bottle weight after exp.",][-1])

alkali_conc <- as.numeric(meta_data[meta_data[1] == "Alkali, N",][-1])
acid_conc <- as.numeric(meta_data[meta_data[1] == "Acid, N",][-1])


inok_time <- as.numeric(meta_data[meta_data[1] == "Inoculation",][-1])
initial_od <- as.numeric(meta_data[meta_data[1] == "initial_OD_manual",][-1])
final_od <- as.numeric(meta_data[meta_data[1] == "final_OD_manual",][-1])
time_final_od <- as.numeric(meta_data[meta_data[1] == "time_final_OD",][-1])
reactor_volume <- as.numeric(meta_data[meta_data[1] == "Bioreactor volume, mL",][-1])


raw_files <- as.character(meta_data[meta_data[1] == "RAW file name",][-1])

out_filename = paste(
  project_id,
  exp_nr,
  "BioreactorResults_R.xlsx",
  sep = "_"
)


######################################
#File parsing, converts time to [h]  #
######################################
ferm_files  <- list() # create empty list for storing data frames
for (i in raw_files) {
  ferm_files[[i]] <- read_table(       #This function reads raw bioxpert tables into R data frames and saves it to ferm_files, which is a list of data frames
    paste0("RAW/",i),
    col_types = cols(
      Time   = col_character(),
      pH     = col_double(),
      Vc     = col_double(),
      T      = col_double(),
      aOD    = col_double(),
      Stir   = col_double(),
      Redox  = col_double(),
      Alkali = col_double(),
      Acid   = col_double(),
      In     = col_double(),
      pmp4   = col_double()
    )
  )[-1,]
  ferm_files[[i]]$Time <-  sapply(                          #This function converts time from hh:mm to [h]
    strsplit(as.character(ferm_files[[i]]$Time), ":"),
    function(x) {
      x <- as.numeric(x)
      x[1] + x[2] / 60
    }
  )
  ferm_files[[i]] <- as.data.frame(ferm_files[[i]]) %>%
    filter(!is.na(Time))      %>%   #Remove NA values
    mutate(Experiment = experiment_ids[which (raw_files == i)])    %>%   # Adds exp number name from experiment_id
    filter(c(diff(Acid),NA) <  2) %>%  #Removes garbage lines
    filter(c(diff(Acid),NA) > -2) %>%  #Removes garbage lines %>% 
    filter(!is.na(Alkali)) %>% 
    filter(!is.na(Acid)) %>% 
    select(Experiment, Time, pH, Vc, T, aOD, Stir, Redox, Alkali, Acid, In, pmp4) # Selects relevant columns for output file
}

#Rescale Alkali to bottle weight, substract Alkali used for neutralization:
for (i in 1:length(ferm_files)){
  ferm_files[[i]]$Alkali_mL_L <- rescale(
    ferm_files[[i]]$Alkali,
    from = c(
      0,
      last(na.omit(ferm_files[[i]]$Alkali))
    ),
    to = c(
      0,
      alkali_use[i]
    )
  )/0.6
  matched_inok_time <- which.min(abs(ferm_files[[i]]$Time-inok_time[i])) # Finds index of closest time in ferm file to inok_time
  Alkali_neutralization <- ferm_files[[i]]$Alkali_mL_L[matched_inok_time] # Finds alkali value closest to inok time, by which time pH should be neutralized in the reactor  
  ferm_files[[i]]$Alkali_mL_L <- ferm_files[[i]]$Alkali_mL_L - Alkali_neutralization # Substracts the Alkali used for neutralization
  ferm_files[[i]]$Alkali_mL_L[c(0:matched_inok_time)] <- 0 #sets all values before inok to 0
  ferm_files[[i]]$Alkali_mol_L <- ferm_files[[i]]$Alkali_mL_L/1000*alkali_conc[i]
  rm(matched_inok_time, Alkali_neutralization) #Removes temporary variables used in this loop
}



#Rescale Acid to bottle weight, substract Acid used for neutralization:
for (i in 1:length(ferm_files)){
  ferm_files[[i]]$Acid_mL_L <- rescale(
    ferm_files[[i]]$Acid,
    from = c(
      0,
      last(ferm_files[[i]]$Acid)
    ),
    to = c(
      0,
      acid_use[i]
    )
  )/(reactor_volume[i]/1000)
  matched_inok_time <- which.min(abs(ferm_files[[i]]$Time-inok_time[i])) # Finds closest time in ferm file to inok time
  Acid_neutralization <- ferm_files[[i]]$Acid_mL_L[matched_inok_time] # Finds Acid value closest to inok time, by which time pH should be neutralized in the reactor  
  ferm_files[[i]]$Acid_mL_L <- ferm_files[[i]]$Acid_mL_L - Acid_neutralization # Substracts the Acid used for neutralization
  ferm_files[[i]]$Acid_mL_L[ferm_files[[i]]$Acid_mL_L <0 ] <- 0
  ferm_files[[i]]$Acid_mol_L <- ferm_files[[i]]$Acid_mL_L/1000*acid_conc[i]
  rm(matched_inok_time, Acid_neutralization) #Removes temporary variables used in this loop
}






#Rescale OD, adds as aOD_corr column
for (i in 1:length(ferm_files)){
  matched_aod_time <- which.min(abs(ferm_files[[i]]$Time-time_final_od[i])) 
  matched_aOD <- ferm_files[[i]]$aOD[matched_aod_time]
  ferm_files[[i]]$aOD_corr <- NA
  try(ferm_files[[i]]$aOD_corr <- rescale(
    ferm_files[[i]]$aOD,
    
    
    from = c(
      first(ferm_files[[i]]$aOD),
      matched_aOD
    ),
    to = c(
      initial_od[i],
      final_od[i]
    )
  ))
  rm(matched_aOD, matched_aod_time)
}





#Spline Alkali_mL_L and calculate and rAlk_mL_h_L
for (i in 1:length(ferm_files)){
  ferm_files[[i]]$Alkali_ml_L_spl <- NA
  try(ferm_files[[i]]$Alkali_ml_L_spl <- smooth.spline(
    ferm_files[[i]]$Time,
    ferm_files[[i]]$Alkali_mL_L
  )$y)
  ferm_files[[i]]$Alkali_mol_L_spl <- ferm_files[[i]]$Alkali_ml_L_spl/1000*alkali_conc[i]
  ferm_files[[i]]$rAlk_mL_h_L <- NA
  try(ferm_files[[i]] <- ferm_files[[i]] %>%
        mutate(rAlk_mL_h_L = c(0,diff(Alkali_ml_L_spl)/diff(Time))))
}

#Spline Acid_mL_L and calculate and rAci_mL_h_L
for (i in 1:length(ferm_files)){
  ferm_files[[i]]$Acid_ml_L_spl <- NA
  try(ferm_files[[i]]$Acid_ml_L_spl <- smooth.spline(
    ferm_files[[i]]$Time,
    ferm_files[[i]]$Acid_mL_L
  )$y)
  ferm_files[[i]]$Acid_mol_L_spl <- ferm_files[[i]]$Acid_ml_L_spl/1000*alkali_conc[i]
  ferm_files[[i]]$rAci_mL_h_L <- NA
  try(ferm_files[[i]] <- ferm_files[[i]] %>%
        mutate(rAci_mL_h_L = c(0,diff(Acid_ml_L_spl)/diff(Time))))
}

for (i in raw_files){
  aOD_splined <- NA
  try(aOD_splined <- smooth.spline(
    ferm_files[[i]]$Time,
    ferm_files[[i]]$aOD_corr,
    spar = 0.5
  )$y)
  time <- ferm_files[[i]]$Time
  myy <- c()
  try(for (j in c(1:(length(ferm_files[[i]]$aOD_corr)-1))){
    myy <- c(myy,
             summary(lm(log(aOD_splined[c(j:(j+5))])~time[c(j:(j+5))]))$coefficients[2,1]
    )
  })
  myy <- c(myy, NA)
  ferm_files[[i]]$aOD_splined <- aOD_splined
  ferm_files[[i]]$myy_splined_OD <- NA
  ferm_files[[i]]$myy_splined_OD <- myy
  rm(aOD_splined,myy, time)
}



output_data_sheets <- ferm_files 

names(output_data_sheets) <- controller_ids

write_xlsx(output_data_sheets, out_filename)

all_exp <- do.call(rbind, ferm_files) #create dataframe for ggplot, binds ferm_files into one dataframe, columns have to be the same, which the should be if this script worked correctly

alkali_plot <- ggplot(
  data = all_exp,
  aes(
    x = Time,
    y = Alkali_ml_L_spl,
    color = Experiment
  )
) + geom_line() + theme_classic()

rAlk_plot <- ggplot(
  data = all_exp,
  aes(
    x = Time,
    y = rAlk_mL_h_L,
    color = Experiment
  )
) + geom_line() + theme_classic()

aOD_plot <- ggplot(
  data = all_exp,
  aes(
    x = Time,
    y = aOD_corr,
    color = Experiment
  )
) + geom_line() + theme_classic()

myy_plot <- ggplot(
  data = all_exp,
  aes(
    x = Time,
    y = myy_splined_OD,
    color = Experiment
  )
) + geom_line() + theme_classic()




four_plots <- ggarrange(
  alkali_plot,
  rAlk_plot,
  aOD_plot,
  myy_plot,
  ncol = 2,
  nrow = 2
)

X11()
four_plots



svg_name = paste(
  project_id,
  exp_nr,
  "R_plot.svg",
  sep = "_"
)


svg (svg_name, height = 8, width = 10)
four_plots
dev.off()




