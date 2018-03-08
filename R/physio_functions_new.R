# Initialize folders ---------------------------
hrv_initilize <- function() {
  base <- here::here()
  
  if ( !(dir.exists(file.path(base,'data'))) ) {
    dir.create(file.path(base,'data'))
    message('Creating data folder...')
  }
  
  dirs <- list.dirs( file.path(base,'data'), full.names = FALSE )
  need <- c('processed', 'physio_files', 'plots', 'outliers_removed',
            'final', 'hrv')
  
  for (folder in need){
    if (!(folder %in% dirs)) {
      dir.create(file.path(base,'data',folder))
      message('Creating ',folder, '...')
    }
  }
}

# Compress physio data ---------------------------

compress_ecg <- function(file, delete = TRUE) {
  R.utils::gzip(file, remove = delete)
  f <- basename(file)
  message(f, ' succesffuly compressed')
}

# Extract ecg and trigger ---------------------------
extract_ecg <- function(file) {
  basename(file) %>% stringr::str_replace(., '.txt.gz','') -> filename

  data <- data.table::fread(paste0("gunzip -c ", file), header = FALSE)
  ecg <- data$V1
  trigger <- data$V15
  
  if (!(dir.exists('./data/processed'))) {
    dir.create('./data/processed')
  }
  
  # fwrite did not support connections when I wrote this pipeline
  # and it still does not!
  ecg_filename <- file.path('./data/processed',paste0(filename,'_ecg.csv'))
  data.table::fwrite(as.data.frame(ecg), ecg_filename, col.names = FALSE)
  R.utils::gzip(ecg_filename, destname = paste0(ecg_filename,'.gz'))
  
  trigger_filename <- file.path('./data/processed',paste0(filename,'_trigger.csv'))
  data.table::fwrite(as.data.frame(trigger), trigger_filename, col.names = FALSE)
  R.utils::gzip(trigger_filename, destname = paste0(trigger_filename,'.gz'))
  
  message('Raw ECG data and trigger for ', filename, ' written to disk')
}

# Setup ecg data file from raw data ---------------------------
phys_file <- function(file){
  basename(file) %>% stringr::str_replace(., '_ecg.csv.gz','') -> filename
  
  data <- read.csv(file, header = FALSE)
  data$time <- seq(0, ((length(data$V1)-1)*.001), by = .001)
  data <- data[,c(2,1)]
  colnames(data) <- c("time","ecg")
  
  if (!(dir.exists(here('data/physio_files/')))) {
    dir.create(here('data/physio_files/'))
  }
  
  save <- here("data/physio_files", paste0(filename,".phys"))
  
  write.table(data, file = paste0(save,'.csv'), row.names = FALSE, sep =',')
  
  R.utils::gzip(paste0(save,'.csv'), destname = paste0(save,'.gz'))
  
  message('ECG data for ', filename, ' written to disk')
}

# Setup info file ---------------------------
phys_info <- function(file, fs, origin = NA){
  basename(file) %>% stringr::str_replace(., '_ecg.csv.gz','') -> filename
  
  origin <- origin
  fs <- as.numeric(fs)
  data <- cbind.data.frame(origin,fs)
  
  if (!(dir.exists(here('data/physio_files/')))) {
    dir.create(here('data/physio_files/'))
  }
  
  write.table(data, file = here('data/physio_files/', paste0(filename,".info.txt")),
              sep = ",", row.names = FALSE)
}

# Make phys events file ---------------------------
make_events <- function(file){
  basename(file) %>% stringr::str_replace(., '_trigger.csv.gz','') -> filename
  
  events <- read.table(file, header = FALSE)
  
  events$V2 <- seq(1, length(events$V1))
  colnames(events)[2] <- "Time(ms)"
  
  # Set value variable
  colnames(events)[1] <- "Value"
  
  # Make InitTime variable
  op <- options(digits.secs=3)
  events$InitTime <- events$`Time(ms)`/1000
  
  # Make type variable
  values <- c(1.30,1.36)
  condition <- c("Baseline","ER_Task")
  duration <- c(300,600)
  
  match <- cbind.data.frame(values, condition, duration)
  
  events$Type1 <- match(events$Value,match$values)
  events$Type <- match[events$Type1,2]
  events$Duration <- match[events$Type1,3]
  
  # Reorder and save
  events <- events[c(3,5,6,1,2,4)]
  events <- events[,c(-5,-6)]
  events <- events[complete.cases(events$Type),]
  events <- events[!rev(duplicated(rev(events$Type))),]
  
  InitTime <- events$InitTime
  Type <- events$Type
  Duration <- events$Duration
  Value <- events$Value
  
  episodes <- list(InitTime=InitTime,Type=Type,Duration=Duration,Value=Value)
  
  if (!(dir.exists(here('data/processed/')))) {
    dir.create(here('data/processed/'))
  }
  
  save(episodes, file = here('data/processed',paste0(filename,".RData")))
}

# Read and write IBI files ---------------------------
write_ibi <- function(file){
  basename(file) %>% stringr::str_replace(., '.ibi.gz','') -> filename
  
  data <- read.csv(file)
  data <- data[-1]
  
  if (!(dir.exists(here('data/processed/')))) {
    dir.create(here('data/processed/'))
  }
  
  write.table(data, here("data/processed", paste0(filename,".ibi.txt")), sep="\t",
              col.names=FALSE, row.names = FALSE)
}

# Visually inspect IBI series for outliers ---------------------------
outliers_ibi <- function(file, variable= "V1"){
  basename(file) %>% stringr::str_replace(., '.ibi.txt','') -> filename
  
  df <- read.delim(file, header = FALSE)
  df$x <- seq(1, length(df[,variable]))
  
  y <- seq(2, length(df[,variable]))
  d <- diff(df$V1, lag = 1, differences = 1)
  diff <- cbind.data.frame(y,d)
  data <- merge(df,diff, by.x = "x", by.y = "y", all.x = TRUE)
  
  # Flag potential outliers
  data$V2 <- ifelse(data$d >= 300, NA,
                 ifelse(data$d <= -300, NA, data$V1))
  
  if (sum(is.na(data$V2)) >= 1) {
    # Predict values
    data$predict <- forecast::tsclean(data$V2)
    data$V3 <- ifelse(is.na(data$V2)==TRUE, data$predict, data$V2)
  }

  # Plot and compare
  data$colour <- ifelse(is.na(data$V2)==TRUE, "red", "black")
  
  # Plot 1
  plot1 <- ggplot2::ggplot(data) +
    ggplot2::geom_point(ggplot2::aes(y = V1, x = x, colour = colour)) +
    ggplot2::scale_colour_manual(values=c("black","red"), guide = FALSE) +
    ggplot2::ggtitle("Raw Data")
  
  # Plot 2
  plot2 <- ggplot2::ggplot(data) +
    ggplot2::geom_point(ggplot2::aes(y = V3, x = x, colour = colour)) +
    ggplot2::scale_colour_manual(values=c("black","red"), guide = FALSE) +
    ggplot2::ggtitle("Outliers Imputed")
  
  # multiplot and Save
  if (!dir.exists(here('data/plots'))) {
    dir.create(here('data/plots'))
  }
  
  png(filename = here('data/plots', paste(filename,"_plot.png",sep="")), width=1000, height=669,
      units="px")
  gridExtra::grid.arrange(plot1,plot2, ncol=2, top=filename)
  dev.off()
  
  if (!dir.exists(here('data/outliers_removed'))) {
    dir.create(here('data/outliers_removed'))
  }
  
  # Write data and save
  write.table(data$V3, file=here('data/outliers_removed', paste0(filename,"_ecg_clean.txt")), 
              sep="\t", row.names=FALSE, col.names=FALSE)
  
}
