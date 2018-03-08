freq_analysis_batch <- function(file){
  library(RHRV)
  
  basename(file) %>% stringr::str_replace(., '_ecg_clean.txt','') -> filename
  
  hrv.data = CreateHRVData()
  hrv.data = SetVerbose(hrv.data, FALSE)
  hrv.data = LoadBeatRR(hrv.data, RecordName=file, RecordPath=".", scale = .001)
  
  # we add the info about the episodes
  load(here('data/processed',paste0(filename,".RData")))
  InitTime <- episodes$InitTime
  Type <- episodes$Type
  Duration <- episodes$Duration
  Value <- episodes$Value
  
  hrv.data = AddEpisodes(hrv.data, InitTimes = episodes$InitTime, 
                         Tags = episodes$Type,
                         Durations = episodes$Duration,
                         Values = episodes$Value)
  
  hrv.data = BuildNIHR(hrv.data)
  hrv.data = FilterNIHR(hrv.data)
  
  if (!(dir.exists(here('data/plots/')))){
    dir.create(here('data/plots/'))
    plots <- here('data/plots')
  } else{
    plots <- here('data/plots')
  }
  
  # plot all tags
  png(filename = file.path(plots, paste0(filename,"_tagged_plot.png")), width=1000, height=669,
      units="px")
  PlotNIHR(hrv.data, Tags=episodes$Type)
  dev.off()
  
  hrv.data = InterpolateNIHR(hrv.data, freqhr = 4)
  
  # Perform frequency analysis
  # Calculating spectrogram and power per band using wavelet analysis:
  hrv.data = CreateFreqAnalysis(hrv.data)
  hrv.data = CalculatePowerBand(hrv.data, indexFreqAnalysis = 1, type="wavelet",
                                wavelet="d4",bandtolerance=0.1)
  
  # plot powerband for all files
  png(filename = file.path(plots, paste0(filename,"_powerband.png",sep="")), width=1000, height=669,
      units="px")
  PlotPowerBand(hrv.data, normalized = TRUE, hr = TRUE, Tags = "all")
  dev.off()
  
  # Save the data by stimulus type:
    splitting.data = SplitPowerBandByEpisodes(hrv.data,indexFreqAnalysis = 1, Tag = c("Baseline"))
    Baseline <- log(mean(splitting.data$OutEpisodes$HF))
    
    splitting.data = SplitPowerBandByEpisodes(hrv.data,indexFreqAnalysis = 1, Tag = c("ER_Task"))
    Task <- log(mean(splitting.data$OutEpisodes$HF))
    
    subject_nr <- readr::parse_number(file)
    sub <- cbind.data.frame(subject_nr,Baseline,Task)
    
    if (!(dir.exists(here('data/hrv')))){
      dir.create(here('data/hrv'))
      hrv <- here('data/hrv')
    } else {
      hrv <- here('data/hrv')
    }
    
    write.table(sub, file = file.path(hrv, "data_all_physio.csv"), sep = ",", append = TRUE,
                col.names = FALSE, row.names = FALSE)
  }
    
    