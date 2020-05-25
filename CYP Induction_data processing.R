# Load packages 
pacman::p_load(pacman, dplyr, tidyr, stringr, readxl, writexl, outliers)

# Get Rawdata file name
Name_rawdata_file <- list.files(pattern = "*.xls", full.names = T)

# Load raw data
Rawdata <- read_excel(Name_rawdata_file, sheet = "Raw QPCR CT data")

# Cleanup raw data and subset datafrma for individual reporter gene 
Val_replicate <- 3                                                                     # 3 for triplicate
Val_pct_CV_cutoff <- 12                                                                # Cutoff value of %CV to remove outlier
Val_pct_PC_NR <- 2                                                                     # Cutoff value of %PC to set as "NR"
Name_PosCtrls <- c("Ome", "Pheb", "Rifa")                                              # Name of positive control in sample name for 1A2, 2B6 and 3A4
Rawdata <- Rawdata[ , c(4, 2, 5, 7, 9)]                                                # Cleanup to have only essential columns for further data processing.
Rawdata <- drop_na(Rawdata, "Sample Name")                                             # Drop rows containing NA's in Sample Name column
Rawdata$CT <- as.numeric(Rawdata$CT)                                                   # Convert CT from character into numeric
Val_reporter_1 <- Rawdata$Reporter[1]                                                  # Define reporter as variable to use as filter criteria
Val_reporter_2 <- Rawdata$Reporter[2]                                                  # Define reporter as variable to use as filter criteria

Df_reporter_1 <- filter(Rawdata, Rawdata$Reporter == Val_reporter_1)                   # Select all samples with reporter 1
Df_reporter_2 <- filter(Rawdata, Rawdata$Reporter == Val_reporter_2)                   # Select all samples wiht reporter 2
Df_all <- left_join(Df_reporter_1, Df_reporter_2, "Well Position")                     # Join dataframe by well position
Df_all <- mutate(Df_all, "Delta_CT" = Df_all$CT.y - Df_all$CT.x)                       # Add a new column and calculate Delta CT values
Df_all <- Df_all[ , -6]                                                                # Remove duplicate sample name column

Row_num_DMSO <- str_which(Df_all$`Sample Name.x`, "DMSO")                              # Get row number of all DMSO samples
Row_num_PC <- str_which(Df_all$`Sample Name.x`, "PC_")                                 # Get row number of all PosCtrl samples

Df_all$`Sample Name.x`[Row_num_DMSO] <- str_replace(Df_all$`Sample Name.x`[Row_num_DMSO], "_", "_Vehicle_")   
Df_all$`Sample Name.x`[Row_num_PC] <- str_replace(Df_all$`Sample Name.x`[Row_num_PC], "PC_", "")

Df_all <- Df_all %>% separate(`Sample Name.x`, c("ID", "Replicate"), "_0")

# Create function to calculate avg, sd and %CV without the value most differing from mean, and fold of change
mean.rm.o <- function(Data) { zval_1 <- rm.outlier(Data)
                              zval_2 <- mean(zval_1, na.rm = TRUE)
                              return(zval_2)} 

pct_CV <- function(Data) { zval_1 <- mean(Data, na.rm = TRUE)
                           zval_2 <- sd(Data, na.rm = TRUE)
                           zval_3 <- zval_2/zval_1*100
                           zval_4 <- round(zval_3, 1)
                           return(zval_4)} 

pct_CV.rm.o <- function(Data) { zval_1 <- rm.outlier(Data)
                                zval_2 <- mean(zval_1, na.rm = TRUE)
                                zval_3 <- sd(zval_1, na.rm = TRUE)
                                zval_4 <- zval_3/zval_2*100
                                zval_5 <- round(zval_4, 1)
                                return(zval_5)} 

fold.of.change <- function(ddCT) { zval_1 <- 2^(-ddCT)
                                   return(zval_1)}

# Data processing for 1A2
Df_1A2 <- Df_all %>% filter(`Target Name.y` == "CYP1A2") %>%
                     group_by(ID) %>%
                     summarize(dCT_avg = mean(Delta_CT, na.rm = TRUE),
                               Pct_CV = pct_CV(Delta_CT),
                               dCT_avg_rm_o = mean.rm.o(Delta_CT),
                               Pct_CV_rm_o = pct_CV.rm.o(Delta_CT))

Df_1A2 <- Df_1A2 %>% mutate("dCT_avg_1A2" = ifelse(Df_1A2$Pct_CV < Val_pct_CV_cutoff, Df_1A2$dCT_avg, 
                                                 ifelse(Df_1A2$Pct_CV_rm_o < Val_pct_CV_cutoff, Df_1A2$dCT_avg_rm_o, NA))
                            )                                                                            

Val_dCT_DMSO_1A2 <- Df_1A2$`dCT_avg_1A2`[str_which(Df_1A2$ID, "DMSO")]                                  # Avg dCT of DMSO (vehicle control)
Df_1A2 <- Df_1A2 %>% mutate("ddCT_avg_1A2" = dCT_avg_1A2 - Val_dCT_DMSO_1A2)                            # Calculate delta delta CT (avg)
Df_1A2 <- Df_1A2 %>% mutate("Fold_of_change_1A2" = fold.of.change(Df_1A2$ddCT_avg_1A2))                 # Calculate fold of change
Val_fold_change_PC_1A2 <- Df_1A2$`Fold_of_change_1A2`[str_which(Df_1A2$ID, Name_PosCtrls[1])]           # Fold of change of positive control
Df_1A2 <- Df_1A2 %>% mutate("%PC_1A2" = (Fold_of_change_1A2-1)/(Val_fold_change_PC_1A2-1)*100)          # Calculate %PC
Df_1A2_final <- Df_1A2[, c(1,8:9)]                                                                      # Get a dataframe of only ID, fold of change and %PC

# Data processing for 2B6
Df_2B6 <- Df_all %>% filter(`Target Name.y` == "CYP2B6") %>%
  group_by(ID) %>%
  summarize(dCT_avg = mean(Delta_CT, na.rm = TRUE),
            Pct_CV = pct_CV(Delta_CT),
            dCT_avg_rm_o = mean.rm.o(Delta_CT),
            Pct_CV_rm_o = pct_CV.rm.o(Delta_CT))

Df_2B6 <- Df_2B6 %>% mutate("dCT_avg_2B6" = ifelse(Df_2B6$Pct_CV < Val_pct_CV_cutoff, Df_2B6$dCT_avg, 
                                                   ifelse(Df_2B6$Pct_CV_rm_o < Val_pct_CV_cutoff, Df_2B6$dCT_avg_rm_o, NA))
)

Val_dCT_DMSO_2B6 <- Df_2B6$`dCT_avg_2B6`[str_which(Df_2B6$ID, "DMSO")]
Df_2B6 <- Df_2B6 %>% mutate("ddCT_avg_2B6" = dCT_avg_2B6 - Val_dCT_DMSO_2B6)
Df_2B6 <- Df_2B6 %>% mutate("Fold_of_change_2B6" = fold.of.change(Df_2B6$ddCT_avg_2B6))
Val_fold_change_PC_2B6 <- Df_2B6$`Fold_of_change_2B6`[str_which(Df_2B6$ID, Name_PosCtrls[2])]  
Df_2B6 <- Df_2B6 %>% mutate("%PC_2B6" = (Fold_of_change_2B6-1)/(Val_fold_change_PC_2B6-1)*100)
Df_2B6_final <- Df_2B6[, c(1,8:9)]

# Data processing for 3A4
Df_3A4 <- Df_all %>% filter(`Target Name.y` == "CYP3A4") %>%
  group_by(ID) %>%
  summarize(dCT_avg = mean(Delta_CT, na.rm = TRUE),
            Pct_CV = pct_CV(Delta_CT),
            dCT_avg_rm_o = mean.rm.o(Delta_CT),
            Pct_CV_rm_o = pct_CV.rm.o(Delta_CT))

Df_3A4 <- Df_3A4 %>% mutate("dCT_avg_3A4" = ifelse(Df_3A4$Pct_CV < Val_pct_CV_cutoff, Df_3A4$dCT_avg, 
                                                   ifelse(Df_3A4$Pct_CV_rm_o < Val_pct_CV_cutoff, Df_3A4$dCT_avg_rm_o, NA))
                            )

Val_dCT_DMSO_3A4 <- Df_3A4$`dCT_avg_3A4`[str_which(Df_3A4$ID, "DMSO")]
Df_3A4 <- Df_3A4 %>% mutate("ddCT_avg_3A4" = dCT_avg_3A4 - Val_dCT_DMSO_3A4)
Df_3A4 <- Df_3A4 %>% mutate("Fold_of_change_3A4" = fold.of.change(Df_3A4$ddCT_avg_3A4))
Val_fold_change_PC_3A4 <- Df_3A4$`Fold_of_change_3A4`[str_which(Df_3A4$ID, Name_PosCtrls[3])]  
Df_3A4 <- Df_3A4 %>% mutate("%PC_3A4" = (Fold_of_change_3A4-1)/(Val_fold_change_PC_3A4-1)*100)
Df_3A4_final <- Df_3A4[, c(1,8:9)]

# Compile fold of change and %PC for all three target genes
Summary_all <- full_join(Df_1A2_final, Df_2B6_final, "ID")
Summary_all <- full_join(Summary_all, Df_3A4_final, "ID")
Summary_all[ ,2:7] <- round(Summary_all[ ,2:7], 2)
Summary_all <- Summary_all %>% separate(ID, c("ID", "Concentration"), "_")
Summary_all$Concentration <- str_replace(Summary_all$Concentration, "uM", "")

Val_PC_1A2_fold_change <- Summary_all$Fold_of_change_1A2[str_which(Summary_all$ID, Name_PosCtrls[1])]
Val_PC_1A2_Pct_PC <- Summary_all$`%PC_1A2`[str_which(Summary_all$ID, Name_PosCtrls[1])]
Val_PC_2B6_fold_change <- Summary_all$Fold_of_change_2B6[str_which(Summary_all$ID, Name_PosCtrls[2])]
Val_PC_2B6_Pct_PC <- Summary_all$`%PC_2B6`[str_which(Summary_all$ID, Name_PosCtrls[2])]
Val_PC_3A4_fold_change <- Summary_all$Fold_of_change_3A4[str_which(Summary_all$ID, Name_PosCtrls[3])]
Val_PC_3A4_Pct_PC <- Summary_all$`%PC_3A4`[str_which(Summary_all$ID, Name_PosCtrls[3])]
 
List_PosCtrl_final <- c("Positive Control", "PC", Val_PC_1A2_fold_change, Val_PC_1A2_Pct_PC,                  # Compile values of positive control of each target gene
                                                  Val_PC_2B6_fold_change, Val_PC_2B6_Pct_PC, 
                                                  Val_PC_3A4_fold_change, Val_PC_3A4_Pct_PC)

Summary_all <- rbind(Summary_all, List_PosCtrl_final)                                           
Summary_all <- Summary_all[!(Summary_all$ID %in% Name_PosCtrls), ]                                            # Remove rows of individual positive control
Summary_all[ ,3:8] <- lapply(Summary_all[ ,3:8], as.numeric)
Summary_all_2 <- Summary_all                                                                                  # Duplicate summary for export

Summary_all$`%PC_1A2`[Summary_all$`%PC_1A2` < Val_pct_PC_NR] <- "NR"  
Summary_all$`%PC_2B6`[Summary_all$`%PC_2B6` < Val_pct_PC_NR] <- "NR" 
Summary_all$`%PC_3A4`[Summary_all$`%PC_3A4` < Val_pct_PC_NR] <- "NR" 

# Export
Val_current_date <- Sys.Date()                                                                                # Get the current date to attach to the file name

List_Summary <- list("Summary" = Summary_all)
write_xlsx(List_Summary,                                                                                      # Export summary to an excel file
           paste(Val_current_date, " CYP Induction_Summary", ".xlsx", sep = ""))  

List_processed_data <- list("1A2 Processed data" = Df_1A2, 
                            "2B6 Processed data" = Df_2B6, 
                            "3A4 Processed data" = Df_3A4,
                            "Summary_no NR" = Summary_all_2,
                            "Summary_final" = Summary_all)
write_xlsx(List_processed_data,                                                                               # Export processed data to an excel file
           paste(Val_current_date, " CYP Induction_Processed data", ".xlsx", sep = ""))  

##### CLEAN UP ######

rm(list = ls())   # Clear environment
p_unload(all)     # Remove all add-ons
cat("\014")       # ctrl+L # Clear console