source(here::here('packages.R'))
source(here::here('functions.R'))

#Data from the database#####
con <- dbConnect(PostgreSQL(), dbname = "fishecu_complex-survey_db", user = "fishecuuser", host = "172.21.3.20", password = "f1sh3cuus3r!")

#Latin names
# specs <- data.table(dbGetQuery(conn = con, statement =  "SELECT * FROM fishecu.species;"))#to see the whole table (reservoir). Important to see the reservoir ID number.
# specs <- specs[,c(1,2,10)]
# write.xlsx(specs, here::here('specs.xlsx'))
specs <- setDT(readxl::read_xlsx(here::here('specs.xlsx')))

#Localities
# all_locality <- data.table(dbGetQuery(conn = con, statement = "SELECT *
# FROM fishecu.locality;"))
# all_locality <- all_locality[,c(1,3)]
# write.xlsx(all_locality, here::here('all_locality.xlsx'))
all_locality <- setDT(readxl::read_xlsx(here::here('all_locality.xlsx')))

#Sampling
# all_sampling <- data.table(dbGetQuery(conn = con, statement = "SELECT *
# FROM fishecu.sampling;"))
# all_sampling <- all_sampling[grep("LIP", all_sampling$sa_samplingid), ]
# write.xlsx(all_sampling, here::here('all_sampling.xlsx'))
all_sampling <- setDT(readxl::read_xlsx(here::here('all_sampling.xlsx')))
all_sampling <- merge(all_sampling, all_locality, by = "lo_localityid")
all_sampling[, year := year(sa_date_start)]

# Gillnet deployments####
# all_gill_sto <- data.table(dbGetQuery(conn = con, statement = "SELECT *
# FROM fishecu.gillnet_sampling;"))
# all_gill_sto <- all_gill_sto[grep("LIP", all_gill_sto$sa_samplingid), ]
# write.xlsx(all_gill_sto, here::here('all_gill_sto.xlsx'))
all_gill_sto <- setDT(readxl::read_xlsx(here::here('all_gill_sto.xlsx')))
all_gill_sto <- merge(all_gill_sto, all_sampling[, .(sa_samplingid, sa_date_start, lo_nameczech)], by = "sa_samplingid")

#Catches from gillnet####
# catches_gillnet <- data.table(dbGetQuery(conn = con, statement = paste("SELECT * FROM fishecu.catch
#                                                                   WHERE  sa_samplingid IN ('",paste(all_gill_sto$sa_samplingid, collapse = "','"), "')
#                                                                   ;", sep = "")))
# write.xlsx(catches_gillnet, here::here('catches_gillnet.xlsx'))
catches_gillnet <- setDT(readxl::read_xlsx(here::here('catches_gillnet.xlsx')))


#beachsein####
# all_seining <- data.table(dbGetQuery(conn = con, statement = "SELECT *
# FROM fishecu.beachsein_sampling;"))
# all_seining <- all_seining[grep("LIP", all_seining$sa_samplingid), ]
write.xlsx(all_seining, here::here('all_seining.xlsx'))
all_seining <- setDT(readxl::read_xlsx(here::here('all_seining.xlsx')))
all_seining_lo <- merge(all_seining, all_sampling[, .(sa_samplingid, lo_nameczech, sa_date_start, year)], by = "sa_samplingid")

# Catches from beachsein
# catches_seining <- data.table(dbGetQuery(conn = con, statement = paste("SELECT * FROM fishecu.catch
#                                                                   WHERE  sa_samplingid IN ('",paste(all_seining$sa_samplingid, collapse = "','"), "')
#                                                                   ;", sep = "")))
# write.xlsx(catches_seining, here::here('catches_seining.xlsx'))
catches_seining <- setDT(readxl::read_xlsx(here::here('catches_seining.xlsx')))



#trawling####
# all_trawling <- data.table(dbGetQuery(conn = con, statement = paste("SELECT * FROM fishecu.trawl_sampling;")))
# all_trawling <- all_trawling[grep("LIP", all_trawling$sa_samplingid), ]
# write.xlsx(all_trawling, here::here('all_trawling.xlsx'))
all_trawling <- setDT(readxl::read_xlsx(here::here('all_trawling.xlsx')))
all_trawling <- merge(all_trawling, all_sampling[, .(sa_samplingid, lo_nameczech, sa_date_start, year)], by = "sa_samplingid")

#Catches from trawling
# catches_trawling <- data.table(dbGetQuery(conn = con, statement = paste("SELECT * FROM fishecu.catch
#                                                                   WHERE  sa_samplingid IN ('",paste(all_trawling$sa_samplingid, collapse = "','"), "')
#                                                                   ;", sep = "")))
# write.xlsx(catches_trawling, here::here('catches_trawling.xlsx'))
catches_trawling <- setDT(readxl::read_xlsx(here::here('catches_trawling.xlsx')))

#pas####
# all_pas <- data.table(dbGetQuery(conn = con, statement = paste("SELECT * FROM fishecu.pas_sampling;")))
# all_pas <- all_pas[grep("LIP", all_pas$sa_samplingid), ]
# write.xlsx(all_pas, here::here('all_pas.xlsx'))
all_pas <- setDT(readxl::read_xlsx(here::here('all_pas.xlsx')))

#Catches from pas
# catches_pas <- data.table(dbGetQuery(conn = con, statement = paste("SELECT * FROM fishecu.catch
#                                                                   WHERE  sa_samplingid IN ('",paste(all_pas$sa_samplingid, collapse = "','"), "')
#                                                                   ;", sep = "")))
# write.xlsx(catches_pas, here::here('catches_pas.xlsx'))
catches_pas <- setDT(readxl::read_xlsx(here::here('catches_pas.xlsx')))

#All catch####
# all_catch <- data.table(dbGetQuery(conn = con, statement = paste("SELECT * FROM fishecu.catch;")))
# write.xlsx(all_catch, here::here('all_catch.xlsx'))
# all_catch_lip <- all_catch[grep("LIP", all_catch$sa_samplingid), ]
# write.xlsx(all_catch_lip, here::here('all_catch_lip.xlsx'))
all_catch_lip <- setDT(readxl::read_xlsx(here::here('all_catch_lip.xlsx')))
all_catch_lip <- merge(all_catch_lip, all_sampling[, .(sa_samplingid, year)], by = "sa_samplingid")
# dbDisconnect(con)

#data stomach####
Stomach_content_all <- setDT(readxl::read_xlsx(here::here('Stomach.content.all.xlsx'), sheet = "stomach"))
#Replace NA to 0
Stomach_content_all[, 13:69][is.na(Stomach_content_all[, 13:69])] <- 0
#standardization of sex
Stomach_content_all[["Sex"]][is.na(Stomach_content_all[["Sex"]])] <- "X"
Stomach_content_all$Sex[Stomach_content_all$Sex == "Male"] <- "M"
Stomach_content_all$Sex[Stomach_content_all$Sex == "Female"] <- "F"
Stomach_content_all$Sex[Stomach_content_all$Sex == "?"] <- "X"
#rename for easy coding
colnames(Stomach_content_all)[which(colnames(Stomach_content_all) %in% c("SL (mm)","Weight (g)", "Stomach content") )] <- c("SL","Wg", "Stomach_content")
Stomach_content_all <- Stomach_content_all[, c("Date","Zooplancton", "Chydoridae (Alona)"):=NULL]
#add age group
Stomach_content_all <- merge(Stomach_content_all, all_catch_lip[,.(ct_catchid,ct_agegroup)], by = "ct_catchid", all.x = T)
#add unique id
n_na <- sum(is.na(Stomach_content_all$ct_catchid))
Stomach_content_all$ct_catchid[is.na(Stomach_content_all$ct_catchid)] <- paste0("Fish_", seq_len(n_na))
#Add invertebrate category
Stomach_content_all <- Stomach_content_all %>% 
  mutate(invertebrate = case_when(Chironomidae != 0 ~ 1, #condition 1
                                  Trichoptera != 0 ~ 1,
                                  Koretra != 0 ~ 1,
                                  Jepice != 0 ~ 1,
                                  terr.hmyz != 0 ~ 1,
                                  crayfish != 0 ~ 1,
                                  TRUE ~ 0)) #all other cases
#Add zooplankton category
Stomach_content_all <- Stomach_content_all %>% 
  mutate(zooplankton = case_when(Daphnia != 0 ~ 1, #condition 1
                                 Leptodora != 0 ~ 1,
                                 Diaphanosoma != 0 ~ 1,
                                 Cyclopidae != 0 ~ 1,
                                 TRUE ~ 0)) #all other cases

#separating data####
#candat
Stomach_content_candat <- Stomach_content_all[Species == "candat"]
# size_lip_yoy <- all_catch_lip[sp_speciesid == "candat" & ct_sl > 1,.(Max = max(ct_sl),
#                                  Min = min(ct_sl)),
#                               by =.(year, ct_agegroup)]

Stomach_content_candat <- Stomach_content_candat %>%
  mutate(size_class = case_when(SL <= 160 ~ "YOY",
                                SL %in% 160:300 ~ "older<300",
                                SL > 300 ~ "older>300"))
Stomach_content_candat$size_class <- factor(Stomach_content_candat$size_class, levels = c("YOY", "older<300", "older>300"))
# Stomach_content_candat$ct_agegroup <- dplyr::coalesce(Stomach_content_candat$ct_agegroup, Stomach_content_candat$size_class)

#stomach content analises
Stomach_content_candat_melt <- setDT(melt(Stomach_content_candat[, .(ct_catchid,Year, Species, Fish, invertebrate, zooplankton, size_class)], 
                                       id.vars = c("ct_catchid", "Year", "Species", "size_class"), variable.name = "prey", value.name = "prey_n"))
Stomach_content_candat_melt <- Stomach_content_candat_melt[!prey_n==0,]

ggplot(Stomach_content_candat_melt, aes(x = as.factor(Year), y = prey_n, fill = prey)) +
  geom_col(position = "stack") +
  facet_wrap(~size_class, scales = "free") + 
  labs(x="Year", y="prey n")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#fish in the stomach
Stomach_fish_candat <- Stomach_content_candat[Fish==1, ]
Stomach_fish_candat <- Stomach_fish_candat %>% 
  mutate(cannibal = case_when(candat != 0 ~ 1, #condition 1
                              TRUE ~ 0)) #all other cases

Stomach_fish_candat <- merge(Stomach_fish_candat, fry_abundance[prey_sp == "candat"], by = "Year", all.x = T)
set.seed(3333)
summary(glm(data = Stomach_fish_candat, formula = cannibal ~ size_class+Year+Fry_abundance, family = binomial(link = "logit")))

Stomach_fish_candat2 <- Stomach_fish_candat
Stomach_fish_candat2$cannibal[Stomach_fish_candat2$cannibal == 0] <- "no"
Stomach_fish_candat2$cannibal[Stomach_fish_candat2$cannibal == 1] <- "yes"
fry_abundance <- setDT(readxl::read_xlsx(here::here('Stomach.content.all.xlsx'), sheet = "fry"))
Stomach_fish_candat2 <- merge(Stomach_fish_candat2, fry_abundance[prey_sp == "candat"], by = "Year", all.x = T)


g.yoy <- ggplot(Stomach_fish_candat2[size_class == "YOY"], aes(x = as.factor(Year), y = Fish, fill = forcats::fct_rev(cannibal))) +
  geom_col(position = "stack") +
  geom_line(data=Stomach_fish_candat2[size_class == "YOY"], group=1, mapping=aes(x=as.factor(Year),y=Fry_abundance/20), linewidth = 3)+
  scale_y_continuous(name = "Fish n", 
                     sec.axis = sec_axis(~.*20, name = "Fry n per hectare"))+
  labs(x="Year", fill = "cannibalism", linewidth = NULL, title="YOY")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y.right = element_text(vjust=2),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom")

 g.old <- ggplot(Stomach_fish_candat2[size_class == "older<300"], aes(x = as.factor(Year), y = Fish, fill = forcats::fct_rev(cannibal))) +
  geom_col(position = "stack") +
  geom_line(data=Stomach_fish_candat2[size_class == "older<300"], group=1, mapping=aes(x=as.factor(Year),y=Fry_abundance/50), linewidth = 3)+
  scale_y_continuous(name = "Fish n", 
                     sec.axis = sec_axis(~.*50, name = "Fry n per hectare"))+
  labs(x="Year", fill = "cannibalism", linewidth = NULL, title="older<300")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y.right = element_text(vjust=2),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom")

 g.older <- ggplot(Stomach_fish_candat2[size_class == "older>300"], aes(x = as.factor(Year), y = Fish, fill = forcats::fct_rev(cannibal))) +
  geom_col(position = "stack") +
  geom_line(data=Stomach_fish_candat2[size_class == "older>300"], group=1, mapping=aes(x=as.factor(Year),y=Fry_abundance/20), linewidth = 3)+
  scale_y_continuous(name = "Fish n", 
                     sec.axis = sec_axis(~.*20, name = "Fry n per hectare"))+
  labs(x="Year", fill = "cannibalism", linewidth = NULL, title="older>300")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y.right = element_text(vjust=2),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom")

ggpubr::ggarrange(g.yoy, g.old, g.older,
                  ncol = 2, nrow = 2,
                  common.legend = TRUE, legend = "bottom")

Stomach_fish_candat_melt <- setDT(melt(Stomach_fish_candat[, .(ct_catchid, SL, Year, Species, size_class, cannibal, candat, ouklej, 
                                                               okoun, plotice, jezdik, cejn, cejnek, kaprovitka, okounovitá, Unknown)], 
                                          id.vars = c("ct_catchid", "SL", "Year", "Species", "size_class", "cannibal"), 
                                       variable.name = "prey", value.name = "prey_n"))
Stomach_fish_candat_melt <- Stomach_fish_candat_melt[!prey_n==0,]
summary(glm(data = Stomach_fish_candat_melt, formula = cannibal ~ size_class+Year+prey_n, family = binomial(link = "logit")))

Stomach_fish_candat_melt2 <- Stomach_fish_candat_melt
Stomach_fish_candat_melt2$cannibal[Stomach_fish_candat_melt2$cannibal == 0] <- "no"
Stomach_fish_candat_melt2$cannibal[Stomach_fish_candat_melt2$cannibal == 1] <- "yes"

ggplot(Stomach_fish_candat_melt2, aes(x = as.factor(Year), y = prey_n)) +
  geom_col(position = "stack") +
  facet_wrap(~size_class, scales = "free") + 
  labs(x="Year", y="fish n", fill="Cannibalism")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

Stomach_fish_candat_melt2 <- merge(Stomach_fish_candat_melt2, fry_abundance[prey_sp == "candat"], by = "Year", all.x = T)

#fish size in the stomach
Stomach_fish_candat_size <- setDT(melt(Stomach_fish_candat[, .(ct_catchid, SL, Year, Species, size_class, candat_1,candat_2,candat_3,candat_4, candat_5,
                                                              candat_6, candat_7, candat_8, candat_9, ouklej_1, ouklej_2, ouklej_3, okoun_1,
                                                              okoun_2, okoun_3, okoun_4, plotice_1, plotice_2, jezdik_1, jezdik_2, jezdik_3,
                                                              cejn_1, cejnek_1, cejnek_2, kaprovitka_1, kaprovitka_2,Unknown_1, Unknown_2,
                                                              Unknown_3, Unknown_4, cannibal)], 
                                     id.vars = c("ct_catchid", "Year", "Species", "SL", "size_class", "cannibal"), variable.name = "prey_sp", value.name = "prey_size"))
Stomach_fish_candat_size <- Stomach_fish_candat_size[!prey_size == 0]
Stomach_fish_candat_size[, ':='(ratio_prey = prey_size/SL)]
Stomach_fish_candat_size$prey_sp <- sub("\\_.*", "", as.character(Stomach_fish_candat_size$prey_sp))
summary(glm(data = Stomach_fish_candat_size, formula = cannibal ~ size_class+Year+ratio_prey, family = binomial(link = "logit")))

ggplot(Stomach_fish_candat_size, aes(as.factor(Year), ratio_prey)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2)+
  facet_wrap(~size_class) + 
  labs(x="Year", y="Size ratio")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(Stomach_fish_candat_size, aes(as.factor(Year), ratio_prey)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, aes(color = prey_sp, shape = prey_sp), size = 3)+
  labs(x="Year", y="Size ratio", color="Prey sp", shape="Prey sp")+
  facet_wrap(~size_class, scales = "free_x") + 
  scale_shape_manual(values = rep(15:18, 2)) +
  scale_color_manual(values=c(rep("blue3",1), rep("black",3), rep( "green4", 1), rep( "red", 3))) +
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#okoun
Stomach_content_okoun <- Stomach_content_all[Species == "okoun"]
size_lip_ok <- all_catch_lip[sp_speciesid == "okoun" & ct_sl > 1,.(Max = max(ct_sl),
                                                                     Min = min(ct_sl)),
                              by =.(year, ct_agegroup)]

Stomach_content_okoun <- Stomach_content_okoun %>%
  mutate(size_class = case_when(SL <= 69 ~ "YOY",
                                SL %in% 69:200 ~ "older<200",
                                SL > 200 ~ "older>200"))
Stomach_content_okoun$size_class <- factor(Stomach_content_okoun$size_class, levels = c("YOY", "older<200", "older>200"))
# Stomach_content_candat$ct_agegroup <- dplyr::coalesce(Stomach_content_candat$ct_agegroup, Stomach_content_candat$size_class)

Stomach_content_okoun_melt <- setDT(melt(Stomach_content_okoun[, .(ct_catchid,Year, Species, Fish, invertebrate, zooplankton, size_class)], 
                                          id.vars = c("ct_catchid", "Year", "Species", "size_class"), variable.name = "prey", value.name = "prey_n"))
Stomach_content_okoun_melt <- Stomach_content_okoun_melt[!prey_n==0,]

ggplot(Stomach_content_okoun_melt, aes(x = as.factor(Year), y = prey_n, fill = prey)) +
  geom_col(position = "stack") +
  facet_wrap(~size_class, scales = "free") + 
  labs(x="Year", y="prey n")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#fish in the stomach
Stomach_fish_okoun <- Stomach_content_okoun[Fish==1, ]
Stomach_fish_okoun <- Stomach_fish_okoun %>%
  mutate(cannibal = case_when(okoun != 0 ~ 1, #condition 1
                              TRUE ~ 0)) #all other cases

Stomach_fish_okoun2 <- Stomach_fish_okoun
Stomach_fish_okoun2$cannibal[Stomach_fish_okoun2$cannibal == 0] <- "no"
Stomach_fish_okoun2$cannibal[Stomach_fish_okoun2$cannibal == 1] <- "yes"
# fry_abundance <- setDT(readxl::read_xlsx(here::here('Stomach.content.all.xlsx'), sheet = "fry"))
Stomach_fish_okoun2 <- merge(Stomach_fish_okoun2, fry_abundance[prey_sp == "okoun"], by = "Year", all.x = T)


c.old <- ggplot(Stomach_fish_okoun2[size_class == "older<200"], aes(x = as.factor(Year), y = Fish, fill = forcats::fct_rev(cannibal))) +
  geom_col(position = "stack") +
  geom_line(data=Stomach_fish_okoun2[size_class == "older<200"], group=1, mapping=aes(x=as.factor(Year),y=Fry_abundance/180), linewidth = 3)+
  scale_y_continuous(name = "Fish n",
                     sec.axis = sec_axis(~.*180, name = "Fry n per hectare"))+
  labs(x="Year", fill = "cannibalism", linewidth = NULL, title="older<200")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y.right = element_text(vjust=2),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom")

c.older <- ggplot(Stomach_fish_okoun2[size_class == "older>200"], aes(x = as.factor(Year), y = Fish, fill = forcats::fct_rev(cannibal))) +
  geom_col(position = "stack") +
  geom_line(data=Stomach_fish_okoun2[size_class == "older>200"], group=1, mapping=aes(x=as.factor(Year),y=Fry_abundance/300), linewidth = 3)+
  scale_y_continuous(name = "Fish n",
                     sec.axis = sec_axis(~.*300, name = "Fry n per hectare"))+
  labs(x="Year", fill = "cannibalism", linewidth = NULL, title="older>200")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y.right = element_text(vjust=2),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom")

ggpubr::ggarrange(c.old, c.older,
                  ncol = 2, nrow = 1,
                  common.legend = TRUE, legend = "bottom")

Stomach_fish_okoun_melt <- setDT(melt(Stomach_fish_okoun[, .(ct_catchid, SL, Year, Species, size_class, cannibal, candat, ouklej,
                                                               okoun, plotice, jezdik, cejn, cejnek, kaprovitka, okounovitá, Unknown)],
                                       id.vars = c("ct_catchid", "SL", "Year", "Species", "size_class", "cannibal"),
                                       variable.name = "prey", value.name = "prey_n"))
Stomach_fish_okoun_melt <- Stomach_fish_okoun_melt[!prey_n==0,]
Stomach_fish_okoun_melt2 <- Stomach_fish_okoun_melt
Stomach_fish_okoun_melt2$cannibal[Stomach_fish_okoun_melt2$cannibal == 0] <- "no"
Stomach_fish_okoun_melt2$cannibal[Stomach_fish_okoun_melt2$cannibal == 1] <- "yes"

ggplot(Stomach_fish_okoun_melt, aes(x = as.factor(Year), y = prey_n, fill = prey)) +
  geom_col(position = "stack") +
  facet_wrap(~size_class, scales = "free") +
  labs(x="Year", y="fish n", fill="Prey sp")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#fish size in the stomach
Stomach_fish_okoun_size <- setDT(melt(Stomach_fish_okoun[, .(ct_catchid, SL, Year, Species, size_class, candat_1,candat_2,candat_3,candat_4, candat_5,
                                                               candat_6, candat_7, candat_8, candat_9, ouklej_1, ouklej_2, ouklej_3, okoun_1,
                                                               okoun_2, okoun_3, okoun_4, plotice_1, plotice_2, jezdik_1, jezdik_2, jezdik_3,
                                                               cejn_1, cejnek_1, cejnek_2, kaprovitka_1, kaprovitka_2,Unknown_1, Unknown_2,
                                                               Unknown_3, Unknown_4, cannibal)],
                                       id.vars = c("ct_catchid", "Year", "Species", "SL", "size_class", "cannibal"), variable.name = "prey_sp", value.name = "prey_size"))
Stomach_fish_okoun_size <- Stomach_fish_okoun_size[!prey_size == 0]
Stomach_fish_okoun_size[, ':='(ratio_prey = prey_size/SL)]
Stomach_fish_okoun_size$prey_sp <- sub("\\_.*", "", as.character(Stomach_fish_okoun_size$prey_sp))
summary(glm(data = Stomach_fish_okoun_size, formula = cannibal ~ size_class+Year+ratio_prey, family = binomial(link = "logit")))


ggplot(Stomach_fish_okoun_size, aes(as.factor(Year), ratio_prey)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2)+
  facet_wrap(~size_class) +
  labs(x="Year", y="Size ratio")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(Stomach_fish_okoun_size, aes(as.factor(Year), ratio_prey)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, aes(color = prey_sp, shape = prey_sp), size = 3)+
  labs(x="Year", y="Size ratio", color="Prey sp", shape="Prey sp")+
  facet_wrap(~size_class, scales = "free_x") +
  scale_shape_manual(values = 15:18) +
  scale_color_manual(values=c(rep("blue3",1), rep("green4",1), rep( "black", 1), rep( "red", 1))) +
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
