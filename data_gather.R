source(here::here('packages.R'))
source(here::here('functions.R'))

#Data from the database#####
# con <- dbConnect(PostgreSQL(), dbname = "fishecu_complex-survey_db", user = "fishecuuser", host = "172.21.3.20", password = "f1sh3cuus3r!")

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
Stomach_content_all <- setDT(readxl::read_xlsx(here::here('Stomach.content.all.xlsx'), sheet = "candat"))
#Replace NA to 0
Stomach_content_all[, 8:59][is.na(Stomach_content_all[, 8:59])] <- 0
#standardization of sex
Stomach_content_all[["Sex"]][is.na(Stomach_content_all[["Sex"]])] <- "X"
Stomach_content_all$Sex[Stomach_content_all$Sex == "Male"] <- "M"
Stomach_content_all$Sex[Stomach_content_all$Sex == "Female"] <- "F"
Stomach_content_all$Sex[Stomach_content_all$Sex == "?"] <- "X"
#add age group
Stomach_content_all <- merge(Stomach_content_all, all_catch_lip[,.(ct_catchid,ct_agegroup)], by = "ct_catchid", all.x = T)
#add unique id
n_na <- sum(is.na(Stomach_content_all$ct_catchid))
Stomach_content_all$ct_catchid[is.na(Stomach_content_all$ct_catchid)] <- paste0("Fish_", seq_len(n_na))

Stomach_content_all <- Stomach_content_all %>%
  mutate(size_class = case_when(SL <= 160 ~ "YOY",
                                SL %in% 160:300 ~ "older<300",
                                SL > 300 ~ "older>300"))
Stomach_content_all$size_class <- factor(Stomach_content_all$size_class, levels = c("YOY", "older<300", "older>300"))
Stomach_content_all <- Stomach_content_all %>%
  mutate(Dataset = case_when(Year <= 1970 ~ "19's",
                            Year > 2009 ~ "20's"))

ggplot(Stomach_content_all, aes(x = SL)) + 
  geom_histogram(stat="count")


#fish in the stomach
Stomach_fish_candat_melt <- setDT(melt(Stomach_content_all[, .(ct_catchid, SL, Year, size_class, candat, ouklej, 
                                                               okoun, plotice, jezdik, cejn, cejnek, kaprovitka, okounovitá, Unknown, Dataset)], 
                                          id.vars = c("ct_catchid", "SL", "Year", "size_class", "Dataset"), 
                                       variable.name = "prey", value.name = "prey_n"))
# Stomach_fish_candat_melt <- Stomach_fish_candat_melt[!prey_n==0,]
# Stomach_fish_candat_melt[, sp_grouped := fct_lump(f = prey, prop = 0.05, w = prey_n)]
Stomach_fish_candat_melt <- merge(Stomach_fish_candat_melt, specs[,.(sp_speciesid,sp_taxonomicorder, sp_scientificname)], by.x = "prey", by.y = "sp_speciesid", all.x = T)
Stomach_fish_candat_melt$sp_taxonomicorder[Stomach_fish_candat_melt$prey == "kaprovitka"] <- "Cypriniformes"
Stomach_fish_candat_melt$sp_taxonomicorder[Stomach_fish_candat_melt$prey == "okounovitá"] <- "Perciformes"
Stomach_fish_candat_melt$sp_taxonomicorder[Stomach_fish_candat_melt$prey == "Unknown"] <- "Unknown"
Stomach_fish_candat_melt$sp_scientificname[Stomach_fish_candat_melt$prey == "kaprovitka"] <- "Cyprinid"
Stomach_fish_candat_melt$sp_scientificname[Stomach_fish_candat_melt$prey == "okounovitá"] <- "Percid"
Stomach_fish_candat_melt$sp_scientificname[Stomach_fish_candat_melt$prey == "Unknown"] <- "Unknown"
Stomach_fish_candat_melt$sp_name[Stomach_fish_candat_melt$prey == "kaprovitka"] <- "Cyprinid"
Stomach_fish_candat_melt$sp_name[Stomach_fish_candat_melt$prey == "Unknown"] <- "Unknown"
Stomach_fish_candat_melt$sp_name[Stomach_fish_candat_melt$prey == "okounovitá"] <- "Percid"
Stomach_fish_candat_melt$sp_name[Stomach_fish_candat_melt$sp_scientificname == "Sander lucioperca"] <- "Pikeperch"
Stomach_fish_candat_melt$sp_name[Stomach_fish_candat_melt$sp_scientificname == "Abramis brama"] <- "Bream"
Stomach_fish_candat_melt$sp_name[Stomach_fish_candat_melt$sp_scientificname == "Gymnocephalus cernua"] <- "Ruffe"
Stomach_fish_candat_melt$sp_name[Stomach_fish_candat_melt$sp_scientificname == "Perca fluviatilis"] <- "Perch"
Stomach_fish_candat_melt$sp_name[Stomach_fish_candat_melt$sp_scientificname == "Alburnus alburnus"] <- "Bleak"
Stomach_fish_candat_melt$sp_name[Stomach_fish_candat_melt$sp_scientificname == "Rutilus rutilus"] <- "Roach"

ggplot(Stomach_fish_candat_melt, aes(x = prey_n)) + 
  geom_histogram(stat="count")

#prey n sp ####
sum_sp_dec <- Stomach_fish_candat_melt[!sp_scientificname == "Unknown",.(prey_n = sum(prey_n),
                                          predator_n = uniqueN(ct_catchid)),
                                         by =.(Dataset, sp_name)]
sum_sp_dec$predato_r <- sum_sp_dec$prey_n/sum_sp_dec$predator_n
ggplot(sum_sp_dec, aes(x = sp_name, y = predato_r, fill = Dataset)) +
  geom_col(position="dodge") +
  labs(x = "Species", y = "Mean prey per predator", fill = "Dataset")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 45, hjust =1.1, vjust = 1.05, face = "italic"),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

sum_dec_y <- Stomach_fish_candat_melt[,.(prey_n = sum(prey_n),
                                            predator_n = uniqueN(ct_catchid), 
                                            predato_r = sum(prey_n)/uniqueN(ct_catchid)),
                                       by =.(Dataset)]

ggplot(sum_dec_y, aes(x = predato_r)) + 
  geom_histogram(stat="count")

skewness(sum_dec_y$predato_r)
shapiro.test(sum_dec_y$predato_r)
hist(sum_dec_y$predato_r, breaks = 100)
set.seed(3333)
mzyf <- glm(predato_r ~ Year, data = sum_dec_y)
mzy <-glm(predato_r ~ Year + Dataset, data = sum_dec_y)
anova(mzy, mzyf)
summary(mzy)
summary(mzyf)
dfun(mzyf)
dfun(mzy)
with(summary(mzy), 1 - deviance/null.deviance)
with(summary(mzyf), 1 - deviance/null.deviance)

sum_sp_dec_y <- Stomach_fish_candat_melt[!sp_scientificname == "Unknown",.(prey_n = sum(prey_n),
                                            predator_n = uniqueN(ct_catchid), 
                                            predato_r = sum(prey_n)/uniqueN(ct_catchid)),
                                         by =.(Dataset, Year, sp_name)]
teste <- sum_sp_dec_y[,.(predato_m = mean(predato_r), 
                         predator_se = plotrix::std.error(predato_r)), 
                      by =.(Dataset, sp_name)]


ggplot(sum_sp_dec_y[!sp_name %in% c("Percid","Cyprinid")], #cpue island/control
       aes(x = sp_name, y = predato_r, fill = Dataset)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.9), outlier.shape = NA) +
  geom_jitter(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.1)) +
  # facet_wrap(~Time, scale = 'free_y')+
  theme(strip.text = element_text(face = "italic")) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 5, col = 'black', position = position_dodge(.9)) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 4, col = 'white', position = position_dodge(.9)) +
  labs(x = "Species", y = 'Mean N prey per Pikeperch')+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 45, hjust =1.1, vjust = 1.05, face = "italic"),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

shapiro.test(sum_sp_dec_y$predato_r)
skewness(sum_sp_dec_y$predato_r)
mzy_qp_sp <- glm(predato_r ~ Dataset + sp_name, data = sum_sp_dec_y, family = "quasipoisson")
summary(mzy_qp_sp)
dfun(mzy_qp_sp)
with(summary(mzy_qp_sp), 1 - deviance/null.deviance)

dcast_sp_dec_y <- dcast(data = sum_sp_dec_y, formula = Dataset + Year ~ sp_name, value.var = "predato_r")

lapply(dcast_sp_dec_y[, c(3:11)],shapiro.test)
lapply(dcast_sp_dec_y[, c(3:11)],skewness)

set.seed(9999)
mzy_zi_sp <- lapply(dcast_sp_dec_y[, c(3:11)], 
                          function(x) glm(formula = x ~ Dataset, 
                                                   data = dcast_sp_dec_y, family = "quasipoisson"))
summary(mzy_zi_sp$Bream)
dfun(mzy_zi_sp$Bream)
summary(mzy_zi_sp$Bleak)
dfun(mzy_zi_sp$Bleak)
summary(mzy_zi_sp$`Silver bream`)
dfun(mzy_zi_sp$`Silver bream`)
summary(mzy_zi_sp$Cyprinid)
dfun(mzy_zi_sp$Cyprinid)
summary(mzy_zi_sp$Percid)
dfun(mzy_zi_sp$Percid)
summary(mzy_zi_sp$Ruffe)
dfun(mzy_zi_sp$Ruffe)
with(summary(mzy_zi_sp$Ruffe), 1 - deviance/null.deviance)
summary(mzy_zi_sp$Perch)
dfun(mzy_zi_sp$Perch)
with(summary(mzy_zi_sp$Perch), 1 - deviance/null.deviance)
summary(mzy_zi_sp$Roach)
dfun(mzy_zi_sp$Roach)
summary(mzy_zi_sp$Pikeperch)
dfun(mzy_zi_sp$Pikeperch)
with(summary(mzy_zi_sp$Pikeperch), 1 - deviance/null.deviance)

#family analises
sum_or_dec_y <- Stomach_fish_candat_melt[!sp_taxonomicorder == "Unknown",.(prey_n = sum(prey_n),
                                            predator_n = uniqueN(ct_catchid), 
                                            predato_r = sum(prey_n)/uniqueN(ct_catchid)),
                                         by =.(Dataset, Year, sp_taxonomicorder)]

ggplot(sum_or_dec_y, #cpue island/control
       aes(x = sp_taxonomicorder, y = predato_r, fill = Dataset)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.9), outlier.shape = NA) +
  geom_jitter(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.1)) +
  scale_x_discrete(labels=c('Cyprinid', 'Percid'))+
  # facet_wrap(~Time, scale = 'free_y')+
  theme(strip.text = element_text(face = "italic")) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 5, col = 'black', position = position_dodge(.9)) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 4, col = 'white', position = position_dodge(.9)) +
  labs(x = "Family", y = 'Mean N prey per Pikeperch')+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 45, hjust =1.1, vjust = 1.05, face = "italic"),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

shapiro.test(sum_or_dec_y$predato_r)
skewness(sum_or_dec_y$predato_r)
mzy_qp_or <- glm(predato_r ~ Dataset + sp_taxonomicorder, data = sum_or_dec_y, family = "quasipoisson")
summary(mzy_qp_or)
dfun(mzy_qp_or)
with(summary(mzy_qp_or), 1 - deviance/null.deviance)

dcast_or_dec_y <- dcast(data = sum_or_dec_y, formula = Dataset + Year ~ sp_taxonomicorder, value.var = "predato_r")
lapply(dcast_or_dec_y[, c(3:4)],shapiro.test)
lapply(dcast_or_dec_y[, c(3:4)],skewness)

set.seed(9999)
mzy_zi_or <- lapply(dcast_or_dec_y[, c(3:4)], 
                    function(x) glm(formula = x ~ Dataset, 
                                    data = dcast_or_dec_y, family = "quasipoisson"))

summary(mzy_zi_or$Cypriniformes)
dfun(mzy_zi_or$Cypriniformes)
with(summary(mzy_zi_or$Cypriniformes), 1 - deviance/null.deviance)
summary(mzy_zi_or$Perciformes)
dfun(mzy_zi_or$Perciformes)
with(summary(mzy_zi_or$Perciformes), 1 - deviance/null.deviance)


#fish size in the stomach
Stomach_fish_candat_size <- setDT(melt(Stomach_content_all[, .(ct_catchid, SL, Year, Species, candat_1,candat_2,candat_3,candat_4, candat_5,
                                                              candat_6, candat_7, candat_8, candat_9, ouklej_1, ouklej_2, ouklej_3, okoun_1,
                                                              okoun_2, okoun_3, okoun_4, okoun_5, okoun_6, okoun_7, okoun_8, okoun_9, okoun_10, okoun_11, 
                                                              plotice_1, plotice_2, jezdik_1, jezdik_2, jezdik_3, jezdik_4, jezdik_5, jezdik_6, jezdik_7,
                                                              cejn_1, cejnek_1, cejnek_2, kaprovitka_1, kaprovitka_2, kaprovitka_3, Unknown_1, Unknown_2,
                                                              Unknown_3, Unknown_4, Dataset)], 
                                     id.vars = c("ct_catchid", "Year", "Species", "SL", "Dataset"), variable.name = "prey_sp", value.name = "prey_size"))
Stomach_fish_candat_size <- Stomach_fish_candat_size[!prey_size == 0]
Stomach_fish_candat_size[, ':='(ratio_prey = prey_size/SL)]
Stomach_fish_candat_size$prey_sp <- sub("\\_.*", "", as.character(Stomach_fish_candat_size$prey_sp))
Stomach_fish_candat_size$prey_sp <- factor(Stomach_fish_candat_size$prey_sp, levels = c("candat", "okoun", "jezdik", "cejn", "cejnek", "kaprovitka", "ouklej", "plotice", "Unknown"))
Stomach_fish_candat_size <- merge(Stomach_fish_candat_size, specs[,.(sp_speciesid,sp_taxonomicorder, sp_scientificname)], by.x = "prey_sp", by.y = "sp_speciesid", all.x = T)
Stomach_fish_candat_size$sp_taxonomicorder[Stomach_fish_candat_size$prey_sp == "kaprovitka"] <- "Cypriniformes"
Stomach_fish_candat_size$sp_scientificname[Stomach_fish_candat_size$prey_sp == "kaprovitka"] <- "Cypriniformes"
Stomach_fish_candat_size$sp_taxonomicorder[Stomach_fish_candat_size$prey_sp == "Unknown"] <- "Unknown"
Stomach_fish_candat_size$sp_scientificname[Stomach_fish_candat_size$prey_sp == "Unknown"] <- "Unknown"

ggplot(Stomach_fish_candat_size, aes(x = prey_size)) + 
  geom_histogram(stat="count")

Stomach_percid_size <- Stomach_fish_candat_size[sp_taxonomicorder == "Perciformes"]
Stomach_cyprinid_size <- Stomach_fish_candat_size[sp_taxonomicorder == "Cypriniformes"]

#table?
# Stomach_fish_candat_size$n_prey <- 1
# sum_size_dec <- Stomach_fish_candat_size[,.(Mean = round(mean(prey_size, na.rm = T), 2),
#                                              SE = round(plotrix::std.error(prey_size), 2),
#                                              Max = max(prey_size),
#                                              Min = min(prey_size),
#                                             N = length(prey_size)),
#                                by =.(Dataset, sp_taxonomicorder)]
# sum_size_dec
# sum_size_dec[is.na(sum_size_dec)] <- 0
# # write.xlsx(sum_size_dec, here::here('sum_size_dec.xlsx'))
# sum_size_d <- Stomach_fish_candat_size[,.(prey_n = sum(n_prey),
#                                           predator_n = length(unique(ct_catchid))),
#                                        by =.(Dataset)]
# sum_size_d$mean_p <- sum_size_d$prey_n/sum_size_d$predator_n

shapiro.test(Stomach_fish_candat_size$ratio_prey)
skewness(Stomach_fish_candat_size$ratio_prey)
ggplot(Stomach_fish_candat_size, aes(x = ratio_prey)) + 
  geom_histogram()
shapiro.test(Stomach_cyprinid_size$ratio_prey)
skewness(Stomach_cyprinid_size$ratio_prey)
ggplot(Stomach_cyprinid_size, aes(x = ratio_prey)) + 
  geom_histogram()
shapiro.test(Stomach_percid_size$ratio_prey)
skewness(Stomach_percid_size$ratio_prey)
ggplot(Stomach_percid_size, aes(x = ratio_prey)) + 
  geom_histogram()

set.seed(3333)
summary(glm(data = Stomach_cyprinid_size, formula = ratio_prey ~ SL * Dataset))
summary(glm(data = Stomach_percid_size, formula = ratio_prey ~ SL * Dataset, family = "quasipoisson"))
summary(stats::aov(data = Stomach_fish_candat_size, formula = ratio_prey ~ SL+Dataset+sp_taxonomicorder))

#cyprinid
set.seed(3333)
m1 <- lme4::lmer(data = Stomach_fish_candat_size, formula = ratio_prey ~ SL + Dataset + (1|sp_taxonomicorder), REML = FALSE)
m2 <- glm(data = Stomach_fish_candat_size, formula = ratio_prey ~ SL + Dataset, family = "quasipoisson")
# with(summary(m1), 1 - deviance/null.deviance)
anova(m1, m2)
car::Anova(m1)
summary(m1)
summary(m2)
#percid
set.seed(3333)
m4 <- lme4::lmer(data = Stomach_fish_candat_size[!sp_taxonomicorder == "Unknown"], formula = ratio_prey ~ SL + Dataset + (1|sp_taxonomicorder), REML = FALSE)
m5 <- lme4::lmer(data = Stomach_fish_candat_size[!sp_taxonomicorder == "Unknown"], formula = ratio_prey ~ SL + (1|sp_taxonomicorder), REML = FALSE)
# with(summary(m1), 1 - deviance/null.deviance)
anova(m4, m5)
car::Anova(m4)
summary(m4)
sjPlot::tab_model(m4)


ggplot(Stomach_fish_candat_size, aes(as.factor(Year), ratio_prey)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2)+
  facet_wrap(~size_class) + 
  labs(x="Year", y="Size ratio")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggplot(Stomach_fish_candat_size, aes(Dataset, ratio_prey)) +
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
  geom_jitter(width = 0.2, aes(color = sp_taxonomicorder, shape = sp_taxonomicorder), size = 3)+
  labs(x="Year", y="Size ratio", color="Prey order", shape="Prey order")+
  # facet_wrap(~size_class, scales = "free") + 
  scale_shape_manual(values = rep(15:18, )) +
  scale_color_manual(values=c(rep("blue3",1), rep("black",1), rep("green4", 1))) +
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(Stomach_fish_candat_size, aes(Dataset, ratio_prey)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, aes(color = sp_taxonomicorder, shape = sp_taxonomicorder), size = 3)+
  labs(x="Year", y="Size ratio", color="Prey order", shape="Prey order")+
  facet_wrap(~size_class, scales = "free_x") + 
  scale_shape_manual(values = rep(15:18, )) +
  scale_color_manual(values=c(rep("blue3",1), rep("black",1), rep("green4", 1))) +
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(Stomach_fish_candat_size[prey_sp %in% c("okoun", "jezdik", "plotice")], aes(SL, ratio_prey)) +
  geom_jitter(width = 0.2, aes(color = Dataset))+
  facet_wrap(~prey_sp, scales = "free_y", ncol = 1) + 
  geom_smooth(method='glm', formula= y~x, method.args = list(family = "quasipoisson"), aes(color = Dataset, fill = Dataset))+
  labs(x="Pikeperch SL (mm)", y="Prey size")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(Stomach_fish_candat_size[!size_class == "YOY"& prey_sp %in% c("okoun", "jezdik", "plotice")], aes(SL, ratio_prey)) +
  geom_jitter(width = 0.2, aes(color = Dataset))+
  facet_wrap(~prey_sp, scales = "free_y", ncol = 1) + 
  geom_smooth(method='lm', formula= y~x, aes(color = Dataset, fill = Dataset))+
  labs(x="SL (mm)", y="Size ratio")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#order
ggplot(Stomach_fish_candat_size[!size_class == "YOY"& !sp_taxonomicorder == "Unknown"], aes(SL, prey_size)) +
  geom_jitter(width = 0.2, aes(color = Dataset))+
  facet_wrap(~prey_sp, scales = "free_y", ncol = 1) + 
  geom_smooth(method='lm', formula= y~x, aes(color = Dataset, fill = Dataset))+
  labs(x="SL (mm)", y="Prey size")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(Stomach_fish_candat_size, aes(SL, ratio_prey, colour = Dataset)) +
  geom_jitter(width = 0.2)+
  # facet_wrap(~size_class, scales = "free", ncol = 1) + 
  geom_smooth(method='glm', formula= y~x, aes(fill = Dataset), method.args = list(family = "quasipoisson"))+
 stat_poly_eq(use_label(c("eq", "R2"))) +
  labs(x="SL (mm)", y="PPR")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 26),
        legend.title = element_text(size=28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#authors comparison
author <- setDT(readxl::read_xlsx(here::here('Author_ratio.xlsx')))
Stomach_fish_candat_size$author <- "HBU"
Stomach_fish_comparison_size <- rbind(author, Stomach_fish_candat_size[,.(SL,ratio_prey, sp_taxonomicorder, author, prey_sp, prey_size)])
Stomach_fish_comparison_size$author <- factor(Stomach_fish_comparison_size$author, levels = c("HBU", "Keskinen_2004", "Dorner_2007", "Nolan_2018", "Bousseba_2020"))
shapiro.test(Stomach_fish_comparison_size$ratio_prey)
set.seed(3333)
m1 <- lme4::lmer(data = Stomach_fish_comparison_size, formula = ratio_prey ~ SL + (1 + SL | author), REML = F)
summary(m1)
coef(m1)

set.seed(3333)
model <- caret::train(ratio_prey ~ SL + author, Stomach_fish_comparison_size,
               method = "glm",
               family = "quasipoisson",
               trControl = trainControl(method= "repeatedcv",
                                        number = 10, 
                                        verboseIter = T,
                                        repeats = 10))
summary(model)
model

set.seed(3333)
model2 <- caret::train(ratio_prey ~ author, Stomach_fish_comparison_size,
                       method = "glm",
                       family = "quasipoisson",
                      trControl = trainControl(method= "repeatedcv",
                                               number = 10, 
                                               verboseIter = T,
                                               repeats = 10))
summary(model2)
model2
anova(model$finalModel, model2$finalModel)
car::Anova(model$finalModel)

#sp analises
Stomach_fish_comparison_size$ct_catchid <- paste0("Fish_", seq_len(2198))
Stomach_fish_comparison_sp <- dcast(Stomach_fish_comparison_size[prey_sp %in% c("okoun") & author%in% c("HBU","Dorner_2007")],
                                    author + ct_catchid + SL ~ prey_sp, value.var = "ratio_prey")
Stomach_fish_comparison_sp$author <- factor(Stomach_fish_comparison_sp$author, levels = c("HBU","Dorner_2007"))

set.seed(3333)
model_sp <- caret::train(okoun ~ author + SL + author:SL, Stomach_fish_comparison_sp,
                      method = "lm",
                      trControl = trainControl(method= "repeatedcv",
                                               number = 10, 
                                               verboseIter = T,
                                               repeats = 10))
summary(model_sp)
model_sp

set.seed(3333)
model_sp2 <- caret::train(okoun ~ author + SL, Stomach_fish_comparison_sp,
                         method = "lm",
                         trControl = trainControl(method= "repeatedcv",
                                                  number = 10, 
                                                  verboseIter = T,
                                                  repeats = 10))
summary(model_sp2)
model_sp2
anova(model_sp$finalModel, model_sp2$finalModel)

#graphs
ggplot(Stomach_fish_comparison_size, aes(SL, ratio_prey)) +
  geom_jitter(width = 0.2, aes(color = author)) +
  geom_smooth(method='glm', formula= y~x, aes(color = author, fill = author), method.args = list(family = "quasipoisson")) +
  scale_color_manual(values=c(rep("red3",1), rep("black",1), rep("green4", 1), rep("blue1", 1), rep("purple4", 1), rep("grey",1))) +
  scale_fill_manual(values=c(rep("red3",1), rep("black",1), rep("green4", 1), rep("blue1", 1), rep("purple4", 1), rep("grey",1))) +
  labs(x="SL (mm)", y="PPR", fill = "Authors", color = "Authors") + 
  theme(plot.title = element_text(size = 32, face = "bold"),
                      axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
                      axis.text.y = element_text(size = 28),
                      strip.text = element_text(size = 20),
                      axis.title.x = element_text(size = 20),
                      axis.title.y = element_text(size = 26),
                      legend.title = element_text(size=28),
                      legend.text = element_text(size = 26, face = "italic")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom")

ggplot(Stomach_fish_comparison_size[sp_taxonomicorder %in% c("Cypriniformes", "Perciformes")], aes(SL, ratio_prey)) +
  geom_jitter(width = 0.2, aes(color = author))+
  facet_wrap(~sp_taxonomicorder, scales = "free_y", ncol = 1) +
  geom_smooth(method='lm', formula= y~x, aes(color = author, fill = author))+
  scale_color_manual(values=c(rep("red3",1), rep("green4",1), rep("blue1", 1)))+
  scale_fill_manual(values=c(rep("red3",1), rep("green4",1), rep("blue1", 1)))+
  labs(x="SL (mm)", y="PPR", fill = "Authors", color = "Authors")+
  theme(plot.title = element_text(size = 32, face = "bold"),
                      axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
                      axis.text.y = element_text(size = 28),
                      strip.text = element_text(size = 20),
                      axis.title.x = element_text(size = 20),
                      axis.title.y = element_text(size = 26),
                      legend.title = element_text(size=28),
                      legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom")

ggplot(Stomach_fish_comparison_size[prey_sp %in% c("plotice", "jezdik", "candat", "okoun")], aes(SL, ratio_prey)) +
  geom_jitter(width = 0.2, aes(color = author))+
  facet_wrap(~prey_sp, scales = "free_y", ncol = 2) +
  geom_smooth(method='lm', formula= y~x, aes(color = author, fill = author))+
  scale_color_manual(values=c(rep("red3",1), rep("green4",1)))+
  scale_fill_manual(values=c(rep("red3",1), rep("green4",1)))+
  labs(x="SL (mm)", y="PPR", fill = "Authors", color = "Authors")+
  theme(plot.title = element_text(size = 32, face = "bold"),
                      axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
                      axis.text.y = element_text(size = 28),
                      strip.text = element_text(size = 20),
                      axis.title.x = element_text(size = 20),
                      axis.title.y = element_text(size = 26),
                      legend.title = element_text(size=28),
                     legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom")

#prey per 100 fish
stats_sum_prey <- Stomach_fish_candat_melt[,.(fish_n = length(unique(ct_catchid)),
                                              prey_n = sum(prey_n)),
                                              by = .(size_class, Year, sp_taxonomicorder, prey)]
stats_sum_prey[, ':='(prey100 = (prey_n/fish_n)*100)]
stats_sum_prey[, sp_grouped := fct_lump(f = prey, prop = 0.05, w = prey100)]
stats_sum_prey

dec_sum_prey <- Stomach_fish_candat_melt[,.(fish_n = length(unique(ct_catchid)),
                                              prey_n = sum(prey_n)),
                                           by = .(size_class, Dataset, sp_taxonomicorder, prey)]
dec_sum_prey$yv[dec_sum_prey$Dataset == "19's"] <- 5
dec_sum_prey$yv[dec_sum_prey$Dataset == "20's"] <- 9
dec_sum_prey[, ':='(prey100 = (prey_n/fish_n)*100)]
dec_sum_prey[, sp_grouped := fct_lump(f = prey, prop = 0.05, w = prey100)]
dec_sum_prey

ggplot(stats_sum_prey, aes(x = as.factor(Year), y = prey100, fill = sp_taxonomicorder)) +
  geom_col(position = "stack") +
  facet_wrap(~size_class, scales = "free", ncol = 2) + 
  labs(x="Year", y="Prey n per 100 Candat", fill="Prey order")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(dec_sum_prey, aes(x = Dataset, y = prey100, fill = sp_taxonomicorder)) +
  geom_col(position = "stack") +
  facet_wrap(~size_class, scales = "free") + 
  labs(x="Dataset", y="Prey n per 100 Camdat", fill="Prey order")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(stats_sum_prey) +
  geom_col(position = "stack", aes(x = as.factor(Year), y = prey100, fill = prey)) +
  facet_wrap(~size_class, scales = "free", ncol = 2) + 
  scale_fill_manual(values=c(rep("red1",1), rep("black",1), rep("green1", 1), rep("blue1", 1), rep("purple1", 1),
                              rep("green4",1), rep("grey",1), rep("red4", 1), rep("blue4", 1), rep("purple4", 1))) +
  labs(x="Year", y="Prey n per 100 Camdat", fill="Prey sp")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(dec_sum_prey, aes(x = Dataset, y = prey100, fill = prey)) +
  geom_col(position = "stack") +
  facet_wrap(~size_class, scales = "free") + 
  scale_fill_manual(values=c(rep("red1",1), rep("black",1), rep("green1", 1), rep("blue1", 1), rep("purple1", 1),
                             rep("green4",1), rep("grey",1), rep("red4", 1), rep("blue4", 1), rep("purple4", 1))) +
  labs(x="Dataset", y="Prey n per 100 Camdat", fill="Prey sp")+
  theme(plot.title = element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 90, hjust =.1, vjust = .5),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 26),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 26, face = "italic"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



