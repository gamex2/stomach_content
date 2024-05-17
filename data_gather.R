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
  mutate(decade = case_when(Year <= 1970 ~ "19's",
                            Year > 2009 ~ "20's"))

ggplot(Stomach_content_all, aes(x = as.factor(Year))) + 
  geom_histogram(stat="count")


#fish in the stomach
Stomach_fish_candat_melt <- setDT(melt(Stomach_content_all[, .(ct_catchid, SL, Year, Species, size_class, candat, ouklej, 
                                                               okoun, plotice, jezdik, cejn, cejnek, kaprovitka, okounovitá, Unknown, decade)], 
                                          id.vars = c("ct_catchid", "SL", "Year", "Species", "size_class", "decade"), 
                                       variable.name = "prey", value.name = "prey_n"))
# Stomach_fish_candat_melt <- Stomach_fish_candat_melt[!prey_n==0,]
Stomach_fish_candat_melt[, sp_grouped := fct_lump(f = prey, prop = 0.05, w = prey_n)]
Stomach_fish_candat_melt <- merge(Stomach_fish_candat_melt, specs[,.(sp_speciesid,sp_taxonomicorder)], by.x = "prey", by.y = "sp_speciesid", all.x = T)
Stomach_fish_candat_melt$sp_taxonomicorder[Stomach_fish_candat_melt$prey == "kaprovitka"] <- "Cypriniformes"
Stomach_fish_candat_melt$sp_taxonomicorder[Stomach_fish_candat_melt$prey == "okounovitá"] <- "Perciformes"
Stomach_fish_candat_melt$sp_taxonomicorder[Stomach_fish_candat_melt$prey == "Unknown"] <- "Unknown"

set.seed(3333)
m1 <- lm(data = Stomach_fish_candat_melt, formula = prey_n ~ sp_taxonomicorder+SL+Year+sp_taxonomicorder:SL+sp_taxonomicorder:Year)
m2 <- lm(data = Stomach_fish_candat_melt, formula = prey_n ~ sp_taxonomicorder+SL+Year)
anova(m2)
summary(m2)
with(summary(m2), 1 - deviance/null.deviance)
anova(m1, m2)
m3 <- lm(data = Stomach_fish_candat_melt, formula = prey_n ~ SL+Year)
anova(m2, m3)

summary(aov(data = Stomach_fish_candat_melt, formula = prey_n ~ (size_class+Year)*sp_taxonomicorder))

ggplot(Stomach_fish_candat_melt, aes(x = as.factor(Year), y = prey_n, fill = sp_taxonomicorder)) +
  geom_col(position = "stack") +
  facet_wrap(~size_class, scales = "free", ncol = 2) + 
  labs(x="Year", y="Prey n", fill="Prey order")+
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

ggplot(Stomach_fish_candat_melt, aes(x = decade, y = (sum(prey_n)/length(unique(ct_catchid)))*100, fill = sp_taxonomicorder)) +
  geom_col(position = "stack") +
  facet_wrap(~size_class, scales = "free") + 
  labs(x="Year", y="Prey n", fill="Prey order")+
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


#fish size in the stomach
Stomach_fish_candat_size <- setDT(melt(Stomach_content_all[, .(ct_catchid, SL, Year, Species, size_class, candat_1,candat_2,candat_3,candat_4, candat_5,
                                                              candat_6, candat_7, candat_8, candat_9, ouklej_1, ouklej_2, ouklej_3, okoun_1,
                                                              okoun_2, okoun_3, okoun_4, okoun_5, okoun_6, okoun_7, okoun_8, okoun_9, okoun_10, okoun_11, 
                                                              plotice_1, plotice_2, jezdik_1, jezdik_2, jezdik_3, jezdik_4, jezdik_5, jezdik_6, jezdik_7,
                                                              cejn_1, cejnek_1, cejnek_2, kaprovitka_1, kaprovitka_2, decade)], 
                                     id.vars = c("ct_catchid", "Year", "Species", "SL", "size_class", "decade"), variable.name = "prey_sp", value.name = "prey_size"))
Stomach_fish_candat_size <- Stomach_fish_candat_size[!prey_size == 0]
Stomach_fish_candat_size[, ':='(ratio_prey = prey_size/SL)]
Stomach_fish_candat_size$prey_sp <- sub("\\_.*", "", as.character(Stomach_fish_candat_size$prey_sp))
Stomach_fish_candat_size$prey_sp <- factor(Stomach_fish_candat_size$prey_sp, levels = c("candat", "okoun", "jezdik", "cejn", "cejnek", "kaprovitka", "ouklej", "plotice", "Unknown"))
Stomach_fish_candat_size <- merge(Stomach_fish_candat_size, specs[,.(sp_speciesid,sp_taxonomicorder)], by.x = "prey_sp", by.y = "sp_speciesid", all.x = T)
Stomach_fish_candat_size$sp_taxonomicorder[Stomach_fish_candat_size$prey_sp == "kaprovitka"] <- "Cypriniformes"
Stomach_fish_candat_size$sp_taxonomicorder[Stomach_fish_candat_size$prey_sp == "okounovitá"] <- "Perciformes"
Stomach_percid_size <- Stomach_fish_candat_size[sp_taxonomicorder == "Perciformes" & !size_class == "YOY"]
Stomach_cyprinid_size <- Stomach_fish_candat_size[sp_taxonomicorder == "Cypriniformes" & !size_class == "YOY"]

set.seed(3333)
summary(lm(data = Stomach_cyprinid_size, formula = prey_size ~ SL + decade))
summary(lm(data = Stomach_percid_size, formula = ratio_prey ~ SL + decade))
summary(stats::aov(data = Stomach_fish_candat_size, formula = ratio_prey ~ SL+decade+sp_taxonomicorder))

#cyprinid
set.seed(3333)
m1 <- lme4::lmer(data = Stomach_cyprinid_size, formula = ratio_prey ~ SL + Year + (1|ct_catchid), REML = FALSE)
m2 <- lme4::lmer(data = Stomach_cyprinid_size, formula = ratio_prey ~ SL + (1|ct_catchid), REML = FALSE)
# with(summary(m1), 1 - deviance/null.deviance)
anova(m1, m2)
m3 <- lme4::lmer(data = Stomach_cyprinid_size, formula = ratio_prey ~ SL + (1 + SL|ct_catchid), REML = FALSE)
anova(m2, m3)
car::Anova(m1)
summary(m1)

#percid
set.seed(3333)
m4 <- lme4::lmer(data = Stomach_percid_size, formula = ratio_prey ~ SL + Year + (1|ct_catchid)+(1|size_class), REML = FALSE)
m5 <- lme4::lmer(data = Stomach_percid_size, formula = ratio_prey ~ SL + (1|ct_catchid)+(1|size_class), REML = FALSE)
# with(summary(m1), 1 - deviance/null.deviance)
anova(m4, m5)
m6 <- lme4::lmer(data = Stomach_percid_size, formula = ratio_prey ~ Year + (1|ct_catchid), REML = FALSE)
anova(m4, m6)
car::Anova(m6)
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


ggplot(Stomach_fish_candat_size, aes(decade, ratio_prey)) +
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

ggplot(Stomach_fish_candat_size, aes(decade, ratio_prey)) +
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

ggplot(Stomach_fish_candat_size[!size_class == "YOY"& prey_sp %in% c("okoun", "jezdik", "plotice")], aes(SL, prey_size)) +
  geom_jitter(width = 0.2, aes(color = decade))+
  facet_wrap(~prey_sp, scales = "free_y", ncol = 1) + 
  geom_smooth(method='lm', formula= y~x, aes(color = decade, fill = decade))+
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

ggplot(Stomach_fish_candat_size[!size_class == "YOY"& prey_sp %in% c("okoun", "jezdik", "plotice")], aes(SL, ratio_prey)) +
  geom_jitter(width = 0.2, aes(color = decade))+
  facet_wrap(~prey_sp, scales = "free_y", ncol = 1) + 
  geom_smooth(method='lm', formula= y~x, aes(color = decade, fill = decade))+
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
  geom_jitter(width = 0.2, aes(color = decade))+
  facet_wrap(~prey_sp, scales = "free_y", ncol = 1) + 
  geom_smooth(method='lm', formula= y~x, aes(color = decade, fill = decade))+
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

ggplot(Stomach_fish_candat_size[!size_class == "YOY"& !sp_taxonomicorder == "Unknown"], aes(SL, ratio_prey)) +
  geom_jitter(width = 0.2, aes(color = decade))+
  facet_wrap(~sp_taxonomicorder, scales = "free_y", ncol = 1) + 
  geom_smooth(method='lm', formula= y~x, aes(color = decade, fill = decade))+
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

#authors comparison
author <- setDT(readxl::read_xlsx(here::here('Author_ratio.xlsx')))
Stomach_fish_candat_size$author <- "HBU"
Stomach_fish_comparison_size <- rbind(author, Stomach_fish_candat_size[,.(SL,ratio_prey, sp_taxonomicorder, author, prey_sp, prey_size)])
Stomach_fish_comparison_size$author <- factor(Stomach_fish_comparison_size$author, levels = c("HBU", "Bousseba_2020", "Keskinen_2004", "Nolan_2018", "Dorner_2007"))

set.seed(3333)
m1 <- lm(data = Stomach_fish_comparison_size, formula = ratio_prey ~ SL + author + SL:author)
m2 <- lme4::lmer(data = Stomach_fish_comparison_size, formula = ratio_prey ~ 1)
# with(summary(m1), 1 - deviance/null.deviance)
anova(m1, m2)
car::Anova(m1)
summary(m2)

set.seed(3333)
model <- caret::train(ratio_prey ~ SL + author, Stomach_fish_comparison_size,
               method = "lm",
               trControl = trainControl(method= "repeatedcv",
                                        number = 10, 
                                        verboseIter = T,
                                        repeats = 10))
summary(model)
model

set.seed(3333)
model2 <- caret::train(ratio_prey ~ author, Stomach_fish_comparison_size,
                      method = "lm",
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
  geom_smooth(method='loess', formula= y~x, aes(color = author, fill = author)) +
  scale_color_manual(values=c(rep("red3",1), rep("black",1), rep("green4", 1), rep("blue1", 1), rep("purple4", 1))) +
  scale_fill_manual(values=c(rep("red3",1), rep("black",1), rep("green4", 1), rep("blue1", 1), rep("purple4", 1))) +
  labs(x="SL (mm)", y="Prey ratio", fill = "Authors", color = "Authors") + 
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
  scale_color_manual(values=c(rep("black",1), rep("green4",1), rep("purple4", 1)))+
  scale_fill_manual(values=c(rep("black",1), rep("green4",1), rep("purple4", 1)))+
  labs(x="SL (mm)", y="Prey ratio", fill = "Authors", color = "Authors")+
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
  scale_color_manual(values=c(rep("black",1), rep("green4",1)))+
  scale_fill_manual(values=c(rep("black",1), rep("green4",1)))+
  labs(x="SL (mm)", y="Prey ratio", fill = "Authors", color = "Authors")+
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
                                           by = .(size_class, decade, sp_taxonomicorder, prey)]
dec_sum_prey$yv[dec_sum_prey$decade == "19's"] <- 5
dec_sum_prey$yv[dec_sum_prey$decade == "20's"] <- 9
dec_sum_prey[, ':='(prey100 = (prey_n/fish_n)*100)]
dec_sum_prey[, sp_grouped := fct_lump(f = prey, prop = 0.05, w = prey100)]
dec_sum_prey

ggplot(stats_sum_prey, aes(x = as.factor(Year), y = prey100, fill = sp_taxonomicorder)) +
  geom_col(position = "stack") +
  facet_wrap(~size_class, scales = "free", ncol = 2) + 
  labs(x="Year", y="Prey n per 100 Camdat", fill="Prey order")+
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

ggplot(dec_sum_prey, aes(x = decade, y = prey100, fill = sp_taxonomicorder)) +
  geom_col(position = "stack") +
  facet_wrap(~size_class, scales = "free") + 
  labs(x="Decade", y="Prey n per 100 Camdat", fill="Prey order")+
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

ggplot(dec_sum_prey, aes(x = decade, y = prey100, fill = prey)) +
  geom_col(position = "stack") +
  facet_wrap(~size_class, scales = "free") + 
  scale_fill_manual(values=c(rep("red1",1), rep("black",1), rep("green1", 1), rep("blue1", 1), rep("purple1", 1),
                             rep("green4",1), rep("grey",1), rep("red4", 1), rep("blue4", 1), rep("purple4", 1))) +
  labs(x="Decade", y="Prey n per 100 Camdat", fill="Prey sp")+
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



