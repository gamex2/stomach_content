# #Add invertebrate category
# Stomach_content_all <- Stomach_content_all %>% 
#   mutate(invertebrate = case_when(Chironomidae != 0 ~ 1, #condition 1
#                                   Trichoptera != 0 ~ 1,
#                                   Koretra != 0 ~ 1,
#                                   Jepice != 0 ~ 1,
#                                   terr.hmyz != 0 ~ 1,
#                                   crayfish != 0 ~ 1,
#                                   TRUE ~ 0)) #all other cases
# #Add zooplankton category
# Stomach_content_all <- Stomach_content_all %>% 
#   mutate(zooplankton = case_when(Daphnia != 0 ~ 1, #condition 1
#                                  Leptodora != 0 ~ 1,
#                                  Diaphanosoma != 0 ~ 1,
#                                  Cyclopidae != 0 ~ 1,
#                                  TRUE ~ 0)) #all other cases
# 
# #separating data####
# #candat
# Stomach_content_candat <- Stomach_content_all[Species == "candat"]
# # size_lip_yoy <- all_catch_lip[sp_speciesid == "candat" & ct_sl > 1,.(Max = max(ct_sl),
# #                                  Min = min(ct_sl)),
# #                               by =.(year, ct_agegroup)]

# #okoun####
# Stomach_content_okoun <- Stomach_content_all[Species == "okoun"]
# size_lip_ok <- all_catch_lip[sp_speciesid == "okoun" & ct_sl > 1,.(Max = max(ct_sl),
#                                                                      Min = min(ct_sl)),
#                               by =.(year, ct_agegroup)]
# 
# Stomach_content_okoun <- Stomach_content_okoun %>%
#   mutate(size_class = case_when(SL <= 69 ~ "YOY",
#                                 SL %in% 69:200 ~ "older<200",
#                                 SL > 200 ~ "older>200"))
# Stomach_content_okoun$size_class <- factor(Stomach_content_okoun$size_class, levels = c("YOY", "older<200", "older>200"))
# # Stomach_content_candat$ct_agegroup <- dplyr::coalesce(Stomach_content_candat$ct_agegroup, Stomach_content_candat$size_class)
# 
# Stomach_content_okoun_melt <- setDT(melt(Stomach_content_okoun[, .(ct_catchid,Year, Species, Fish, invertebrate, zooplankton, size_class)], 
#                                           id.vars = c("ct_catchid", "Year", "Species", "size_class"), variable.name = "prey", value.name = "prey_n"))
# Stomach_content_okoun_melt <- Stomach_content_okoun_melt[!prey_n==0,]
# 
# ggplot(Stomach_content_okoun_melt, aes(x = as.factor(Year), y = prey_n, fill = prey)) +
#   geom_col(position = "stack") +
#   facet_wrap(~size_class, scales = "free") + 
#   labs(x="Year", y="prey n")+
#   theme(plot.title = element_text(size = 32, face = "bold"),
#         axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
#         axis.text.y = element_text(size = 28),
#         strip.text = element_text(size = 26),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20),
#         legend.title = element_text(size=28),
#         legend.text = element_text(size = 26, face = "italic"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# #fish in the stomach
# Stomach_fish_okoun <- Stomach_content_okoun[Fish==1, ]
# Stomach_fish_okoun <- Stomach_fish_okoun %>%
#   mutate(cannibal = case_when(okoun != 0 ~ 1, #condition 1
#                               TRUE ~ 0)) #all other cases
# # fry_abundance <- setDT(readxl::read_xlsx(here::here('Stomach.content.all.xlsx'), sheet = "fry"))
# Stomach_fish_okoun <- merge(Stomach_fish_okoun, fry_abundance[prey_sp == "okoun"], by = "Year", all.x = T)
# 
# summary(glm(data = Stomach_fish_okoun, formula = cannibal ~ Year+size_class+Fry_abundance, family = binomial(link = "logit")))
# Stomach_fish_okoun2 <- Stomach_fish_okoun
# Stomach_fish_okoun2$cannibal[Stomach_fish_okoun2$cannibal == 0] <- "no"
# Stomach_fish_okoun2$cannibal[Stomach_fish_okoun2$cannibal == 1] <- "yes"
# 
# 
# 
# c.old <- ggplot(Stomach_fish_okoun2[size_class == "older<200"], aes(x = as.factor(Year), y = Fish, fill = forcats::fct_rev(cannibal))) +
#   geom_col(position = "stack") +
#   geom_line(data=Stomach_fish_okoun2[size_class == "older<200"], group=1, mapping=aes(x=as.factor(Year),y=Fry_abundance/180), linewidth = 3)+
#   scale_y_continuous(name = "Fish n",
#                      sec.axis = sec_axis(~.*180, name = "Fry n per hectare"))+
#   labs(x="Year", fill = "cannibalism", linewidth = NULL, title="older<200")+
#   theme(plot.title = element_text(size = 32, face = "bold"),
#         axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
#         axis.text.y = element_text(size = 28),
#         strip.text = element_text(size = 26),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20),
#         legend.title = element_text(size=28),
#         legend.text = element_text(size = 26, face = "italic"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y.right = element_text(vjust=2),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom")
# 
# c.older <- ggplot(Stomach_fish_okoun2[size_class == "older>200"], aes(x = as.factor(Year), y = Fish, fill = forcats::fct_rev(cannibal))) +
#   geom_col(position = "stack") +
#   geom_line(data=Stomach_fish_okoun2[size_class == "older>200"], group=1, mapping=aes(x=as.factor(Year),y=Fry_abundance/300), linewidth = 3)+
#   scale_y_continuous(name = "Fish n",
#                      sec.axis = sec_axis(~.*300, name = "Fry n per hectare"))+
#   labs(x="Year", fill = "cannibalism", linewidth = NULL, title="older>200")+
#   theme(plot.title = element_text(size = 32, face = "bold"),
#         axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
#         axis.text.y = element_text(size = 28),
#         strip.text = element_text(size = 26),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20),
#         legend.title = element_text(size=28),
#         legend.text = element_text(size = 26, face = "italic"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y.right = element_text(vjust=2),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom")
# 
# ggpubr::ggarrange(c.old, c.older,
#                   ncol = 2, nrow = 1,
#                   common.legend = TRUE, legend = "bottom")
# 
# Stomach_fish_okoun_melt <- setDT(melt(Stomach_fish_okoun[, .(ct_catchid, SL, Year, Species, size_class, candat, ouklej,
#                                                                okoun, plotice, jezdik, cejn, cejnek, kaprovitka, okounovitá, Unknown)],
#                                        id.vars = c("ct_catchid", "SL", "Year", "Species", "size_class"),
#                                        variable.name = "prey", value.name = "prey_n"))
# Stomach_fish_okoun_melt <- Stomach_fish_okoun_melt[!prey_n==0,]
# Stomach_fish_okoun_melt <- Stomach_fish_okoun_melt %>% 
#   mutate(cannibal = case_when(prey == "okoun" ~ 1, #condition 1
#                               TRUE ~ 0)) #all other cases
# summary(glm(data = Stomach_fish_okoun_melt, formula = cannibal ~ Year+size_class+prey_n, family = binomial(link = "logit")))
# 
# Stomach_fish_okoun_melt2 <- Stomach_fish_okoun_melt
# Stomach_fish_okoun_melt2$cannibal[Stomach_fish_okoun_melt2$cannibal == 0] <- "no"
# Stomach_fish_okoun_melt2$cannibal[Stomach_fish_okoun_melt2$cannibal == 1] <- "yes"
# 
# ggplot(Stomach_fish_okoun_melt2, aes(x = as.factor(Year), y = prey_n, fill = prey)) +
#   geom_col(position = "stack") +
#   facet_wrap(~size_class, scales = "free") +
#   labs(x="Year", y="Prey n", fill="Cannibalism")+
#   theme(plot.title = element_text(size = 32, face = "bold"),
#         axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
#         axis.text.y = element_text(size = 28),
#         strip.text = element_text(size = 26),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20),
#         legend.title = element_text(size=28),
#         legend.text = element_text(size = 26, face = "italic"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# #fish size in the stomach
# Stomach_fish_okoun_size <- setDT(melt(Stomach_fish_okoun[, .(ct_catchid, SL, Year, Species, size_class, candat_1,candat_2,candat_3,candat_4, candat_5,
#                                                                candat_6, candat_7, candat_8, candat_9, ouklej_1, ouklej_2, ouklej_3, okoun_1,
#                                                                okoun_2, okoun_3, okoun_4, plotice_1, plotice_2, jezdik_1, jezdik_2, jezdik_3,
#                                                                cejn_1, cejnek_1, cejnek_2, kaprovitka_1, kaprovitka_2,Unknown_1, Unknown_2,
#                                                                Unknown_3, Unknown_4)],
#                                        id.vars = c("ct_catchid", "Year", "Species", "SL", "size_class"), variable.name = "prey_sp", value.name = "prey_size"))
# Stomach_fish_okoun_size <- Stomach_fish_okoun_size[!prey_size == 0]
# Stomach_fish_okoun_size[, ':='(ratio_prey = prey_size/SL)]
# Stomach_fish_okoun_size$prey_sp <- sub("\\_.*", "", as.character(Stomach_fish_okoun_size$prey_sp))
# Stomach_fish_okoun_size <- Stomach_fish_okoun_size %>% 
#   mutate(cannibal = case_when(prey_sp == "okoun" ~ 1, #condition 1
#                               TRUE ~ 0)) #all other cases
# summary(glm(data = Stomach_fish_okoun_size, formula = cannibal ~ size_class+Year+ratio_prey, family = binomial(link = "logit")))
# 
# 
# ggplot(Stomach_fish_okoun_size, aes(as.factor(Year), ratio_prey)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.2)+
#   facet_wrap(~size_class) +
#   labs(x="Year", y="Size ratio")+
#   theme(plot.title = element_text(size = 32, face = "bold"),
#         axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
#         axis.text.y = element_text(size = 28),
#         strip.text = element_text(size = 20),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 26),
#         legend.title = element_text(size=28),
#         legend.text = element_text(size = 26, face = "italic"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# ggplot(Stomach_fish_okoun_size, aes(as.factor(Year), ratio_prey)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.2, aes(color = prey_sp, shape = prey_sp), size = 3)+
#   labs(x="Year", y="Size ratio", color="Prey sp", shape="Prey sp")+
#   facet_wrap(~size_class, scales = "free_x") +
#   scale_shape_manual(values = 15:18) +
#   scale_color_manual(values=c(rep("blue3",1), rep("green4",1), rep( "black", 1), rep( "red", 1))) +
#   theme(plot.title = element_text(size = 32, face = "bold"),
#         axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
#         axis.text.y = element_text(size = 28),
#         strip.text = element_text(size = 20),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 26),
#         legend.title = element_text(size=28),
#         legend.text = element_text(size = 26, face = "italic"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))

#stomach content analises
# Stomach_content_candat_melt <- setDT(melt(Stomach_content_all[, .(ct_catchid,Year, Species, Fish, invertebrate, zooplankton, size_class)], 
#                                           id.vars = c("ct_catchid", "Year", "Species", "size_class"), variable.name = "prey", value.name = "prey_n"))
# Stomach_content_candat_melt <- Stomach_content_candat_melt[!prey_n==0,]
# 
# ggplot(Stomach_content_candat_melt, aes(x = as.factor(Year), y = prey_n, fill = prey)) +
#   geom_col(position = "stack") +
#   facet_wrap(~size_class, scales = "free") + 
#   labs(x="Year", y="prey n")+
#   theme(plot.title = element_text(size = 32, face = "bold"),
#         axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
#         axis.text.y = element_text(size = 28),
#         strip.text = element_text(size = 26),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20),
#         legend.title = element_text(size=28),
#         legend.text = element_text(size = 26, face = "italic"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Stomach_fish_candat <- Stomach_content_candat[Fish==1, ]
# Stomach_fish_candat <- Stomach_fish_candat %>% 
#   mutate(cannibal = case_when(candat != 0 ~ 1, #condition 1
#                               TRUE ~ 0)) #all other cases
# 
# Stomach_fish_candat <- merge(Stomach_fish_candat, fry_abundance[prey_sp == "candat"], by = "Year", all.x = T)
# set.seed(3333)
# summary(glm(data = Stomach_fish_candat, formula = cannibal ~ Year+size_class+Fry_abundance, family = binomial(link = "logit")))
# 
# Stomach_fish_candat2 <- Stomach_fish_candat
# Stomach_fish_candat2$cannibal[Stomach_fish_candat2$cannibal == 0] <- "no"
# Stomach_fish_candat2$cannibal[Stomach_fish_candat2$cannibal == 1] <- "yes"
# fry_abundance <- setDT(readxl::read_xlsx(here::here('Stomach.content.all.xlsx'), sheet = "fry"))
# Stomach_fish_candat2 <- merge(Stomach_fish_candat2, fry_abundance[prey_sp == "candat"], by = "Year", all.x = T)
# 
# 
# g.yoy <- ggplot(Stomach_fish_candat2[size_class == "YOY"], aes(x = as.factor(Year), y = Fish, fill = forcats::fct_rev(cannibal))) +
#   geom_col(position = "stack") +
#   geom_line(data=Stomach_fish_candat2[size_class == "YOY"], group=1, mapping=aes(x=as.factor(Year),y=Fry_abundance/20), linewidth = 3)+
#   scale_y_continuous(name = "Fish n", 
#                      sec.axis = sec_axis(~.*20, name = "Fry n per hectare"))+
#   labs(x="Year", fill = "cannibalism", linewidth = NULL, title="YOY")+
#   theme(plot.title = element_text(size = 32, face = "bold"),
#         axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
#         axis.text.y = element_text(size = 28),
#         strip.text = element_text(size = 26),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20),
#         legend.title = element_text(size=28),
#         legend.text = element_text(size = 26, face = "italic"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y.right = element_text(vjust=2),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom")
# 
# g.old <- ggplot(Stomach_fish_candat2[size_class == "older<300"], aes(x = as.factor(Year), y = Fish, fill = forcats::fct_rev(cannibal))) +
#   geom_col(position = "stack") +
#   geom_line(data=Stomach_fish_candat2[size_class == "older<300"], group=1, mapping=aes(x=as.factor(Year),y=Fry_abundance/50), linewidth = 3)+
#   scale_y_continuous(name = "Fish n", 
#                      sec.axis = sec_axis(~.*50, name = "Fry n per hectare"))+
#   labs(x="Year", fill = "cannibalism", linewidth = NULL, title="older<300")+
#   theme(plot.title = element_text(size = 32, face = "bold"),
#         axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
#         axis.text.y = element_text(size = 28),
#         strip.text = element_text(size = 26),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20),
#         legend.title = element_text(size=28),
#         legend.text = element_text(size = 26, face = "italic"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y.right = element_text(vjust=2),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom")
# 
# g.older <- ggplot(Stomach_fish_candat2[size_class == "older>300"], aes(x = as.factor(Year), y = Fish, fill = forcats::fct_rev(cannibal))) +
#   geom_col(position = "stack") +
#   geom_line(data=Stomach_fish_candat2[size_class == "older>300"], group=1, mapping=aes(x=as.factor(Year),y=Fry_abundance/20), linewidth = 3)+
#   scale_y_continuous(name = "Fish n", 
#                      sec.axis = sec_axis(~.*20, name = "Fry n per hectare"))+
#   labs(x="Year", fill = "cannibalism", linewidth = NULL, title="older>300")+
#   theme(plot.title = element_text(size = 32, face = "bold"),
#         axis.text.x = element_text(size = 28,angle = 45, hjust = 1.05, vjust = 1.05),
#         axis.text.y = element_text(size = 28),
#         strip.text = element_text(size = 26),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20),
#         legend.title = element_text(size=28),
#         legend.text = element_text(size = 26, face = "italic"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y.right = element_text(vjust=2),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom")
# 
# ggpubr::ggarrange(g.yoy, g.old, g.older,
#                   ncol = 2, nrow = 2,
#                   common.legend = TRUE, legend = "bottom")

#data stomach####
Stomach_content_sp <- setDT(readxl::read_xlsx(here::here('Stomach.content.all.xlsx'), sheet = "stomach"))
#Replace NA to 0
Stomach_content_sp[, 30:69][is.na(Stomach_content_sp[, 30:69])] <- 0
#standardization of sex
Stomach_content_sp[["Sex"]][is.na(Stomach_content_sp[["Sex"]])] <- "X"
Stomach_content_sp$Sex[Stomach_content_sp$Sex == "Male"] <- "M"
Stomach_content_sp$Sex[Stomach_content_sp$Sex == "Female"] <- "F"
Stomach_content_sp$Sex[Stomach_content_sp$Sex == "?"] <- "X"
Stomach_content_sp <- Stomach_content_sp[!Species == "candat"]
#add age group
Stomach_content_sp <- merge(Stomach_content_sp, all_catch_lip[,.(ct_catchid,ct_agegroup)], by = "ct_catchid", all.x = T)
#add unique id
n_na2 <- sum(is.na(Stomach_content_sp$ct_catchid))
Stomach_content_sp$ct_catchid[is.na(Stomach_content_sp$ct_catchid)] <- paste0("Fish_", seq_len(n_na2))
Stomach_content_sp <- subset( Stomach_content_sp, select = -c(3:8, 13:28) )
colnames(Stomach_content_sp)[4] ="SL"
Stomach_content_bo <- Stomach_content_sp[Species == "bolen"]
Stomach_content_bo <- Stomach_content_bo %>%
  mutate(size_class = case_when(SL <= 120 ~ "YOY",
                                SL %in% 120:300 ~ "older<300",
                                SL > 300 ~ "older>300"))
Stomach_content_bo$size_class <- factor(Stomach_content_bo$size_class, levels = c("YOY", "older<300", "older>300"))

Stomach_bo_candat_melt <- setDT(melt(Stomach_content_bo[, .(ct_catchid, SL, Year, size_class, candat, ouklej, 
                                                               okoun, plotice, jezdik, cejn, cejnek, kaprovitka, okounovitá, Unknown, Fish)], 
                                       id.vars = c("ct_catchid", "SL", "Year", "size_class", "Fish"), 
                                       variable.name = "prey", value.name = "prey_n"))
sum_bo_dec <- Stomach_bo_candat_melt[stomach_c == 0,.(predator_n = uniqueN(ct_catchid)),
                    by =.(Year, size_class)]
sum_bo_dec <- Stomach_bo_candat_melt[,.(prey_n = sum(prey_n),
                                          predator_n = uniqueN(ct_catchid), 
                                          predato_r = sum(prey_n)/uniqueN(ct_catchid)),
                                       by =.(Year, prey, size_class)]
Stomach_bo_candat_size <- setDT(melt(Stomach_content_bo[, .(ct_catchid, SL, Year, Species, candat_1,candat_2,candat_3,candat_4, candat_5,
                                                               candat_6, candat_7, candat_8, candat_9, ouklej_1, ouklej_2, ouklej_3, okoun_1,
                                                               okoun_2, okoun_3, okoun_4, 
                                                               plotice_1, plotice_2, jezdik_1, jezdik_2, jezdik_3,
                                                               cejn_1, cejnek_1, cejnek_2, kaprovitka_1, kaprovitka_2, size_class)],
                                       id.vars = c("ct_catchid", "Year", "Species", "SL", "size_class"), variable.name = "prey_sp", value.name = "prey_size"))
Stomach_bo_candat_size <- Stomach_bo_candat_size[!prey_size == 0]
Stomach_bo_candat_size$prey_sp <- sub("\\_.*", "", as.character(Stomach_bo_candat_size$prey_sp))
sum_bo_dec <- Stomach_bo_candat_size[!prey_size == 0,.(Mean = round(mean(prey_size, na.rm = T), 2),
                                                           SE = round(plotrix::std.error(prey_size), 2),
                                                           Max = max(prey_size),
                                                           Min = min(prey_size),
                                                           N = length(prey_size)),
                                         by =.(size_class, prey_sp)]
write.xlsx(sum_bo_dec, here::here('jan_table.xlsx'))
Stomach_content_su <- Stomach_content_sp[Species == "sumec"]
Stomach_content_su <- Stomach_content_su %>%
  mutate(size_class = case_when(SL %in% 120:300 ~ "older<300",
                                SL > 300 ~ "older>300"))
Stomach_su_candat_melt <- setDT(melt(Stomach_content_su[, .(ct_catchid, SL, Year, size_class, candat, ouklej, 
                                                            okoun, plotice, jezdik, cejn, cejnek, kaprovitka, okounovitá, Unknown, Fish)], 
                                     id.vars = c("ct_catchid", "SL", "Year", "size_class", "Fish"), 
                                     variable.name = "prey", value.name = "prey_n"))
sum_su_dec <- Stomach_su_candat_melt[Fish == 0,.(predator_n = uniqueN(ct_catchid)),
                                     by =.(Year, size_class)]
sum_su_dec <- Stomach_su_candat_melt[,.(prey_n = sum(prey_n),
                                        predator_n = uniqueN(ct_catchid), 
                                        predato_r = sum(prey_n)/uniqueN(ct_catchid)),
                                     by =.(Year, prey, size_class)]
Stomach_su_candat_size <- setDT(melt(Stomach_content_su[, .(ct_catchid, SL, Year, Species, candat_1,candat_2,candat_3,candat_4, candat_5,
                                                            candat_6, candat_7, candat_8, candat_9, ouklej_1, ouklej_2, ouklej_3, okoun_1,
                                                            okoun_2, okoun_3, okoun_4, 
                                                            plotice_1, plotice_2, jezdik_1, jezdik_2, jezdik_3,
                                                            cejn_1, cejnek_1, cejnek_2, kaprovitka_1, kaprovitka_2, size_class)],
                                     id.vars = c("ct_catchid", "Year", "Species", "SL", "size_class"), variable.name = "prey_sp", value.name = "prey_size"))
Stomach_su_candat_size <- Stomach_su_candat_size[!prey_size == 0]
Stomach_su_candat_size$prey_sp <- sub("\\_.*", "", as.character(Stomach_su_candat_size$prey_sp))
sum_su_dec <- Stomach_su_candat_size[!prey_size == 0,.(Mean = round(mean(prey_size, na.rm = T), 2),
                                                       SE = round(plotrix::std.error(prey_size), 2),
                                                       Max = max(prey_size),
                                                       Min = min(prey_size),
                                                       N = length(prey_size)),
                                     by =.(size_class, prey_sp)]
write.xlsx(sum_su_dec, here::here('jan_table.xlsx'))
Stomach_content_su$size_class <- factor(Stomach_content_su$size_class, levels = c("YOY", "older<300", "older>300"))

Stomach_content_ok <- Stomach_content_sp[Species == "okoun"]
Stomach_content_ok <- Stomach_content_ok %>%
  mutate(size_class = case_when(SL <= 70 ~ "YOY",
                                SL %in% 70:300 ~ "older<300",
                                SL > 300 ~ "older>300"))
Stomach_content_ok$size_class <- factor(Stomach_content_ok$size_class, levels = c("YOY", "older<300", "older>300"))
Stomach_ok_candat_melt <- setDT(melt(Stomach_content_ok[, .(ct_catchid, SL, Year, size_class, candat, ouklej, 
                                                            okoun, plotice, jezdik, cejn, cejnek, kaprovitka, okounovitá, Unknown, Fish)], 
                                     id.vars = c("ct_catchid", "SL", "Year", "size_class", "Fish"), 
                                     variable.name = "prey", value.name = "prey_n"))
sum_ok_dec <- Stomach_ok_candat_melt[Fish == 0,.(predator_n = uniqueN(ct_catchid)),
                                     by =.(Year, size_class)]
sum_ok_dec <- Stomach_ok_candat_melt[,.(prey_n = sum(prey_n),
                                        predator_n = uniqueN(ct_catchid), 
                                        predato_r = sum(prey_n)/uniqueN(ct_catchid)),
                                     by =.(Year, prey, size_class)]
Stomach_ok_candat_size <- setDT(melt(Stomach_content_ok[, .(ct_catchid, SL, Year, Species, candat_1,candat_2,candat_3,candat_4, candat_5,
                                                            candat_6, candat_7, candat_8, candat_9, ouklej_1, ouklej_2, ouklej_3, okoun_1,
                                                            okoun_2, okoun_3, okoun_4, 
                                                            plotice_1, plotice_2, jezdik_1, jezdik_2, jezdik_3,
                                                            cejn_1, cejnek_1, cejnek_2, kaprovitka_1, kaprovitka_2, size_class)],
                                     id.vars = c("ct_catchid", "Year", "Species", "SL", "size_class"), variable.name = "prey_sp", value.name = "prey_size"))
Stomach_ok_candat_size <- Stomach_ok_candat_size[!prey_size == 0]
Stomach_ok_candat_size$prey_sp <- sub("\\_.*", "", as.character(Stomach_ok_candat_size$prey_sp))
sum_ok_dec <- Stomach_ok_candat_size[!prey_size == 0,.(Mean = round(mean(prey_size, na.rm = T), 2),
                                                       SE = round(plotrix::std.error(prey_size), 2),
                                                       Max = max(prey_size),
                                                       Min = min(prey_size),
                                                       N = length(prey_size)),
                                     by =.(size_class, prey_sp)]
write.xlsx(sum_ok_dec, here::here('jan_table.xlsx'))

Stomach_content_st <- Stomach_content_sp[Species == "stika"]
Stomach_content_st$size_class <- "older>300"
Stomach_st_candat_melt <- setDT(melt(Stomach_content_st[, .(ct_catchid, SL, Year, size_class, candat, ouklej, 
                                                            okoun, plotice, jezdik, cejn, cejnek, kaprovitka, okounovitá, Unknown, Fish)], 
                                     id.vars = c("ct_catchid", "SL", "Year", "size_class", "Fish"), 
                                     variable.name = "prey", value.name = "prey_n"))
sum_st_dec <- Stomach_st_candat_melt[Fish == 0,.(predator_n = uniqueN(ct_catchid)),
                                     by =.(Year, size_class)]
sum_st_dec <- Stomach_st_candat_melt[,.(prey_n = sum(prey_n),
                                        predator_n = uniqueN(ct_catchid), 
                                        predato_r = sum(prey_n)/uniqueN(ct_catchid)),
                                     by =.(Year, prey, size_class)]
Stomach_st_candat_size <- setDT(melt(Stomach_content_st[, .(ct_catchid, SL, Year, Species, candat_1,candat_2,candat_3,candat_4, candat_5,
                                                            candat_6, candat_7, candat_8, candat_9, ouklej_1, ouklej_2, ouklej_3, okoun_1,
                                                            okoun_2, okoun_3, okoun_4, 
                                                            plotice_1, plotice_2, jezdik_1, jezdik_2, jezdik_3,
                                                            cejn_1, cejnek_1, cejnek_2, kaprovitka_1, kaprovitka_2, size_class)],
                                     id.vars = c("ct_catchid", "Year", "Species", "SL", "size_class"), variable.name = "prey_sp", value.name = "prey_size"))
Stomach_st_candat_size <- Stomach_st_candat_size[!prey_size == 0]
Stomach_st_candat_size$prey_sp <- sub("\\_.*", "", as.character(Stomach_st_candat_size$prey_sp))
sum_st_dec <- Stomach_st_candat_size[!prey_size == 0,.(Mean = round(mean(prey_size, na.rm = T), 2),
                                                       SE = round(plotrix::std.error(prey_size), 2),
                                                       Max = max(prey_size),
                                                       Min = min(prey_size),
                                                       N = length(prey_size)),
                                     by =.(Year, prey_sp, size_class)]
write.xlsx(sum_st_dec, here::here('jan_table.xlsx'))

