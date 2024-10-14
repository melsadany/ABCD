################################################################################
#                       variables predicting ASD status/PGS                    #
################################################################################
rm(list = ls())
gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
################################################################################
################################################################################
project.dir <- paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                      "/jmichaelson-wdata/msmuhammad/projects/ABCD/ASD")
setwd(project.dir)
################################################################################
################################################################################
################################################################################
# load data and reformat per category for both ABCD and RPOE
# 

########
# demo #
abcd.demo <- read_csv("../../../data/ABCD/abcd5/age-sex-by-eventname.csv") %>%
  filter(grepl("baseline", eventname))
abcd.dx <- read_csv(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                           "/jmichaelson-sdata/ABCD/abcd_release_5_0/core/abcd-general/abcd_p_screen.csv")) %>%
  select(IID = src_subject_id, ASD_dx = scrn_asd)
rpoe.demo <- read_csv("../../RPOE/shared_data/data/demo-full.csv")

#######
# cog #
abcd.cog <- read_csv(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                            "/jmichaelson-sdata/ABCD/abcd_release_5_0/core/neurocognition/nc_y_nihtb.csv")) %>%
  filter(grepl("baseline", eventname)) %>%
  select(IID = src_subject_id, eventname, ends_with("agecorrected")) %>%
  inner_join(abcd.demo)

rpoe.cog <- read_rds("../../RPOE/shared_data/data/m1m2-sex-corrected.rds") %>%
  select(te_id,
         nihtbx_picvocab_agecorrected = picture_vocabulary_age_corrected_standard_score,
         nihtbx_flanker_agecorrected = flanker_inhibitory_control_age_corrected_standard_score,
         nihtbx_list_agecorrected = list_sorting_wm_age_corrected_standard_score,
         nihtbx_cardsort_agecorrected = dimensional_change_card_sort_age_corrected_standard_score,
         nihtbx_pattern_agecorrected = pattern_comparison_PS_age_corrected_standard_score,
         nihtbx_picture_agecorrected = picture_sequence_memory_test_age_corrected_standard_score,
         nihtbx_reading_agecorrected = oral_reading_recognition_age_corrected_standard_score)
abcd.cogm <- abcd.cog %>%
  mutate_at(.vars = colnames(abcd.cog)[-c(1,2,13,14)],
            .funs = function(x){
              df <- abcd.cog %>% 
                mutate(y = x) %>%
                select(y, sex)
              z_from_lm(y = df$y, x = df[,-1])
            })
abcd.cog2 <- abcd.cogm %>%
  select(IID, eventname, colnames(rpoe.cog)[-1]) %>%
  drop_na()
rm(abcd.cog)

#######
# MRI #

#   structural
abcd.s.mri <- read_csv(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                              "/jmichaelson-sdata/ABCD/abcd_release_5_0/core/imaging/mri_y_smr_vol_aseg.csv")) %>% 
  filter(grepl("baseline", eventname)) %>%
  rowwise() %>%
  mutate(total_volume = sum(c_across(-c(1, 2)), na.rm = TRUE)) %>%
  ungroup()
rpoe.s.mri <- read_rds(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                              "/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data",
                              "/derivatives/brain-vols-synthreg-age-sex-tot-vol-corrected.rds"))
colnames(rpoe.s.mri) <- c("te_id",
                          "smri_vol_scs_cbwmatterlh", colnames(rpoe.s.mri)[3],
                          "smri_vol_scs_ltventriclelh","smri_vol_scs_inflatventlh","smri_vol_scs_crbwmatterlh",
                          "smri_vol_scs_crbcortexlh","smri_vol_scs_tplh","smri_vol_scs_caudatelh",    
                          "smri_vol_scs_putamenlh","smri_vol_scs_pallidumlh","smri_vol_scs_3rdventricle", 
                          "smri_vol_scs_4thventricle","smri_vol_scs_bstem","smri_vol_scs_hpuslh",   
                          "smri_vol_scs_amygdalalh","smri_vol_scs_csf",
                          "smri_vol_scs_aal","smri_vol_scs_vedclh",
                          
                          "smri_vol_scs_cbwmatterrh", colnames(rpoe.s.mri)[21],
                          "smri_vol_scs_ltventriclerh","smri_vol_scs_inflatventrh","smri_vol_scs_crbwmatterrh",
                          "smri_vol_scs_crbcortexrh","smri_vol_scs_tprh","smri_vol_scs_caudaterh",    
                          "smri_vol_scs_putamenrh","smri_vol_scs_pallidumrh","smri_vol_scs_hpusrh",
                          "smri_vol_scs_amygdalarh","smri_vol_scs_aar",
                          "smri_vol_scs_vedcrh",
                          
                          colnames(rpoe.s.mri)[c(34:40)])
rpoe.s.mri2 <- rpoe.s.mri %>% select(te_id, any_of(colnames(abcd.s.mri)))
abcd.s.mri2 <- abcd.s.mri %>%
  mutate_at(.vars = -c(1,2, ncol(.data)), 
            .funs = function(x) x/abcd.s.mri$total_volume) %>%
  select(IID = src_subject_id, any_of(colnames(rpoe.s.mri2))) %>%
  inner_join(abcd.demo)
abcd.s.mri3 <- abcd.s.mri2 %>%
  mutate_at(.vars = vars(starts_with("smri")),
            .funs = function(x){
              df <- abcd.s.mri2 %>% 
                mutate(y = x) %>%
                select(y, interview_age, sex)
              z_from_lm(y = df$y, x = df[,-1])
            })
rm(rpoe.s.mri);rm(abcd.s.mri);rm(abcd.s.mri2)


#   diffusion
#       FA
abcd.d.mri.fa <- read_csv(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                                 "/jmichaelson-sdata/ABCD/abcd_release_5_0/core/imaging/mri_y_dti_fa_fs_at.csv")) %>% 
  filter(grepl("baseline", eventname)) %>%
  rename(ProjectionBasalGanglia_ThalamicRadiationL = dmdtifp1_10,
         ProjectionBasalGanglia_ThalamicRadiationR = dmdtifp1_9,
         Commissure_CorpusCallosum = dmdtifp1_19,
         ProjectionBrainstem_CorticospinalTractL = dmdtifp1_8,
         ProjectionBrainstem_CorticospinalTractR = dmdtifp1_7,
         ProjectionBasalGanglia_FornixL = dmdtifp1_2,
         ProjectionBasalGanglia_FornixR = dmdtifp1_1,
         Association_InferiorFrontoOccipitalFasciculusL = dmdtifp1_16,
         Association_InferiorFrontoOccipitalFasciculusR = dmdtifp1_15,
         Association_InferiorLongitudinalFasciculusL = dmdtifp1_14,
         Association_InferiorLongitudinalFasciculusR = dmdtifp1_13,
         ProjectionBasalGanglia_CorticostriatalTractL = dmdtifp1_27,
         ProjectionBasalGanglia_CorticostriatalTractR = dmdtifp1_26,
         Association_SuperiorLongitudinalFasciculusL = dmdtifp1_21,
         Association_SuperiorLongitudinalFasciculusR = dmdtifp1_20,
         Association_UncinateFasciculusL = dmdtifp1_12,
         Association_UncinateFasciculusR = dmdtifp1_11) %>%
  select(IID = src_subject_id, starts_with(c("Projection", "Comm", "Asso"))) %>%
  inner_join(abcd.demo)
abcd.d.mri.fa2 <- abcd.d.mri.fa %>%
  mutate_at(.vars = vars(starts_with(c("Projection", "Comm", "Asso"))),
            .funs = function(x){
              df <- abcd.d.mri.fa %>% 
                mutate(y = x) %>%
                select(y, interview_age, sex)
              z_from_lm(y = df$y, x = df[,-1])
            })

rpoe.d.mri <- read_rds(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                              "/jmichaelson-wdata/msmuhammad/projects/RPOE/dwi/data",
                              "/derivatives/recognized-tracts-summary-stats-DSI_RAW.rds"))
rpoe.d.mri.fa <- rpoe.d.mri %>%
  filter(roi %in% colnames(abcd.d.mri.fa2),
         feature %in% c("dti_fa")) %>%
  mutate(new_name = paste0(feature, "__", roi)) %>%
  pivot_wider(names_from = new_name, values_from = value, id_cols = c("te_id")) %>%
  select(te_id, any_of(paste0(c("dti_fa__"), colnames(abcd.d.mri.fa2)))) %>%
  left_join(rpoe.demo %>% select(te_id, MRI_age, sex)) %>%
  drop_na(MRI_age) 
fa.keep <- names(colSums(rpoe.d.mri.fa %>% select(starts_with("dti")))[!is.na(colSums(rpoe.d.mri.fa %>% select(starts_with("dti"))))])
rpoe.d.mri.fa2 <- rpoe.d.mri.fa %>%
  select(te_id, MRI_age, sex, fa.keep)

rpoe.d.mri.fa2 %>%
  ggplot(aes(x = dti_fa__Commissure_CorpusCallosum, y = MRI_age)) +
  geom_point(shape = 1) +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(method = "spearman") +
  bw.theme

rpoe.d.mri.fa3 <- rpoe.d.mri.fa2 %>%
  mutate_at(.vars = fa.keep,
            .funs = function(x) {
              df <- rpoe.d.mri.fa2 %>%
                mutate(y = x) %>%
                select(y, MRI_age, sex)
              z_from_lm(y = df$y, x = df[,-1])
            })


abcd.d.mri.fa3 <- abcd.d.mri.fa2 %>% 
  filter(grepl("baseline", eventname)) %>%
  rename_at(.vars = vars(starts_with(c("Projection", "Comm", "Asso"))), 
            .funs = function(x) paste0("dti_fa__", x)) %>%
  select(IID, any_of(fa.keep))

#       MD
abcd.d.mri.md <- read_csv(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                                 "/jmichaelson-sdata/ABCD/abcd_release_5_0/core/imaging/mri_y_dti_md_fs_at.csv")) %>% 
  filter(grepl("baseline", eventname)) %>%
  rename(ProjectionBasalGanglia_ThalamicRadiationL = dmdtifp1_52,
         ProjectionBasalGanglia_ThalamicRadiationR = dmdtifp1_51,
         Commissure_CorpusCallosum = dmdtifp1_61,
         ProjectionBrainstem_CorticospinalTractL = dmdtifp1_50,
         ProjectionBrainstem_CorticospinalTractR = dmdtifp1_49,
         ProjectionBasalGanglia_FornixL = dmdtifp1_44,
         ProjectionBasalGanglia_FornixR = dmdtifp1_43,
         Association_InferiorFrontoOccipitalFasciculusL = dmdtifp1_58,
         Association_InferiorFrontoOccipitalFasciculusR = dmdtifp1_57,
         Association_InferiorLongitudinalFasciculusL = dmdtifp1_56,
         Association_InferiorLongitudinalFasciculusR = dmdtifp1_55,
         ProjectionBasalGanglia_CorticostriatalTractL = dmdtifp1_69,
         ProjectionBasalGanglia_CorticostriatalTractR = dmdtifp1_68,
         Association_SuperiorLongitudinalFasciculusL = dmdtifp1_63,
         Association_SuperiorLongitudinalFasciculusR = dmdtifp1_62,
         Association_UncinateFasciculusL = dmdtifp1_54,
         Association_UncinateFasciculusR = dmdtifp1_53) %>%
  select(IID = src_subject_id, starts_with(c("Projection", "Comm", "Asso"))) %>%
  inner_join(abcd.demo)
abcd.d.mri.md2 <- abcd.d.mri.md %>%
  mutate_at(.vars = vars(starts_with(c("Projection", "Comm", "Asso"))),
            .funs = function(x){
              df <- abcd.d.mri.md %>% 
                mutate(y = x) %>%
                select(y, interview_age, sex)
              z_from_lm(y = df$y, x = df[,-1])
            })

rpoe.d.mri.md <- rpoe.d.mri %>%
  filter(roi %in% colnames(abcd.d.mri.md2),
         feature %in% c("md")) %>%
  mutate(new_name = paste0(feature, "__", roi)) %>%
  pivot_wider(names_from = new_name, values_from = value, id_cols = c("te_id")) %>%
  select(te_id, any_of(paste0(c("md__"), colnames(abcd.d.mri.md2)))) %>%
  left_join(rpoe.demo %>% select(te_id, MRI_age, sex)) %>%
  drop_na(MRI_age) 
md.keep <- names(colSums(rpoe.d.mri.md %>% select(starts_with("md")))[!is.na(colSums(rpoe.d.mri.md %>% select(starts_with("md"))))])
rpoe.d.mri.md2 <- rpoe.d.mri.md %>%
  select(te_id, MRI_age, sex, md.keep)

rpoe.d.mri.md2 %>%
  ggplot(aes(x = md__Commissure_CorpusCallosum, y = MRI_age)) +
  geom_point(shape = 1) +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(method = "spearman") +
  bw.theme

rpoe.d.mri.md3 <- rpoe.d.mri.md2 %>%
  mutate_at(.vars = md.keep,
            .funs = function(x) {
              df <- rpoe.d.mri.md2 %>%
                mutate(y = x) %>%
                select(y, MRI_age, sex)
              z_from_lm(y = df$y, x = df[,-1])
            })

abcd.d.mri.md3 <- abcd.d.mri.md2 %>% 
  filter(grepl("baseline", eventname)) %>%
  rename_at(.vars = vars(starts_with(c("Projection", "Comm", "Asso"))), 
            .funs = function(x) paste0("md__", x)) %>%
  select(IID, any_of(md.keep))

rm(abcd.d.mri.fa);rm(abcd.d.mri.fa2);rm(abcd.d.mri.md);rm(abcd.d.mri.md2);rm(rpoe.d.mri)
rm(rpoe.d.mri.fa2);rm(rpoe.d.mri.md2);rm(rpoe.d.mri.fa);rm(rpoe.d.mri.md)


#   fMRI
#     fALFF
roi.meta <- read_csv(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                            "/jmichaelson-wdata/msmuhammad/refs/Schaefer2018/Parcellations/MNI/Centroid_coordinates/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv")) %>%
  mutate(full_network = sub("7Networks_.H_", "", `ROI Name`),
         network = sub("_.*", "", full_network),
         hemisphere = ifelse(grepl("_LH_", `ROI Name`), "LH", "RH"),
         h_network = paste0(hemisphere, "_", network),
         roi_name=sub("7Networks_", "", `ROI Name`),
         net9 = sub("[ABC]$", "", network),
         h_net9 = paste0(hemisphere, "_", net9))
rpoe.f.mri.falff <- read_rds("../../RPOE/mri/data/derivatives/R-func/REST-roi-fALFF-del5-100-7.rds") %>%
  pivot_longer(cols = -c(1), names_to = "roi_name") %>%
  left_join(roi.meta) %>%
  group_by(te_id, network) %>%
  dplyr::summarise(value = mean(value, na.rm = T)) %>% ungroup() %>%
  pivot_wider(names_from = "network", values_from = "value", id_cols = "te_id") %>%
  left_join(rpoe.demo %>% select(te_id, MRI_age, sex)) 
rpoe.f.mri.falff2 <- rpoe.f.mri.falff %>%
  mutate_at(.vars = roi.meta$network, .funs = function(x) {
    df <- rpoe.f.mri.falff %>% mutate(y = x) %>%
      select(y, MRI_age, sex)
    z_from_lm(y = df$y, x = df[,-1])
  }) %>% select(-c(MRI_age, sex)) %>%
  rename_at(.vars = vars(-c(te_id)), .funs = function(x) paste0("fALFF__", x))
rm(rpoe.f.mri.falff)


#     rest ReHo
rpoe.f.mri.reho <- read_rds("../../RPOE/mri/data/derivatives/R-func/participants-REST-roi-raw-ReHo-100-7.rds") %>%
  group_by(te_id, network) %>%
  dplyr::summarise(value = mean(ReHo, na.rm = T)) %>% ungroup() %>%
  pivot_wider(names_from = "network", values_from = "value", id_cols = "te_id") %>%
  left_join(rpoe.demo %>% select(te_id, MRI_age, sex)) 
rpoe.f.mri.reho2 <- rpoe.f.mri.reho %>%
  mutate_at(.vars = roi.meta$network, .funs = function(x) {
    df <- rpoe.f.mri.reho %>% mutate(y = x) %>%
      select(y, MRI_age, sex)
    z_from_lm(y = df$y, x = df[,-1])
  }) %>% select(-c(MRI_age, sex)) %>%
  rename_at(.vars = vars(-c(te_id)), .funs = function(x) paste0("ReHo__", x))
rm(rpoe.f.mri.reho)



#################
# mental health #
abcd.mh <- read_csv(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                           "/jmichaelson-sdata/ABCD/abcd_release_5_0/core/mental-health/mh_p_cbcl.csv")) %>%
  select(IID = src_subject_id,2, contains("syn"), contains("dsm5"))
abcd.mh2 <- abcd.mh[,c(T,T,grepl("_r$", colnames(abcd.mh)[-c(1:2)]))] %>%
  filter(grepl("baseline", eventname)) %>% drop_na() %>%
  inner_join(abcd.demo) %>%
  rename_all(.funs = function(x) sub("cbcl_scr_", "", x))
abcd.mh3 <- abcd.mh2 %>%
  mutate_at(.vars = vars(starts_with(c("syn", "dsm5"))),
            .funs = function(x){
              df <- abcd.mh2 %>% 
                mutate(y = x) %>%
                select(y, interview_age, sex)
              z_from_lm(y = df$y, x = df[,-1])
            })

rpoe.mh <- read_csv(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                           "/jmichaelson-wdata/msmuhammad/projects/RPOE/behavior/data/derivatives/all-behavior-data-age-sex-corrected.csv"))  %>%
  select(devGenes_id,
         syn_anxdep_r = syn_anxious,
         syn_withdep_r = syn_withdrawn,
         syn_somatic_r = syn_somatic, 
         syn_social_r = syn_social,
         syn_thought_r = syn_thought,
         syn_attention_r = syn_attention,
         syn_rulebreak_r = syn_rulebreaking,
         syn_aggressive_r = syn_aggressive,
         dsm5_depress_r = dsm5_depressive,
         dsm5_anxdisord_r = dsm5_anxiety,
         dsm5_somaticpr_r = dsm5_somatic, 
         dsm5_adhd_r = dsm5_ADHD)
abcd.mh4 <- abcd.mh3 %>% 
  filter(grepl("baseline", eventname)) %>%
  select(IID, any_of(colnames(rpoe.mh)))

rm(abcd.mh);rm(abcd.mh2);rm(abcd.mh3)

#######
# PGS #
abcd.pgs <- read_tsv(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
                            "/jmichaelson-wdata/trthomas/array/merged_2022_ABCD_iWES1_WGS_2-4/PGS/all_PGS.tsv")) %>%
  filter(grepl("autism|cog_performance-SSGAC-2018", PGS_name), IID %in% abcd.demo$IID) %>%
  pivot_wider(names_from = PGS_name, values_from = PGS_std_PCresid, id_cols = "IID") %>%
  select(IID, ASD = `autism-PGC-2019`, CP = `cog_performance-SSGAC-2018`)

abcd.pgs %>%
  inner_join(abcd.dx) %>%
  mutate(ASD_dx = as.factor(ASD_dx)) %>%
  drop_na() %>%
  ggplot(aes(x=ASD, color = ASD_dx)) +
  geom_density()
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# combine all data

abcd.all <- inner_join(abcd.cog2,
                       abcd.d.mri.fa3) %>%
  inner_join(abcd.d.mri.md3) %>%
  inner_join(abcd.s.mri3) %>%
  inner_join(abcd.mh4) %>%
  inner_join(abcd.pgs) %>%
  inner_join(abcd.dx %>% mutate(ASD_dx = as.factor(ASD_dx))) %>%
  select(-c(interview_age, sex, eventname)) %>%
  drop_na(ASD_dx)

rpoe.all <- inner_join(rpoe.cog,
                       rpoe.d.mri.fa3) %>%
  inner_join(rpoe.d.mri.md3) %>%
  inner_join(rpoe.s.mri2) %>%
  inner_join(rpoe.mh %>% inner_join(rpoe.demo %>% select(te_id, devGenes_id))) %>%
  inner_join(rpoe.demo %>% select(te_id, ASD_dx)) %>%
  select(-c(devGenes_id, sex, MRI_age)) %>%
  rename(IID = te_id, ASD = ASD_dx)

# all(colnames(abcd.all) == colnames(rpoe.all))


################################################################################
################################################################################
################################################################################
################################################################################
m1m2 <- read_rds("../../RPOE/shared_data/data/m1m2.rds")
rpoe.2 <- inner_join(rpoe.cog, rpoe.demo) %>%
  inner_join(rpoe.mh) %>%
  inner_join(m1m2)
################################################################################
################################################################################
################################################################################
# train on behavior to predict ASD PGS
t1 <- abcd.all %>%
  select(ASD, starts_with(c("syn", "dsm"))) %>% ungroup()

glm.model.asd <- glm(ASD ~ ., data = t1 %>% mutate(ASD = as.numeric(ASD)))


# apply on RPOE
predictions.rpoe.glm <- predict(glm.model.asd, 
                                rpoe.2 %>% select(colnames(t1)[-1]))
p1 <- rpoe.2 %>%
  mutate(pred = predictions.rpoe.glm) %>%
  ggplot(aes(x=ASD_dx, y = pred, fill = ASD_dx)) +
  geom_violin(show.legend = F) +
  geom_boxplot(width = 0.2, fill = "white") +
  ggpubr::stat_compare_means() +
  scale_fill_manual(values = boxplot.colors) +
  labs(y = "predcited ASD") +
  bw.theme


################################################################################
################################################################################
################################################################################
# train on cognition to predict cognition total
t2 <- abcd.all %>%
  select(CP, colnames(abcd.cog2)[-c(1,2)]) %>% ungroup()

glm.model.cp <- glm(CP ~ ., data = t2)
# apply on RPOE
predictions.rpoe.glm.cp <- predict(glm.model.cp,
                                   rpoe.2 %>% select(colnames(t2)[-1]))
p2 <- rpoe.2 %>%
  mutate(pred = predictions.rpoe.glm.cp) %>%
  ggplot(aes(x=FSIQ, y = pred)) +
  geom_point(shape = 1) +
  geom_smooth(method = "lm", color = six.colors[3]) +
  ggpubr::stat_cor() +
  labs(y = "predicted CP") +
  bw.theme
################################################################################
################################################################################
# plot predictions on bothe axes, and color based on category
groups <- c("2e", "ASD, not gifted", "gifted, no dx", "no dx, not gifted")
# groups <- c("2e", "either dx, not gifted", "gifted, no dx", "no dx, not gifted")

t3 <- rpoe.2 %>%
  mutate(pred_ASD = predictions.rpoe.glm,
         pred_CP = predictions.rpoe.glm.cp) %>%
  mutate(category = 
           case_when(ASD_dx == T & FSIQ >= 120 ~ "2e",
                     ASD_dx == T & FSIQ < 120 ~ "ASD, not gifted",
                     ASD_dx == F & FSIQ >= 120 ~ "gifted, no dx",
                     ASD_dx == F & FSIQ < 120 ~ "no dx, not gifted"
                     
                     # c(ASD_dx == T | ADHD_dx == T) & FSIQ >= 120 ~ groups[1],
                     # c(ASD_dx == T | ADHD_dx == T) & FSIQ < 120 ~ groups[2],
                     # !c(ASD_dx == T | ADHD_dx == T) & FSIQ >= 120 ~ groups[3],
                     # !c(ASD_dx == T | ADHD_dx == T) & FSIQ < 120 ~ groups[4]
                     )) %>%
  select(te_id, category, starts_with("pred"), ASD_dx, FSIQ) %>%
  mutate_at(.vars = vars(starts_with("pred")), function(x) scale(x, T, T)[,1]) %>%
  drop_na(pred_ASD, pred_CP) %>%
  mutate(pred_ASD2 = ifelse(pred_ASD >= median(pred_ASD), "High", "Low"),
         pred_CP2 = ifelse(pred_CP >= median(pred_CP), "High", "Low"),
         quadrant = case_when(pred_ASD2 == "High" & pred_CP2 == "High" ~ "High ASD, High CP",
                              pred_ASD2 == "High" & pred_CP2 == "Low" ~ "High ASD, Low CP",
                              pred_ASD2 == "Low" & pred_CP2 == "High" ~ "Low ASD, High CP",
                              pred_ASD2 == "Low" & pred_CP2 == "Low" ~ "Low ASD, Low CP"))
  
p3 <- t3 %>%
  ggplot(aes(x=pred_ASD, y = pred_CP, 
             # color = ASD_dx
             # color = FSIQ > 120
             color = category
             )) +
  geom_point(size = 2, show.legend = F) +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  scale_color_manual(values = six.colors) +
  bw.theme


p4 <- t3 %>%
  ggplot(aes(x = quadrant, fill = category)) +
  geom_bar() +
  annotate("text", x = 2, y = 20, 
           label = paste("Fisher's Test; ", 
                         "p-value:", 
                         round(fisher.test(table(t3$category, t3$quadrant))$p.value, 4))) +
  scale_fill_manual(values = six.colors, name = "") +
  bw.theme +
  labs(y = "count") +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(nrow = 2))


p5 <- t3 %>%
  pivot_longer(cols = c(pred_ASD, pred_CP)) %>%
  ggplot(aes(x = category, fill = category, y = value)) +
  geom_violin(show.legend = F) + geom_boxplot(fill = "white", width = 0.2) +
  ggpubr::stat_compare_means(label = "p.signif", hide.ns = T, 
                             vjust = 0.6,
                             comparisons = combn(groups, 2, simplify = F)) +
  facet_wrap(~name) +
  labs(y = "predicted value",
       caption = paste0("n(", groups[1], "): ", sum(t3$category == groups[1]), "\n",
                        "n(", groups[2], "): ", sum(t3$category == groups[2]), "\n",
                        "n(", groups[3], "): ", sum(t3$category == groups[3]), "\n",
                        "n(", groups[4], "): ", sum(t3$category == groups[4]), "\n")) +
  scale_fill_manual(values = six.colors) +
  bw.theme +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))

patchwork::wrap_plots(patchwork::wrap_plots(p1,p2, p3), 
                      patchwork::wrap_plots(p4, p5, nrow = 1, widths = c(1,2)), ncol = 1)
ggsave("figs/models-predictions-on-RPOE_V2.png", bg = "white",
       width = 12, height = 10, units = "in", dpi = 360)

################################################################################
################################################################################
# check what MRI metrics are correlated with that ASD continuous prediction

# DTI
t3 %>%
  inner_join(rpoe.d.mri.fa3) %>%
  inner_join(rpoe.d.mri.md3) %>%
  inner_join(rpoe.demo) %>%
  ggplot(aes(x=MRI_age)) +
  geom_histogram()
t3 %>%
  inner_join(rpoe.d.mri.fa3) %>%
  inner_join(rpoe.d.mri.md3) %>%
  # inner_join(rpoe.demo) %>% filter(MRI_age < 600) %>%
  pivot_longer(cols = c(colnames(rpoe.d.mri.fa3)[-c(1:3)], colnames(rpoe.d.mri.md3)[-c(1:3)])) %>%
  filter(!grepl("SuperiorLongi|Uncinate", name)) %>% # those are not significant
  mutate(name = sub("ProjectionBasalGanglia_", "", name),
         name = sub("Association_", "", name),
         name = sub("Commissure_", "", name),
         metric = case_when(grepl("dti_fa*",name) ~ "FA",
                            grepl("md_*",name) ~ "MD"),
         name = sub(".*__", "", name)) %>%
  drop_na(pred_ASD, value) %>%
  ggplot(aes(x=pred_ASD, y = value)) +
  geom_point(shape = 1) +
  geom_smooth(method = "lm", color = six.colors[3]) +
  ggpubr::stat_cor(aes(color = ifelse(..r.. > 0, "red", "blue"),
                       label = ifelse(..p.value.. < 0.1, label, "")),
                   show.legend = F, na.rm = T) +
  scale_color_manual(values = redblu.col[c(2,1)]) +
  ggh4x::facet_grid2(rows = vars(name), cols = vars(metric), 
                     scales = "free") +
  labs(x = "predicted ASD", y = "DTI value") +
  bw.theme +
  theme(strip.text.y.right = element_text(angle = 0))
ggsave("figs/predicted-ASD_RPOE_DTI.png", bg = "white",
       width = 7, height = 10, units = "in", dpi = 360)

# resting-state
# no correlations with any of the networks ReHo or fALFF
# t3 %>%
#   left_join(rpoe.f.mri.falff2) %>%
#   left_join(rpoe.f.mri.reho2) %>%
#   pivot_longer(cols = c(colnames(rpoe.f.mri.falff2)[-c(1)], colnames(rpoe.f.mri.reho2)[-c(1)])) %>%
#   mutate(metric = sub("__.*", "", name),
#          name = sub(".*__", "", name)) %>%
#   drop_na(pred_ASD, value) %>%
#   ggplot(aes(x=pred_ASD, y = value)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "lm", color = six.colors[3]) +
#   ggpubr::stat_cor() +
#   ggh4x::facet_grid2(rows = vars(name), cols = vars(metric), 
#                      scales = "free") +
#   labs(x = "predicted ASD", y = "DTI value") +
#   bw.theme +
#   theme(strip.text.y.right = element_text(angle = 0))


################################################################################
################################################################################
# continuous CP

# DTI
# NOTHING
# t3 %>%
#   left_join(rpoe.d.mri.fa3) %>%
#   left_join(rpoe.d.mri.md3) %>%
#   pivot_longer(cols = c(colnames(rpoe.d.mri.fa3)[-c(1:3)], colnames(rpoe.d.mri.md3)[-c(1:3)])) %>%
#   # filter(!grepl("SuperiorLongi|Uncinate", name)) %>% # those are not significant
#   mutate(name = sub("ProjectionBasalGanglia_", "", name),
#          name = sub("Association_", "", name),
#          name = sub("Commissure_", "", name),
#          metric = case_when(grepl("dti_fa*",name) ~ "FA",
#                             grepl("md_*",name) ~ "MD"),
#          name = sub(".*__", "", name)) %>%
#   drop_na(pred_CP, value) %>%
#   ggplot(aes(x=pred_CP, y = value)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "lm", color = six.colors[3]) +
#   ggpubr::stat_cor(aes(color = ifelse(..r.. > 0, "red", "blue"),
#                        label = ifelse(..p.value.. < 0.1, label, "")),
#                    show.legend = F, na.rm = T) +
#   scale_color_manual(values = redblu.col[c(2,1)]) +
#   ggh4x::facet_grid2(rows = vars(name), cols = vars(metric), 
#                      scales = "free") +
#   labs(x = "predicted CP", y = "DTI value") +
#   bw.theme +
#   theme(strip.text.y.right = element_text(angle = 0))
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# train and test on ABCD
rf.model1 <- randomForest::randomForest(ASD ~ .,
                                        data = abcd.all %>%
                                          select(-IID, -ASD_dx))
glm.model1 <- glm(ASD ~ ., data = abcd.all %>%
                    select(-IID, -ASD_dx))
rf.model2 <- randomForest::randomForest(ASD_dx ~ .,
                                        data = abcd.all %>%
                                          select(-IID, -ASD))
glm.model2 <- glm(ASD_dx ~ ., data = abcd.all %>%
                    select(-IID, -ASD) %>%
                    mutate(ASD_dx = as.numeric(ASD_dx)))
##
library(tensorflow)
library(doMC)
use_condaenv("EDL")
library(keras)
rm(k.model1); rm(k.model);gc()
k_clear_session()
set.seed(123)

# Extract features and labels for training and testing
x_train <- abcd.all %>% select(-ASD, -IID, -ASD_dx) %>% as.matrix()
y_train <- abcd.all$ASD

# Normalize 
x_train <- scale(x_train)

# build model
ipt <- layer_input(shape=c(dim(x_train)[2]))
k.model1 <- ipt %>%
  layer_dense(units = 64, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 32, activation = 'relu') %>%
  layer_dropout(rate = 0.3)
out <- k.model1 %>%
  layer_dense(units = 1)
k.model <- keras_model(ipt, out)
# Compile
k.model %>% compile(
  optimizer = optimizer_adam(1e-3),
  loss = "mse",
  weighted_metrics = "mae")
# Train
history <- k.model %>% 
  fit(x = x_train, 
      y = y_train,
      epochs = 20,
      batch_size = 32,
      validation_split = 0.2,
      verbose = 1)
# Make predictions
predictions <- k.model %>% 
  predict(rpoe.all %>% select(-c(IID, ASD)) %>% as.matrix())
write_csv(predictions %>% as.data.frame(), "data/derivatives/keras-predictions-rpoe.csv")
k.predictions <- read_csv("data/derivatives/keras-predictions-rpoe.csv")
##

rf.predictions1 <- predict(rf.model1, 
                           rpoe.all %>% select(-c(IID, ASD)))
g.predictions1 <- predict(glm.model1, 
                          rpoe.all %>% select(-c(IID, ASD)))
rf.predictions2 <- predict(rf.model2, 
                           rpoe.all %>% select(-c(IID, ASD)))
g.predictions2 <- predict(glm.model2, 
                          rpoe.all %>% select(-c(IID, ASD)))
summary(glm(rpoe.all$ASD ~ predictions))
rpoe.all %>%
  mutate(DL_pred = k.predictions$V1,
         glm_pred1 = g.predictions1,
         RF_pred1 = rf.predictions1,
         glm_pred2 = g.predictions2,
         # RF_pred2 = rf.predictions2
         ) %>%
  mutate_at(.vars = vars(contains("_pred")), .funs = function(x) as.numeric(x)) %>%
  pivot_longer(cols = contains("_pred")) %>%
  ggplot(aes(x = ASD, y = value, fill = ASD)) +
  geom_violin(show.legend = F) +
  geom_boxplot(width= 0.2, fill = "white") +
  ggpubr::stat_compare_means(color = "red") +
  facet_wrap(~name, scales = "free") +
  bw.theme
ggsave("figs/models-predictions-on-RPOE.png", bg = "white",
       width = 6, height = 8, units = "in", dpi = 360)



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# all0 <- inner_join(cbcl, pgs %>% select(src_subject_id = IID, 
#                                         PGS_autism = PGS_std_PCresid)) %>%
#   left_join(cog %>% select(src_subject_id, eventname, 
#                            starts_with("nih")&contains("agecorrected"))) %>%
#   left_join(demo %>% rename(src_subject_id = IID))
# all1 <- all0 %>%
#   mutate_at(.vars = vars(contains(c("syn", "dsm", "pattern_"))), .funs = function(x) {
#     df <- all0 %>% mutate(y =x) %>%
#       select(y, interview_age, sex)
#     z_from_lm(y = df$y, x = df[,-1])
#   })
#   
# all1 %>%
#   drop_na(nihtbx_pattern_agecorrected) %>%
#   mutate(quant_ps = case_when(nihtbx_pattern_agecorrected <= quantile(nihtbx_pattern_agecorrected)[2][[1]] ~ "25%",
#                               nihtbx_pattern_agecorrected <= quantile(nihtbx_pattern_agecorrected)[3][[1]] ~ "50%",
#                               nihtbx_pattern_agecorrected <= quantile(nihtbx_pattern_agecorrected)[4][[1]] ~ "75%",
#                               nihtbx_pattern_agecorrected <= quantile(nihtbx_pattern_agecorrected)[5][[1]] ~ "100%")) %>%
#   # select(nihtbx_pattern_agecorrected, quant_ps)
#   pivot_longer(cols = contains(c("syn", "dsm5"))) %>%
#   ggplot(aes(x = PGS_autism, y = value, color = quant_ps)) +
#   geom_point(shape = 1, alpha = 0.4) +
#   geom_smooth(method = "lm") +
#   ggpubr::stat_cor() +
#   facet_wrap(~name, scales = "free")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
