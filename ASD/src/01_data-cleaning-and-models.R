################################################################################
#                       variables predicting ASD status/PGS                    #
################################################################################
rm(list = ls()); gc()
device <- ifelse(grepl("/LSS/", system("cd &pwd", intern = T)), "IDAS", "argon")
source(paste0(ifelse(device == "IDAS", "~/LSS", "/Dedicated"),
              "/jmichaelson-wdata/msmuhammad/workbench/customized-functions/correct_path.R"))
source(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R"))
################################################################################
################################################################################
project.dir <- correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/ABCD/ASD")
setwd(project.dir)
################################################################################
################################################################################
################################################################################
# load data and reformat per category for both ABCD and RPOE

########
# demo #
########
abcd.demo <- read_csv("../../../data/ABCD/abcd5/age-sex-by-eventname.csv") %>%
  filter(grepl("baseline", eventname))
abcd.dx.comp <- read_csv(correct_path("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core/abcd-general/abcd_p_screen.csv")) %>%
  select(IID = src_subject_id, common_dx = scrn_commondx, psych_dx_other = scrn_psychdx_other,ID = scrn_intdisab,
         scz_dx = scrn_schiz, ASD_dx = scrn_asd, neur_other = scrn_medcond_other) %>%
  mutate(any_psych = ifelse(common_dx == 1 |psych_dx_other == 1 |ID == 1 |scz_dx == 1 |ASD_dx == 1 |neur_other == 1,1, 0))
abcd.demo <- inner_join(abcd.demo, abcd.dx.comp);rm(abcd.dx.comp)
rpoe.demo <- read_csv("../../RPOE/shared_data/data/demo-full.csv")

#######
# cog #
#######
abcd.cog.raw <- read_csv(correct_path("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core/neurocognition/nc_y_nihtb.csv")) %>%
  filter(grepl("baseline", eventname)) %>% select(IID = src_subject_id, eventname, ends_with("agecorrected")) %>%
  inner_join(abcd.demo)

rpoe.cog.raw <- read_rds("../../RPOE/shared_data/data/m1m2-sex-corrected.rds")
## renaming columns to match naming from ABCD
rpoe.cog <- rpoe.cog.raw %>%
  rename(nihtbx_picvocab_agecorrected = picture_vocabulary_age_corrected_standard_score,
         nihtbx_flanker_agecorrected = flanker_inhibitory_control_age_corrected_standard_score,
         nihtbx_list_agecorrected = list_sorting_wm_age_corrected_standard_score,
         nihtbx_cardsort_agecorrected = dimensional_change_card_sort_age_corrected_standard_score,
         nihtbx_pattern_agecorrected = pattern_comparison_PS_age_corrected_standard_score,
         nihtbx_picture_agecorrected = picture_sequence_memory_test_age_corrected_standard_score,
         nihtbx_reading_agecorrected = oral_reading_recognition_age_corrected_standard_score)

## correct ABCD cog measures for sex
abcd.cog <- abcd.cog.raw %>%
  mutate_at(.vars = colnames(abcd.cog.raw)[-c(1,2,13,14)],
            .funs = function(x){
              df <- abcd.cog.raw %>% mutate(y = x) %>%select(y, sex)
              z_from_lm(y = df$y, x = df[,-1])
            }) %>%
  select(IID, eventname, colnames(rpoe.cog)[grepl("nihtbx", colnames(rpoe.cog))]) %>%
  drop_na()
# rm(abcd.cog.raw);gc()

#######
# MRI #
#######


#######  structural
# for ABCD, get the total volume
abcd.s.mri.raw <- read_csv(correct_path("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core/imaging/mri_y_smr_vol_aseg.csv")) %>% 
  filter(grepl("baseline", eventname)) %>%
  rowwise() %>% mutate(total_volume = sum(c_across(-c(1, 2)), na.rm = TRUE)) %>% ungroup()
rpoe.s.mri.raw <- read_rds(paste0(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/mri/data"),
                              "/derivatives/brain-vols-synthreg-age-sex-tot-vol-corrected.rds"))
## renaming columns to match naming from ABCD
colnames(rpoe.s.mri.raw) <- c("te_id", "smri_vol_scs_cbwmatterlh", colnames(rpoe.s.mri.raw)[3], "smri_vol_scs_ltventriclelh","smri_vol_scs_inflatventlh","smri_vol_scs_crbwmatterlh",
                          "smri_vol_scs_crbcortexlh","smri_vol_scs_tplh","smri_vol_scs_caudatelh", "smri_vol_scs_putamenlh","smri_vol_scs_pallidumlh","smri_vol_scs_3rdventricle", 
                          "smri_vol_scs_4thventricle","smri_vol_scs_bstem","smri_vol_scs_hpuslh", "smri_vol_scs_amygdalalh","smri_vol_scs_csf",
                          "smri_vol_scs_aal","smri_vol_scs_vedclh",
                          "smri_vol_scs_cbwmatterrh", colnames(rpoe.s.mri.raw)[21], "smri_vol_scs_ltventriclerh","smri_vol_scs_inflatventrh","smri_vol_scs_crbwmatterrh",
                          "smri_vol_scs_crbcortexrh","smri_vol_scs_tprh","smri_vol_scs_caudaterh", "smri_vol_scs_putamenrh","smri_vol_scs_pallidumrh","smri_vol_scs_hpusrh",
                          "smri_vol_scs_amygdalarh","smri_vol_scs_aar","smri_vol_scs_vedcrh",
                          colnames(rpoe.s.mri.raw)[c(34:40)])
rpoe.s.mri <- rpoe.s.mri.raw %>% select(te_id, any_of(colnames(abcd.s.mri.raw)))

## divide the volume in ABCD by total brain volume
abcd.s.mri2 <- abcd.s.mri.raw %>%
  mutate_at(.vars = -c(1,2, ncol(.data)), .funs = function(x) x/abcd.s.mri.raw$total_volume) %>%
  select(IID = src_subject_id, any_of(colnames(rpoe.s.mri))) %>% inner_join(abcd.demo)
## correct ratios for age and sex
abcd.s.mri <- abcd.s.mri2 %>%
  mutate_at(.vars = vars(starts_with("smri")),
            .funs = function(x){
              df <- abcd.s.mri2 %>% mutate(y = x) %>% select(y, interview_age, sex)
              z_from_lm(y = df$y, x = df[,-1])
            })
rm(abcd.s.mri2)
# rm(rpoe.s.mri);rm(abcd.s.mri);rm(abcd.s.mri2)


#######  diffusion
#######    FA
# renamed columns to match naming from RPOE
abcd.d.mri.fa.raw <- read_csv(correct_path("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core/imaging/mri_y_dti_fa_fs_at.csv")) %>% 
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
# correct for age and sex
abcd.d.mri.fa2 <- abcd.d.mri.fa.raw %>%
  mutate_at(.vars = vars(starts_with(c("Projection", "Comm", "Asso"))),
            .funs = function(x){
              df <- abcd.d.mri.fa.raw %>% mutate(y = x) %>% select(y, interview_age, sex)
              z_from_lm(y = df$y, x = df[,-1])
            })
# get RPOE diffusion data, keep regions in ABCD
rpoe.d.mri.raw <- read_rds(paste0(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/dwi/data"),
                              "/derivatives/recognized-tracts-summary-stats-DSI_RAW.rds"))
rpoe.d.mri.fa.raw <- rpoe.d.mri.raw %>%
  filter(roi %in% colnames(abcd.d.mri.fa2),feature %in% c("dti_fa")) %>%
  mutate(new_name = paste0(feature, "__", roi)) %>%
  pivot_wider(names_from = new_name, values_from = value, id_cols = c("te_id")) %>%
  # select(te_id, any_of(paste0(c("dti_fa__"), colnames(abcd.d.mri.fa2)))) %>%
  left_join(rpoe.demo %>% select(te_id, MRI_age, sex)) %>% drop_na(MRI_age) 
fa.keep <- names(colSums(rpoe.d.mri.fa.raw %>% select(starts_with("dti")))[!is.na(colSums(rpoe.d.mri.fa.raw %>% select(starts_with("dti"))))])
rpoe.d.mri.fa2 <- rpoe.d.mri.fa.raw %>% select(te_id, MRI_age, sex, fa.keep)
# checking if diffusion metrics are correlated with age
rpoe.d.mri.fa2 %>%
  ggplot(aes(x = dti_fa__Commissure_CorpusCallosum, y = MRI_age)) +
  geom_point(shape = 1) + geom_smooth(method = "lm") +
  ggpubr::stat_cor(method = "spearman") + bw.theme
# correct RPOE diffusion for age and sex
rpoe.d.mri.fa <- rpoe.d.mri.fa2 %>%
  mutate_at(.vars = fa.keep,
            .funs = function(x) {
              df <- rpoe.d.mri.fa2 %>% mutate(y = x) %>% select(y, MRI_age, sex)
              z_from_lm(y = df$y, x = df[,-1])
            })
rm(rpoe.d.mri.fa2)
# keep regions in RPOE
abcd.d.mri.fa <- abcd.d.mri.fa2 %>% 
  filter(grepl("baseline", eventname)) %>%
  rename_at(.vars = vars(starts_with(c("Projection", "Comm", "Asso"))),  .funs = function(x) paste0("dti_fa__", x)) %>%
  select(IID, any_of(fa.keep))
rm(abcd.d.mri.fa2);rm(fa.keep);gc()

#######    MD
### DTI
# read, rename regions to match RPOE naming
abcd.d.mri.md.raw <- read_csv(correct_path("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core/imaging/mri_y_dti_md_fs_at.csv")) %>% 
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
# correct for age and sex
abcd.d.mri.md2 <- abcd.d.mri.md.raw %>%
  mutate_at(.vars = vars(starts_with(c("Projection", "Comm", "Asso"))),
            .funs = function(x){
              df <- abcd.d.mri.md.raw %>% mutate(y = x) %>% select(y, interview_age, sex)
              z_from_lm(y = df$y, x = df[,-1])
            })
# filter RPOE diffusion data to match ABCD
rpoe.d.mri.md.raw <- rpoe.d.mri.raw %>%
  filter(roi %in% colnames(abcd.d.mri.md2), feature %in% c("md")) %>%
  mutate(new_name = paste0(feature, "__", roi)) %>%
  pivot_wider(names_from = new_name, values_from = value, id_cols = c("te_id")) %>%
  select(te_id, any_of(paste0(c("md__"), colnames(abcd.d.mri.md2)))) %>%
  left_join(rpoe.demo %>% select(te_id, MRI_age, sex)) %>% drop_na(MRI_age) 
md.keep <- names(colSums(rpoe.d.mri.md.raw %>% select(starts_with("md")))[!is.na(colSums(rpoe.d.mri.md.raw %>% select(starts_with("md"))))])
rpoe.d.mri.md2 <- rpoe.d.mri.md.raw %>% select(te_id, MRI_age, sex, md.keep)
# check if age has an effect
rpoe.d.mri.md2 %>%
  ggplot(aes(x = md__Commissure_CorpusCallosum, y = MRI_age)) +
  geom_point(shape = 1) + geom_smooth(method = "lm") +
  ggpubr::stat_cor(method = "spearman") + bw.theme
# correct for age and sex
rpoe.d.mri.md <- rpoe.d.mri.md2 %>%
  mutate_at(.vars = md.keep,
            .funs = function(x) {
              df <- rpoe.d.mri.md2 %>% mutate(y = x) %>% select(y, MRI_age, sex)
              z_from_lm(y = df$y, x = df[,-1])
            })
rm(rpoe.d.mri.md2)
# rename, filter ABCD to match
abcd.d.mri.md <- abcd.d.mri.md2 %>% 
  filter(grepl("baseline", eventname)) %>%
  rename_at(.vars = vars(starts_with(c("Projection", "Comm", "Asso"))), .funs = function(x) paste0("md__", x)) %>%
  select(IID, any_of(md.keep))
rm(abcd.d.mri.md2);rm(md.keep);gc()
# rm(abcd.d.mri.fa);rm(abcd.d.mri.fa2);rm(abcd.d.mri.md);rm(abcd.d.mri.md2);rm(rpoe.d.mri)
# rm(rpoe.d.mri.fa2);rm(rpoe.d.mri.md2);rm(rpoe.d.mri.fa);rm(rpoe.d.mri.md)


#######  fMRI
#######    fALFF
rpoe.f.mri.falff <- read_rds("../../RPOE/mri/data/derivatives/R-func/REST-fALFF-summarized-age-sex-corrected.rds")

#######    rest ReHo
rpoe.f.mri.reho <- read_rds("../../RPOE/mri/data/derivatives/R-func/REST-ReHo-summarized-age-sex-corrected.rds")

#################
# mental health #
#################
abcd.mh.raw <- read_csv(correct_path("/Dedicated/jmichaelson-sdata/ABCD/abcd_release_5_0/core/mental-health/mh_p_cbcl.csv")) %>%
  mutate(critical_1_r = cbcl_q06_p + cbcl_q15_p + cbcl_q18_p + cbcl_q40_p + cbcl_q57_p + cbcl_q59_p + cbcl_q67_p + cbcl_q70_p + cbcl_q72_p + cbcl_q91_p + cbcl_q105_p + cbcl_q107_p, 
         critical_2_r = cbcl_q18_p + cbcl_q91_p) %>%
  select(IID = src_subject_id,2, contains("syn"), contains("dsm5"), contains("critical"))
abcd.mh2 <- abcd.mh.raw[,c(T,T,grepl("_r$", colnames(abcd.mh.raw)[-c(1:2)]))] %>%
  filter(grepl("baseline", eventname)) %>% drop_na() %>%
  inner_join(abcd.demo) %>% rename_all(.funs = function(x) sub("cbcl_scr_", "", x))
abcd.mh3 <- abcd.mh2 %>% mutate_at(.vars = vars(starts_with(c("syn", "dsm5", "critical"))),
                                   .funs = function(x){
                                     df <- abcd.mh2 %>% mutate(y = x) %>% select(y, interview_age, sex)
                                     z_from_lm(y = df$y, x = df[,-1])
                                   })
# get RPOE
rpoe.mh.raw <- read_csv(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/projects/RPOE/behavior/data/derivatives/all-behavior-data-age-sex-corrected.csv"))  %>%
  select(devGenes_id, syn_anxdep_r = syn_anxious, syn_withdep_r = syn_withdrawn, syn_somatic_r = syn_somatic, 
         syn_social_r = syn_social, syn_thought_r = syn_thought, syn_attention_r = syn_attention,
         syn_rulebreak_r = syn_rulebreaking, syn_aggressive_r = syn_aggressive, dsm5_depress_r = dsm5_depressive,
         dsm5_anxdisord_r = dsm5_anxiety, dsm5_somaticpr_r = dsm5_somatic, dsm5_adhd_r = dsm5_ADHD, critical)
rpoe.mh <- rpoe.mh.raw %>% select(-critical)
abcd.mh <- abcd.mh3 %>% 
  filter(grepl("baseline", eventname)) %>%
  select(IID, any_of(colnames(rpoe.mh)))
rm(abcd.mh2);rm(abcd.mh3);gc()


#######
# PGS #
#######
abcd.pgs <- read_tsv(correct_path("/Dedicated/jmichaelson-wdata/msmuhammad/data/TT/merged-pgs.tsv")) %>%
  filter(grepl("autism|cog_performance-SSGAC-2018", PGS_name), IID %in% abcd.demo$IID) %>%
  pivot_wider(names_from = PGS_name, values_from = PGS_std_PCresid, id_cols = "IID") %>%
  select(IID, ASD = `autism-PGC-2019`, CP = `cog_performance-SSGAC-2018`)

abcd.pgs %>%
  inner_join(abcd.demo) %>%
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

abcd.all <- inner_join(abcd.cog, abcd.d.mri.fa) %>%
  inner_join(abcd.d.mri.md) %>%
  inner_join(abcd.s.mri) %>%
  inner_join(abcd.mh) %>%
  inner_join(abcd.pgs) %>%
  mutate(ASD_dx = as.factor(ASD_dx)) %>%
  select(-c(interview_age, sex, eventname)) %>%
  drop_na(ASD_dx)

rpoe.all <- full_join(rpoe.cog %>% select(-devGenes_id), rpoe.d.mri.fa) %>%
  full_join(rpoe.d.mri.md) %>%
  full_join(rpoe.s.mri) %>%
  full_join(rpoe.mh %>% inner_join(rpoe.demo %>% select(devGenes_id,te_id))%>%select(-devGenes_id)) %>%
  full_join(rpoe.f.mri.falff %>% select(-devGenes_id)) %>%
  full_join(rpoe.f.mri.reho %>% select(-devGenes_id)) %>%
  full_join(rpoe.demo %>% select(te_id, ASD_dx)) %>%
  select(-c(sex, MRI_age))

# all(colnames(abcd.all) == colnames(rpoe.all))


################################################################################
################################################################################
################################################################################
################################################################################
# m1m2 <- read_rds("../../RPOE/shared_data/data/m1m2.rds")
# rpoe.2 <- inner_join(rpoe.cog, rpoe.demo) %>%
#   inner_join(rpoe.mh) %>%
#   inner_join(m1m2)
################################################################################
################################################################################
################################################################################
# train on behavior to predict ASD PGS
t1 <- abcd.all %>%
  select(ASD, starts_with(c("syn", "dsm"))) %>% ungroup()
glm.model.asd <- glm(ASD ~ ., 
                     data = t1 %>% mutate(ASD = as.numeric(ASD)))
# apply on RPOE
rpoe.all$pred_ASD <- predict(glm.model.asd, rpoe.all %>% select(colnames(t1)[-1]))
p1 <- rpoe.all %>%
  ggplot(aes(x=ASD_dx, y = pred_ASD, fill = ASD_dx)) +
  geom_violin(show.legend = F) +
  geom_boxplot(width = 0.2, fill = "white") +
  ggpubr::stat_compare_means() +
  scale_fill_manual(values = boxplot.colors) +
  labs(y = "predcited ASD") +
  bw.theme
rm(t1)

################################################################################
################################################################################
################################################################################
# train on cognition to predict cognition total
t2 <- abcd.all %>% select(CP, colnames(abcd.cog)[-c(1,2)]) %>% ungroup()
glm.model.cp <- glm(CP ~ ., data = t2)
## apply on RPOE
# you need raw FSIQ scores here for plots
rpoe.cog.raw.RAW <- read_rds("../../RPOE/shared_data/data/m1m2.rds") %>% select(te_id, RAW_FSIQ=FSIQ)

rpoe.all$pred_CP <- predict(glm.model.cp, rpoe.all %>% select(colnames(t2)[-1]))
rpoe.all <- rpoe.all %>% left_join(rpoe.cog.raw.RAW)
p2 <- rpoe.all %>% 
  ggplot(aes(x=RAW_FSIQ, y = pred_CP)) +
  geom_point(shape = 1) +
  geom_smooth(method = "lm", color = six.colors[3]) +
  ggpubr::stat_cor() +
  labs(y = "predicted CP") +
  bw.theme
rm(t2)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# save
save(list = c("rpoe.all", "abcd.all"), file = "data/derivatives/cleaned-data.rda")
# save the models
save(glm.model.asd, glm.model.cp, file = "data/derivatives/glm-models-ASD-CP.rda")
################################################################################
################################################################################
################################################################################
patchwork::wrap_plots(p1,p2)
ggsave2("figs/models-predictions-on-RPOE.png", width = 7, height =4)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# # train on cognition and behavior to predict ASD PGS
# t8 <- abcd.all %>% select(ASD, colnames(abcd.cog)[-c(1,2)], starts_with(c("syn", "dsm"))) %>% ungroup()
# 
# glm.model.asd2 <- glm(ASD ~ ., data = t8 %>% mutate(ASD = as.numeric(ASD)))
# # apply on RPOE
# predictions.rpoe.glm.asd2 <- predict(glm.model.asd2,rpoe.all %>% select(colnames(t8)[-1]))
# p8 <- rpoe.all %>%
#   mutate(pred = predictions.rpoe.glm.asd2) %>%
#   ggplot(aes(x=ASD_dx, y = pred, fill = ASD_dx)) +
#   geom_violin(show.legend = F) +
#   geom_boxplot(width = 0.2, fill = "white") +
#   ggpubr::stat_compare_means() +
#   scale_fill_manual(values = boxplot.colors) +
#   labs(y = "predcited ASD using NIH-TB") +
#   bw.theme
# rm(t8)
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
# plot predictions on bothe axes, and color based on category
groups <- c("2e", "ASD, not gifted", "gifted, no dx", "no dx, not gifted")
# groups <- c("2e", "either dx, not gifted", "gifted, no dx", "no dx, not gifted")

t3 <- rpoe.all %>%
  mutate(category = case_when(ASD_dx == T & RAW_FSIQ >= 120 ~ "2e",
                              ASD_dx == T & RAW_FSIQ < 120 ~ "ASD, not gifted",
                              ASD_dx == F & RAW_FSIQ >= 120 ~ "gifted, no dx",
                              ASD_dx == F & RAW_FSIQ < 120 ~ "no dx, not gifted")) %>%
  select(te_id, category, starts_with("pred"), ASD_dx, RAW_FSIQ) %>%
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
             color = category)) +
  geom_point(size = 2, show.legend = F) +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  scale_color_manual(values = six.colors) +
  bw.theme


p4 <- t3 %>%
  ggplot(aes(x = quadrant, fill = category)) +
  geom_bar() +
  annotate("text", x = 2, y = 24, 
           label = paste("Fisher's Test; ", 
                         "p-value:", 
                         round(fisher.test(table(t3$category, t3$quadrant), simulate.p.value = T)$p.value, 4))) +
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
ggsave2("figs/models-predictions-on-RPOE_V2.png", width = 12, height = 10)

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
# # check what MRI metrics are correlated with that ASD continuous prediction
# 
# # DTI
# t3 %>% inner_join(rpoe.d.mri.fa) %>% inner_join(rpoe.d.mri.md) %>% inner_join(rpoe.demo) %>%
#   ggplot(aes(x=MRI_age)) + geom_histogram()
# t3 %>% inner_join(rpoe.d.mri.fa) %>% inner_join(rpoe.d.mri.md) %>%
#   # inner_join(rpoe.demo) %>% filter(MRI_age < 600) %>%
#   pivot_longer(cols = c(colnames(rpoe.d.mri.fa)[-c(1:3)], colnames(rpoe.d.mri.md)[-c(1:3)])) %>%
#   filter(!grepl("SuperiorLongi|Uncinate", name)) %>% # those are not significant
#   mutate(name = sub("ProjectionBasalGanglia_", "", name), name = sub("Association_", "", name),name = sub("Commissure_", "", name), 
#          metric = case_when(grepl("dti_fa*",name) ~ "FA", grepl("md_*",name) ~ "MD"), name = sub(".*__", "", name)) %>%
#   drop_na(pred_ASD, value) %>%
#   ggplot(aes(x=pred_ASD, y = value)) +
#   geom_point(shape = 1) +geom_smooth(method = "lm", color = six.colors[3]) +
#   ggpubr::stat_cor(aes(color = ifelse(..r.. > 0, "red", "blue"),label = ifelse(..p.value.. < 0.05, label, "")),show.legend = F, na.rm = T) +
#   scale_color_manual(values = redblu.col[c(2,1)]) +
#   ggh4x::facet_grid2(rows = vars(name), cols = vars(metric), scales = "free") +
#   labs(x = "predicted ASD", y = "DTI value") +bw.theme + theme(strip.text.y.right = element_text(angle = 0))
# ggsave2("figs/predicted-ASD_RPOE_DTI.png", width = 7, height = 10)
# 
# # resting-state FALFF
# tt4 <- corr.table(x = t3%>%left_join(rpoe.f.mri.falff)%>%select(pred_ASD,pred_CP),
#            y = t3%>%left_join(rpoe.f.mri.falff)%>%select(colnames(rpoe.f.mri.falff)[-1])) %>%
#   filter(grepl("pred_",V1),!grepl("pred_",V2),pval<0.05) %>%
#   mutate(hemisphere = case_when(grepl("__LH_", V2) ~ "left",grepl("__RH_", V2) ~ "right",
#                                 grepl("__Left-", V2) ~ "left",grepl("__Right-", V2) ~ "right"),
#          struct = sub("fALFF__", "",V2),
#          struct = sub(".H_", "", struct), struct = sub("Left-", "", struct),struct = sub("Right-", "", struct))
# t3 %>%
#   left_join(rpoe.f.mri.falff) %>%
#   pivot_longer(cols = c(colnames(rpoe.f.mri.falff)[-c(1)])) %>%
#   drop_na(pred_ASD, value) %>% inner_join(tt4 %>% rename(name = V2)) %>%
#   ggplot(aes(x=pred_ASD, y = value, color = hemisphere)) +
#   geom_point(shape = 1) +geom_smooth(method = "lm") +
#   ggpubr::stat_cor(show.legend = F) + facet_wrap(~struct) +
#   scale_color_manual(values = antique.colors) +
#   labs(x = "predicted ASD", y = "fALFF value\n(corrected for age and sex)") +
#   bw.theme + theme(strip.text.y.right = element_text(angle = 0))
# ggsave2("figs/predicted-ASD_RPOE_fALFF.png", width = 12, height = 10)
# rm(tt4)
# 
# # language
# rpoe.psvc <- read_rds("../../RPOE/language/data/derivatives/summarized-metrics-psvc.rds")
# inner_join(t3, rpoe.psvc) %>% pivot_longer(cols = colnames(rpoe.psvc)[-1]) %>%
#   pivot_longer(cols = c(pred_ASD, pred_CP), names_to = "pred", values_to = "pred_val") %>%
#   ggplot(aes(x=pred_val, y = value)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "lm", color = six.colors[3]) +
#   ggpubr::stat_cor(aes(color = ifelse(..r.. > 0, "red", "blue"), label = ifelse(..p.value.. < 0.1, label, "")), show.legend = F, na.rm = T, size = 3) +
#   scale_color_manual(values = redblu.col[c(2,1)]) +
#   ggh4x::facet_grid2(rows = vars(name), cols = vars(pred), scales = "free") +
#   labs(x = "predicted value", y = "language metric value") +
#   bw.theme + theme(strip.text.y.right = element_text(angle = 0))
# ggsave2("figs/predicted-ASD-CP_RPOE_language.png", width = 6, height = 14)

# 
# # falff
# inner_join(t3, rpoe.f.mri.falff2) %>%
#   pivot_longer(cols = colnames(rpoe.f.mri.falff2)[-1]) %>%
#   pivot_longer(cols = c(pred_ASD, pred_CP), names_to = "pred", values_to = "pred_val") %>%
#   ggplot(aes(x=pred_val, y = value)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "lm", color = six.colors[3]) +
#   ggpubr::stat_cor(aes(color = ifelse(..r.. > 0, "red", "blue"),
#                        label = ifelse(..p.value.. < 0.1, label, "")),
#                    show.legend = F, na.rm = T, size = 3) +
#   scale_color_manual(values = redblu.col[c(2,1)]) +
#   ggh4x::facet_grid2(rows = vars(name), 
#                      cols = vars(pred), 
#                      scales = "free") +
#   labs(x = "predicted value", y = "fALFF value") +
#   bw.theme +
#   theme(strip.text.y.right = element_text(angle = 0))
# ggsave("figs/predicted-ASD-CP_RPOE_fALFF.png", bg = "white",
#        width = 5, height = 9, units = "in", dpi = 360)


# # reho
# inner_join(t3, rpoe.f.mri.reho2) %>%
#   pivot_longer(cols = colnames(rpoe.f.mri.reho2)[-1]) %>%
#   pivot_longer(cols = c(pred_ASD, pred_CP), names_to = "pred", values_to = "pred_val") %>%
#   ggplot(aes(x=pred_val, y = value)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "lm", color = six.colors[3]) +
#   ggpubr::stat_cor(aes(color = ifelse(..r.. > 0, "red", "blue"),
#                        label = ifelse(..p.value.. < 0.1, label, "")),
#                    show.legend = F, na.rm = T, size = 3) +
#   scale_color_manual(values = redblu.col[c(2,1)]) +
#   ggh4x::facet_grid2(rows = vars(name), 
#                      cols = vars(pred), 
#                      scales = "free") +
#   labs(x = "predicted value", y = "ReHo value") +
#   bw.theme +
#   theme(strip.text.y.right = element_text(angle = 0))
# ggsave("figs/predicted-ASD-CP_RPOE_ReHo.png", bg = "white",
#        width = 5, height = 9, units = "in", dpi = 360)
# 


################################################################################
################################################################################
# continuous CP

# DTI
# NOTHING
# t3 %>% select(-FSIQ) %>%
#   inner_join(rpoe.d.mri.fa3) %>%
#   inner_join(rpoe.d.mri.md3) %>%
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
# check how cognition and pred_ASD interact
# 
# tmp <- t3 %>% inner_join(rpoe.d.mri.fa) %>% inner_join(rpoe.d.mri.md) %>%inner_join(rpoe.cog.raw)
# summary(glm(dti_fa__Commissure_CorpusCallosum ~ pred_ASD*pred_CP, 
#             data = tmp %>% select(dti_fa__Commissure_CorpusCallosum, pred_ASD, pred_CP,colnames(rpoe.cog.raw)[c(13:16)])))
# summary(glm(FSIQ ~ pred_ASD*dti_fa__Commissure_CorpusCallosum, 
#             data = tmp %>% select(dti_fa__Commissure_CorpusCallosum, pred_ASD, ASD_dx,colnames(rpoe.cog.raw)[c(13:16)])))
# 
# 
# cat2 <- c("High ASD, High FSIQ","High ASD, Low FSIQ","Low ASD, High FSIQ","Low ASD, Low FSIQ")
# t4 <- t3 %>% 
#   inner_join(rpoe.d.mri.fa) %>%
#   inner_join(rpoe.d.mri.md) %>%
#   inner_join(rpoe.cog.raw) %>% 
#   mutate(cat2 = case_when(c(pred_ASD >= median(pred_ASD) & FSIQ >= median(FSIQ)) ~ "High ASD, High FSIQ",
#                           c(pred_ASD >= median(pred_ASD) & FSIQ < median(FSIQ)) ~ "High ASD, Low FSIQ",
#                           c(pred_ASD < median(pred_ASD) & FSIQ >= median(FSIQ)) ~ "Low ASD, High FSIQ",
#                           c(pred_ASD < median(pred_ASD) & FSIQ < median(FSIQ)) ~ "Low ASD, Low FSIQ"))
# t4 %>%
#   pivot_longer(cols = contains("Corpus")) %>%
#   mutate(name = sub("_Commissure_", "", name),
#          name = sub("dti_fa", "FA_", name),
#          name = sub("md_", "MD_", name)) %>%
#   ggplot(aes(x = cat2, y = value, fill = cat2)) +
#   geom_violin(show.legend = F) + geom_boxplot(fill = "white", width = 0.2) +
#   ggpubr::stat_compare_means(vjust = 0.35, size = 3.5,
#                              comparisons = combn(cat2, 2, simplify = F)) +
#   scale_fill_manual(values = six.colors) +
#   facet_wrap(~name, scales = "free_y", nrow = 2) +
#   labs(x= "", y = "DTI value (residualized for age and sex)",
#        caption = paste0("n(", cat2[1], "): ", sum(t4$cat2 == cat2[1]), "\n",
#                         "n(", cat2[2], "): ", sum(t4$cat2 == cat2[2]), "\n",
#                         "n(", cat2[3], "): ", sum(t4$cat2 == cat2[3]), "\n",
#                         "n(", cat2[4], "): ", sum(t4$cat2 == cat2[4]), "\n")) +
#   bw.theme +
#   theme(axis.text.x.bottom =  element_text(angle = 45, hjust = 1))
# ggsave2("figs/predicted-ASD_RPOE_DTI-2e.png", width = 5, height = 8)


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
################################################################################
################################################################################
################################################################################
################################################################################
# theoretical mapping of pred_CP and pred_ASD on one axis
# t3 %>%
#   ggplot(aes(x=pred_ASD, y = pred_CP)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "glm") +
#   ggpubr::stat_cor(method = "spearman")
# 
# pc <- princomp(t3 %>% select(pred_CP, pred_ASD) %>% as.matrix() %>% scale())
# 
# pc2 <- prcomp(t3 %>% select(pred_CP, pred_ASD) %>% as.matrix() %>% scale())
# 
# 
# pc$scores[,1]
# pp1 <- t3 %>%
#   mutate(pc1_score = pc$scores[,1]) %>%
#   ggplot(aes(x=pred_ASD, y = pc1_score)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "lm") +
#   ggpubr::stat_cor() +
#   bw.theme
# pp2 <- t3 %>%
#   mutate(pc1_score = pc$scores[,1]) %>%
#   ggplot(aes(x=pred_CP, y = pc1_score)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "lm") +
#   ggpubr::stat_cor() +
#   bw.theme
# 
# 
# pp3 <- t3 %>%
#   mutate(pc1_score = pc$scores[,1]) %>%
#   ggplot(aes(x=pred_ASD, y = pred_CP, size = abs(pc1_score))) +
#   geom_point(shape = 1, aes(color = as.factor(sign(pc1_score)))) +
#   geom_vline(xintercept = 0, color = "pink", linetype = 2) +
#   geom_hline(yintercept = 0, color = "pink", linetype = 2) +
#   scale_color_manual(values = redblu.col[c(2,1)],
#                      name = "PC-score sign") +
#   scale_size_continuous(name = "abs(PC-score)") +
#   bw.theme
# 
# patchwork::wrap_plots(pp1, pp2, pp3, ncol = 1)
# ggsave("figs/predicted-ASD-CP_PC1-score.png", bg = "white",
#        width = 7, height = 13, units = "in", dpi = 360)
# 
# # try to correlate that combined score with other MRI metrics
# 
# t3 %>%
#   mutate(pc_score = pc$scores[,1]) %>%
#   inner_join(rpoe.d.mri.fa3) %>%
#   inner_join(rpoe.d.mri.md3) %>%
#   pivot_longer(cols = c(colnames(rpoe.d.mri.fa3)[-c(1:3)], colnames(rpoe.d.mri.md3)[-c(1:3)])) %>%
#   filter(grepl("Corpus|TractR|Inferior", name)) %>% # those are not significant
#   mutate(name = sub("ProjectionBasalGanglia_", "", name),
#          name = sub("Association_", "", name),
#          name = sub("Commissure_", "", name),
#          metric = case_when(grepl("dti_fa*",name) ~ "FA",
#                             grepl("md_*",name) ~ "MD"),
#          name = sub(".*__", "", name)) %>%
#   drop_na(pc_score, value) %>%
#   ggplot(aes(x=pc_score, y = value)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "lm", color = six.colors[3]) +
#   ggpubr::stat_cor(aes(color = ifelse(..r.. > 0, "red", "blue"),
#                        label = ifelse(..p.value.. < 0.1, label, "")),
#                    show.legend = F, na.rm = T) +
#   scale_color_manual(values = redblu.col[c(2,1)]) +
#   ggh4x::facet_grid2(rows = vars(name), cols = vars(metric), 
#                      scales = "free") +
#   labs(x = "PC1 score", y = "DTI value") +
#   bw.theme +
#   theme(strip.text.y.right = element_text(angle = 0))
# ggsave("figs/PC1_RPOE_DTI.png", bg = "white",
#        width = 7, height = 8.5, units = "in", dpi = 360)

# 
# m0 <- lm(pred_CP ~ pred_ASD, data = t3)
# slope0 <- coef(m0)[2]  # This is 'm'
# intercept0 <- coef(m0)[1]  # This is 'b'
# projected_x0 <- (t3$pred_ASD + slope0 * (t3$pred_CP - intercept0)) / (1 + slope0^2)
# projected_y0 <- slope0 * projected_x0 + intercept0
# # Calculate distance from origin (0,0) to the projected point
# distances0 <- sqrt(projected_x0^2 + projected_y0^2)
# # Get the sign based on whether the original point is above or below the line
# signed_distances0 <- distances0 * sign(t3$pred_CP - (slope0 * t3$pred_ASD + intercept0))
# 
# t3 %>%
#   mutate(distance = signed_distances0) %>%
#   ggplot(aes(x = pred_CP, y = distance)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "glm") +
#   ggpubr::stat_cor(method = "spearman")
# 
# t3 %>%
#   mutate(distance = signed_distances0) %>%
#   ggplot(aes(x = pred_ASD, y = distance, color = ASD_dx)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "glm") +
#   ggpubr::stat_cor(method = "spearman")
#


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# # correlate ASD PGS with DTI in ABCD
# abcd.d.mri.fa3 %>%
#   inner_join(abcd.pgs) %>%
#   pivot_longer(cols = starts_with("dti")) %>%
#   mutate(name = sub("dti_fa__", "", sub("Association_", "",
#                                       name)),
#          name = sub("ProjectionBasalGanglia_", "", name)) %>%
#   filter(abs(value) <=6) %>%
#   ggplot(aes(x=value, y = ASD)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "lm", color = six.colors[3]) +
#   ggpubr::stat_cor(aes(color = ifelse(..r.. > 0, "red", "blue"),
#                        label = ifelse(..p.value.. < 0.1, label, "")),
#                    show.legend = F, na.rm = T) +
#   scale_color_manual(values = redblu.col[c(2,1)]) +
#   facet_wrap(~name, scales = "free") +
#   labs(x = "DTI value", y = "ASD_PGS") +
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
# correlate FA_CC with critical score
abcd.mh3 %>%
  select(IID, starts_with("critical")) %>%
  inner_join(abcd.d.mri.fa3 %>%
               select(IID, FA_CC = dti_fa__Commissure_CorpusCallosum)) %>%
  inner_join(abcd.d.mri.md3 %>%
               select(IID, MD_CC = md__Commissure_CorpusCallosum)) %>%
  pivot_longer(cols = starts_with("critical"), names_to = "critical_type", values_to = "critical_score") %>%
  group_by(critical_type) %>% mutate(critical_category = case_when(critical_score > median(critical_score) ~ "High",
                                                                   # critical_score == median(critical_score) ~ "Median",
                                                                   critical_score <= median(critical_score) ~ "Low"),
                                     critical_category = factor(critical_category, levels = c("High", "Median", "Low"))) %>% ungroup() %>%
  pivot_longer(cols = c(FA_CC, MD_CC), names_to = "DTI_metric", values_to = "DTI_score") %>%
  filter(abs(DTI_score) <=5) %>%
  ggplot(aes(x=critical_category, y = DTI_score, fill = critical_category)) +
  geom_violin(show.legend = F) +
  geom_boxplot(fill = "white", width = 0.2) +
  ggpubr::stat_compare_means(comparisons = combn(x = c("High", "Low"), 2, simplify = F)) +
  scale_fill_manual(values = six.colors) +
  ggh4x::facet_grid2(rows = vars(DTI_metric), 
                     cols = vars(critical_type), 
                     scales = "free") +
  labs(x = "critical score", y = "DTI score") +
  bw.theme +
  theme(strip.text.y.right = element_text(angle = 0))
ggsave("figs/DTI-critical.png", bg = "white",
       width = 8, height = 8, units = "in", dpi = 360)





# make an SEM to predict the critical score in ABCD
# use age, sex, DTI, dx
# SEM
library(lavaan);library(semPlot)
sem.data <- inner_join(abcd.d.mri.fa %>% rename_at(.vars = c(2:18),
                                                   .funs = function(x) paste0("dti_fa__",x)), 
                       abcd.demo) %>%
  inner_join(abcd.d.mri.md %>% rename_at(.vars = c(2:18),
                                        .funs = function(x) paste0("dti_md__",x))) %>%
  inner_join(abcd.dx.comp) %>%
  inner_join(abcd.mh2) %>%
  inner_join(abcd.cog) %>%
  drop_na() %>%
  mutate(critical_3 = ifelse(critical_2_r > 0, 1, 0))

sem.model.1 <- '
  # Latent variable for cognition
  # cognition =~ nihtbx_picvocab_agecorrected + nihtbx_flanker_agecorrected + nihtbx_list_agecorrected + nihtbx_pattern_agecorrected + nihtbx_reading_agecorrected + nihtbx_cryst_agecorrected + nihtbx_cardsort_agecorrected + nihtbx_picture_agecorrected + nihtbx_fluidcomp_agecorrected + nihtbx_totalcomp_agecorrected
  
  # Latent variable for Behavior
  Behavior =~ dsm5_depress_r + dsm5_somaticpr_r + dsm5_opposit_r + dsm5_anxdisord_r + dsm5_adhd_r + dsm5_conduct_r

  # Latent variable for DTI
  # DTI_FA =~ dti_fa__ProjectionBasalGanglia_ThalamicRadiationL + dti_fa__ProjectionBasalGanglia_ThalamicRadiationR + dti_fa__ProjectionBrainstem_CorticospinalTractL + dti_fa__ProjectionBrainstem_CorticospinalTractR + dti_fa__ProjectionBasalGanglia_FornixL + dti_fa__ProjectionBasalGanglia_FornixR + dti_fa__ProjectionBasalGanglia_CorticostriatalTractL + dti_fa__ProjectionBasalGanglia_CorticostriatalTractR + dti_fa__Commissure_CorpusCallosum + dti_fa__Association_InferiorFrontoOccipitalFasciculusL + dti_fa__Association_InferiorFrontoOccipitalFasciculusR + dti_fa__Association_InferiorLongitudinalFasciculusL + dti_fa__Association_InferiorLongitudinalFasciculusR + dti_fa__Association_SuperiorLongitudinalFasciculusL + dti_fa__Association_SuperiorLongitudinalFasciculusR + dti_fa__Association_UncinateFasciculusL + dti_fa__Association_UncinateFasciculusR
  # DTI_MD =~ dti_md__ProjectionBasalGanglia_ThalamicRadiationL + dti_md__ProjectionBasalGanglia_ThalamicRadiationR + dti_md__ProjectionBrainstem_CorticospinalTractL + dti_md__ProjectionBrainstem_CorticospinalTractR + dti_md__ProjectionBasalGanglia_FornixL + dti_md__ProjectionBasalGanglia_FornixR + dti_md__ProjectionBasalGanglia_CorticostriatalTractL + dti_md__ProjectionBasalGanglia_CorticostriatalTractR + dti_md__Commissure_CorpusCallosum + dti_md__Association_InferiorFrontoOccipitalFasciculusL + dti_md__Association_InferiorFrontoOccipitalFasciculusR + dti_md__Association_InferiorLongitudinalFasciculusL + dti_md__Association_InferiorLongitudinalFasciculusR + dti_md__Association_SuperiorLongitudinalFasciculusL + dti_md__Association_SuperiorLongitudinalFasciculusR + dti_md__Association_UncinateFasciculusL + dti_md__Association_UncinateFasciculusR
  # DTI =~ DTI_FA + DTI_MD
  DTI =~ dti_fa__Commissure_CorpusCallosum + dti_md__Commissure_CorpusCallosum
  
  # Regressions: predicting Critical Behaviors
  critical ~ DTI + Behavior + interview_age + sex + any_psych
  
  # Optional: Correlations between predictors
  # cognition ~~ DTI_FA + DTI_MD + Behavior + interview_age + sex + any_psych
  # cognition ~~ interview_age
  # cognition ~~ sex
  # cognition ~~ DTI_FA + DTI_MD
  # cognition ~~ Behavior
  # cognition ~~ any_psych
  # 
  Behavior ~~ any_psych + interview_age + sex
  # Behavior ~~ interview_age
  # Behavior ~~ sex
  # Behavior ~~ any_psych
  # Behavior ~~ DTI_FA + DTI_MD
  # 
  # DTI ~~ cognition + any_psych + Behavior + interview_age + sex
  # DTI ~~ any_psych
  # DTI ~~ interview_age
  # DTI ~~ sex
'
sem.model.2 <- '
  # Latent variable for cognition
  cognition =~ nihtbx_picvocab_agecorrected + nihtbx_flanker_agecorrected + nihtbx_list_agecorrected + nihtbx_pattern_agecorrected + 
               nihtbx_reading_agecorrected + nihtbx_cryst_agecorrected + nihtbx_cardsort_agecorrected + nihtbx_picture_agecorrected + 
               nihtbx_fluidcomp_agecorrected + nihtbx_totalcomp_agecorrected

  # Latent variable for Behavior
  Behavior =~ dsm5_depress_r + dsm5_somaticpr_r + dsm5_opposit_r + dsm5_anxdisord_r + dsm5_adhd_r + dsm5_conduct_r
  

  # Latent variable for DTI (consider using FA and MD as sub-latents or a composite score)
  DTI_FA =~ dti_fa__ProjectionBasalGanglia_ThalamicRadiationL + dti_fa__ProjectionBasalGanglia_ThalamicRadiationR + dti_fa__ProjectionBrainstem_CorticospinalTractL + dti_fa__ProjectionBrainstem_CorticospinalTractR + dti_fa__ProjectionBasalGanglia_FornixL + dti_fa__ProjectionBasalGanglia_FornixR + dti_fa__ProjectionBasalGanglia_CorticostriatalTractL + dti_fa__ProjectionBasalGanglia_CorticostriatalTractR + dti_fa__Commissure_CorpusCallosum + dti_fa__Association_InferiorFrontoOccipitalFasciculusL + dti_fa__Association_InferiorFrontoOccipitalFasciculusR + dti_fa__Association_InferiorLongitudinalFasciculusL + dti_fa__Association_InferiorLongitudinalFasciculusR + dti_fa__Association_SuperiorLongitudinalFasciculusL + dti_fa__Association_SuperiorLongitudinalFasciculusR + dti_fa__Association_UncinateFasciculusL + dti_fa__Association_UncinateFasciculusR
  DTI_MD =~ dti_md__ProjectionBasalGanglia_ThalamicRadiationL + dti_md__ProjectionBasalGanglia_ThalamicRadiationR + dti_md__ProjectionBrainstem_CorticospinalTractL + dti_md__ProjectionBrainstem_CorticospinalTractR + dti_md__ProjectionBasalGanglia_FornixL + dti_md__ProjectionBasalGanglia_FornixR + dti_md__ProjectionBasalGanglia_CorticostriatalTractL + dti_md__ProjectionBasalGanglia_CorticostriatalTractR + dti_md__Commissure_CorpusCallosum + dti_md__Association_InferiorFrontoOccipitalFasciculusL + dti_md__Association_InferiorFrontoOccipitalFasciculusR + dti_md__Association_InferiorLongitudinalFasciculusL + dti_md__Association_InferiorLongitudinalFasciculusR + dti_md__Association_SuperiorLongitudinalFasciculusL + dti_md__Association_SuperiorLongitudinalFasciculusR + dti_md__Association_UncinateFasciculusL + dti_md__Association_UncinateFasciculusR
  DTI =~ DTI_FA + DTI_MD

  # Regression paths: predicting critical behaviors with direct effects
  critical ~ cognition + DTI + Behavior + interview_age + sex + any_psych

  # Correlations between predictors (optional)
  cognition ~~ Behavior + DTI + interview_age + sex + any_psych
  DTI ~~ Behavior + cognition + interview_age + sex + any_psych
  Behavior ~~ cognition + DTI + interview_age + sex + any_psych
'
# Fit the SEM model
# sem.fit.1 <- sem(sem.model.1, data = sem.data %>%
#                    select(starts_with(c("dti_fa", "dti_md", "nihtbx", "dsm5")),
#                           interview_age, sex, any_psych, critical = critical_1_r) %>%
#                    mutate_at(.vars = vars(-c(any_psych, sex)), .funs = function(x) scale(x,T,T)[,1]))
sem.fit.2 <- sem(sem.model.2, data = sem.data %>%
                   select(starts_with(c("dti_fa", "dti_md", "nihtbx", "dsm5")),
                          interview_age, sex, any_psych, critical = critical_3) %>%
                   mutate_at(.vars = vars(-c(any_psych, sex, critical)), 
                             .funs = function(x) scale(x, TRUE, TRUE)[,1]), 
                 estimator = "DWLS")


# View summary with fit indices
summary(sem.fit.2, fit.measures = TRUE, standardized = TRUE)

sem.parameters.2 <- parameterEstimates(sem.fit.2)
semPaths(sem.fit.2, what = "std")



sem.data %>%
  select(starts_with(c("dti_fa", "dti_md", "nihtbx", "dsm5")),
         interview_age, sex, any_psych, 
         critical = critical_2_r) %>%
  mutate_at(.vars = vars(-c(any_psych, sex, critical)), 
            .funs = function(x) scale(x, TRUE, TRUE)[,1]) %>%
  mutate(critical = as.factor(critical)) %>%
  pivot_longer(cols = starts_with("dsm5")) %>%
  # filter(abs(value) < 10) %>%
  ggplot(aes(x=critical, y = value, fill = critical)) +
  geom_violin() +
  geom_boxplot(width = 0.2, fill = "white") +
  geom_hline(yintercept = 0, linetype=2, color = "pink") +
  ggpubr::stat_compare_means(comparisons = combn(c("0","1","2","3"), 
                                                 m = 2, 
                                                 simplify = F)) +
  scale_fill_manual(values = six.colors) +
  facet_wrap(~name, scales = "free") +
  bw.theme

sem.data %>%
  # mutate(critical_1_r = residuals(glm(critical_1_r ~ interview_age, 
  #                                     data = sem.data, family = poisson()))) %>%
  ggplot(aes(x=critical_1_r, y = nihtbx_pattern_agecorrected)) +
  geom_point(shape = 1) +
  geom_smooth(method = "lm", color = six.colors[3]) +
  ggpubr::stat_cor() +
  bw.theme

p1 <- sem.data %>%
  mutate(nihtbx_pattern_agecorrected = scale(nihtbx_pattern_agecorrected, T,T)[,1]) %>%
  pivot_longer(cols = c(critical_1_r, critical_2_r)) %>%
  mutate(value_2 = ifelse(value == 0, "0", ">=1"),
         value_2 = factor(value_2, levels = c("0", ">=1")),
         name = ifelse(name == "critical_1_r", "ALL_critical_questions", "Suicidal_questions")) %>%
  ggplot(aes(x=value_2, y = nihtbx_pattern_agecorrected, fill = value_2)) +
  geom_violin(show.legend = F) +
  geom_boxplot(width = 0.2, fill = "white") +
  ggpubr::stat_compare_means(label.y.npc =  0.8, 
                             comparisons = combn(x = c("0", ">=1"),
                                                 m = 2, simplify = F), 
                             vjust = 0.1, size = 3) +
  bw.theme +
  facet_wrap(~name) +
  scale_fill_manual(values = six.colors[c(2,1)]) +
  labs(x="problems_score")

p2 <- sem.data %>%
  mutate(nihtbx_pattern_agecorrected = scale(nihtbx_pattern_agecorrected, T,T)[,1]) %>%
  pivot_longer(cols = c(critical_1_r, critical_2_r)) %>%
  mutate(name = ifelse(name == "critical_1_r", "ALL_critical_questions", "Suicidal_questions")) %>%
  ggplot(aes(x=value)) +
  geom_bar() +
  geom_text(stat = "count", 
            aes(label = ..count..), 
            vjust = 0.5, size = 3) +  
  bw.theme +
  facet_wrap(~name, scales = "free")

patchwork::wrap_plots(p2,p1, ncol = 1, heights = c(1,3))
ggsave("figs/citical-vs-PS.png", bg = "white",
       width = 8, height = 9, units = "in", dpi = 360)
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
