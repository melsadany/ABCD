################################################################################
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
load("data/derivatives/cleaned-data.rda")
rpoe.all <- rpoe.all %>% mutate_at(.vars = vars(starts_with("pred_")), .funs = function(x) scale(x,T,T)[,1])
################################################################################
################################################################################
## ASD w Diffusion
dd <- rpoe.all %>%
  select(te_id,pred_ASD,pred_CP, starts_with(c("dti", "md"))) %>%
  pivot_longer(cols = -c(te_id,pred_ASD, pred_CP)) %>%
  mutate(metric = case_when(grepl("dti_fa",name)~"FA",grepl("md__",name)~"MD"),
         tract = sub("dti_fa__.*_", "", name),tract = sub("md__.*_", "", tract),
         hemisphere = case_when(grepl("L$", tract)~"Left",grepl("R$", tract)~"Right")) %>%
  filter(te_id != "2E_012") %>%
  drop_na(value, pred_ASD)
dd %>% 
  ggplot(aes(x = value, y = pred_ASD)) +
  geom_point(shape = 1) + geom_smooth(method = "lm",color =six.colors[3]) +
  ggpubr::stat_cor(aes(color = ifelse(..r.. > 0, "red", "blue"),label = ifelse(..p.value.. < 0.05, label, "")),
                   show.legend = F, na.rm = T,method = "spearman") +
  scale_color_manual(values = c(redblu.col[c(2,1)],"grey")) +
  ggh4x::facet_grid2(rows = vars(tract), cols = vars(metric), scales = "free") +
  labs(y = "predicted ASD", x = "DTI value", 
       caption = paste0("n(samples): ", length(unique(dd$te_id)))) +
  bw.theme + theme(strip.text.y.right = element_text(angle = 0))
ggsave2("figs/predicted-ASD_RPOE_DTI.png", width = 7, height = 10)

dd %>%
  filter(grepl("Corpus",tract),metric=="FA") %>%
  ggplot(aes(x = value, y = pred_ASD)) +
  geom_point(shape = 1) + geom_smooth(method = "lm",color =six.colors[3]) +
  ggpubr::stat_cor(color = "red", show.legend = F, na.rm = T) +
  scale_color_manual(values = c(redblu.col[c(2,1)],"grey")) +
  labs(y = "predicted ASD", x = "Commissure Corpus Callosum FA", 
       caption = paste0("n(samples): ", length(unique(dd$te_id)))) +
  bw.theme
ggsave2("figs/predicted-ASD_RPOE_DTI-CC.png", width = 5, height = 4)

################################################################################
################################################################################
## ASD w fALFF
# get signficant ROIs only
tmp.corr <- corr.table(rpoe.all%>%select(pred_ASD,pred_CP, ASD_dx),
                             rpoe.all%>%select(starts_with(c("fALFF", "ReHo")))) %>%
  filter(grepl("pred_",V1)|V1=="ASD_dx",!grepl("pred_",V2)) %>%
  mutate(hemisphere = case_when(grepl("__LH_", V2) ~ "left",grepl("__RH_", V2) ~ "right",grepl("__Left-", V2) ~ "left",grepl("__Right-", V2) ~ "right"),
         type = sub("__.*", "", V2),
         struct = sub(".*__", "",V2),struct = sub(".H_", "", struct), struct = sub("Left-", "", struct),struct = sub("Right-", "", struct))
dd <- rpoe.all %>% select(te_id,pred_ASD,pred_CP, tmp.corr$V2) %>%
  pivot_longer(cols = -c(te_id, starts_with("pred_"))) %>% 
  drop_na(value, pred_ASD)

# dd %>%
#   inner_join(tmp.corr %>% filter(type == "fALFF",pval < 0.2) %>% rename(name = V2)) %>%
#   ggplot(aes(x=pred_ASD, y = value, color = hemisphere)) +
#   geom_point(shape = 1) +geom_smooth(method = "lm") +
#   ggpubr::stat_cor(show.legend = F, label.y.npc = c(0.95)) + facet_wrap(~struct) +
#   scale_color_manual(values = antique.colors) +
#   labs(x = "predicted ASD", y = "fALFF value\n(corrected for age and sex)") +
#   bw.theme + theme(strip.text.y.right = element_text(angle = 0))
# ggsave2("figs/predicted-ASD_RPOE_fALFF.png", width = 12, height = 10)

dd %>%
  inner_join(tmp.corr %>% filter(type=="ReHo", pval < 0.2) %>% rename(name = V2)) %>%
  ggplot(aes(x=pred_ASD, y = value, color = hemisphere)) +
  geom_point(shape = 1) +geom_smooth(method = "lm") +
  ggpubr::stat_cor(show.legend = F, label.y.npc = c(0.95)) + facet_wrap(~struct) +
  scale_color_manual(values = antique.colors) +
  labs(x = "predicted ASD", y = "ReHo value\n(corrected for age and sex)",
       caption = paste0("n(samples)")) +
  bw.theme + theme(strip.text.y.right = element_text(angle = 0))
ggsave2("figs/predicted-ASD_RPOE_ReHo.png", width = 15, height = 10)

################################################################################
################################################################################
## ASD w language
rpoe.psvc <- read_rds("../../RPOE/language/data/derivatives/summarized-metrics-psvc-per-prompt.rds") %>%
  select(-c(19:23,11:13)) %>% mutate(task = ifelse(nchar(prompt)==1,3,1)) %>%
  group_by(te_id,task) %>% summarise_at(.vars = vars(-c(prompt)),
                                        .funs = function(x) mean(x, na.rm=T))
rpoe.all %>% select(te_id, pred_ASD) %>%
  inner_join(rpoe.psvc %>% filter(task==3) %>% select(-task)) %>%
  pivot_longer(cols = -c(te_id, pred_ASD)) %>%
  ggplot(aes(x = pred_ASD, y = value)) +
  geom_point() + geom_smooth(method = "lm") +
  ggpubr::stat_cor(color = "red") +
  facet_wrap(~name, scales = "free") + bw.theme
### NOTHING
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
