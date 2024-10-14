##################################################################
#                   T1 MRI processing for ABCD                   #
##################################################################

# extract arguments  from the bash call 
args <- commandArgs(trailingOnly = T)
### try this 
# args <- c("/Dedicated/jmichaelson-wdata/msmuhammad/projects/ABCD/MRI",
#           "/Dedicated/jmichaelson-sdata/ABCD/abcd_release_3_0/imaging/data/ABCD_T1/nii",
#           "/Dedicated/jmichaelson-sdata/ABCD/abcd_release_3_0/imaging/data/ABCD_T1/nii/NDAR_INVZZZNB0XC.nii.gz",
#           "/Dedicated/jmichaelson-wdata/msmuhammad/projects/ABCD/MRI/derivatives")
###

MAIN_DIR <- args[1]
RAWDATA_DIR <- args[2]



subjT1 <- args[3]
OUT_DIR <- args[4]
numCores <- args[5]

subject <- sub("^.*/([^.]*).*$", "\\1", subjT1)
print(paste0("started subject: ", subject))
system(paste0("export subject=", subject))

#global command?
subject.dir <- paste0(OUT_DIR, "/", subject)



# parameters dataframe
# this is for skullstripping. higher watershed value means strip less from the skull and vice versa. 
# default value for wsthresh is 25 (alredy included in proposed parameters).
wsthreshold <- c(25)

# the other options other than "recon2" are mainly to be used if you made manual changes to intensity normalization or white matter
# if you didn't make any changes and still running the pipeline for the first time, just use "autorecon2"
# recon2 <- c("autorecon2", "autorecon2 autorecon2-wm", "autorecon2 autorecon2-cp", "autorecon2 autorecon2-pial")
recon2 <- c("autorecon2")
recon3 <- c("autorecon3")
parameters <- expand.grid(wsthresh = wsthreshold, recon2 = recon2, recon3 = recon3)

# cmd dataframe 
# This would be the same dataframe for every subject. 
# it has different columns and each column is representing the commands in order

cmd <- data.frame(row.names = seq(1, nrow(parameters)))
for (i in 1:nrow(parameters)) {
  SUBJECTS_DIR <- paste0(subject.dir, "/parameter-set-", i)
  cmd1 <- paste0("export SUBJECTS_DIR=", SUBJECTS_DIR)
  cmd$cmd1[i] <- cmd1
  logfile <- paste0(SUBJECTS_DIR, "/logfile.log")
  cmd2 <- paste0("export logfile=", logfile)
  cmd$cmd2[i] <- cmd2
  cmd3 <- paste0("rm -rf ", SUBJECTS_DIR, "; mkdir -p ", SUBJECTS_DIR, "/mri/orig")
  cmd$cmd3[i] <- cmd3
  cmd4 <- paste0("mri_convert ", subjT1, " ", SUBJECTS_DIR, "/mri/orig/", "001", ".mgz")
  cmd$cmd4[i] <- cmd4
  cmd$cmd5[i] <- paste0("recon-all -autorecon1 -s parameter-set-", i, " -sd ", 
                        subject.dir, " -skullstrip -gcut -norandomness -deface -wsthresh ",
                        parameters[i,1], " -subcortseg -no-wsgcaatlas -", parameters[i,2], 
                        " -", parameters[i,3])
  # cmd$path[i] <- paste0(subject.dir, "/parameter-set-", i)
  # make 3d model for cortical
  cmd6 <- paste0("mris_convert --combinesurfs ", SUBJECTS_DIR, "/surf/lh.pial ", SUBJECTS_DIR, "/surf/rh.pial ", SUBJECTS_DIR, "/cortical.stl")
  cmd$cmd6[i] <- cmd6
  # make 3d model for subcortical regions
  cmd$cmd7[i] <- paste0("mkdir -p ", SUBJECTS_DIR, "/subcortical")
  cmd$cmd8[i] <- paste0("mri_convert ", SUBJECTS_DIR, "/mri/aseg.mgz ", SUBJECTS_DIR, "/subcortical/subcortical.nii")
  cmd$cmd9[i] <- paste0("mri_binarize --i ", SUBJECTS_DIR, "/subcortical/subcortical.nii ",
                        "--match 2 3 24 31 41 42 63 72 77 51 52 13 12 43 50 4 11 26 58 49 10 17 18 53 54 44 5 80 14 15 30 62 ",
                        "--inv ", "--o ", SUBJECTS_DIR, "/subcortical/bin.nii")
  cmd$cmd10[i] <- paste0("fslmaths ", SUBJECTS_DIR, "/subcortical/subcortical.nii ", 
                         "-mul ", SUBJECTS_DIR, "/subcortical/bin.nii ",
                         SUBJECTS_DIR, "/subcortical/subcortical.nii.gz")
  cmd$cmd11[i] <- paste0("cp ", SUBJECTS_DIR, "/subcortical/subcortical.nii.gz ", 
                         SUBJECTS_DIR, "/subcortical/subcortical_tmp.nii.gz")
  cmd$cmd12[i] <- paste0("gunzip -f ", SUBJECTS_DIR, "/subcortical/subcortical_tmp.nii.gz")
  cmd$cmd13[i] <- paste0("for i in 7 8 16 28 46 47 60 251 252 253 254 255; do mri_pretess ", SUBJECTS_DIR, "/subcortical/subcortical_tmp.nii ",
                         "$i ", SUBJECTS_DIR, "/mri/norm.mgz ", SUBJECTS_DIR, "/subcortical/subcortical_tmp.nii ", "; done")
  cmd$cmd14[i] <- paste0("fslmaths ", SUBJECTS_DIR, "/subcortical/subcortical_tmp.nii -bin ", SUBJECTS_DIR, "/subcortical/subcortical_bin.nii")
  cmd$cmd15[i] <- paste0("mri_tessellate ", SUBJECTS_DIR, "/subcortical/subcortical_bin.nii.gz 1 ", SUBJECTS_DIR, "/subcortical/subcortical")
  cmd$cmd16[i] <- paste0("mris_convert ", SUBJECTS_DIR, "/subcortical/subcortical ", SUBJECTS_DIR, "/subcortical.stl")
  #combine  cortical and subcortical
  cmd$cmd17[i] <- "echo 'solid '$SUBJECTS_DIR'/final.stl' > $SUBJECTS_DIR/final.stl"
  cmd$cmd18[i] <- "sed '/solid vcg/d' ${SUBJECTS_DIR}/cortical.stl >> $SUBJECTS_DIR/final.stl"
  cmd$cmd19[i] <- "sed '/solid vcg/d' ${SUBJECTS_DIR}/subcortical.stl >> $SUBJECTS_DIR/final.stl"
  cmd$cmd20[i] <- "echo 'endsolid '${SUBJECTS_DIR}'/final.stl' >> $SUBJECTS_DIR/final.stl"
  
}


# submitting commands for each parameter set
# library(foreach)
# library(doParallel)
# registerDoParallel(numCores)
# foreach(i = 1:nrow(cmd)) %dopar% {
  print("###################################################################################################################")
  for (j in 1:ncol(cmd)) {
    print(paste0("start of cmd: ", j, " for parameter set: ", i))
    system(cmd[i,j])
    # print(cmd[i,j])
    print(paste0("cmd: ", j, " is completed for parameter set: ", i))
  }
  print("###################################################################################################################")
# }