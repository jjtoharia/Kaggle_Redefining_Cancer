### Inicialización (setwd() y rm() y packages):

# setwd(getwd())
try(setwd('D:/Personal/Dropbox/Musica/Tmp_Kaggle/RedefiningCancer'), silent=TRUE)
try(setwd('C:/Personal/Dropbox/Musica/Tmp_Kaggle/RedefiningCancer'), silent=TRUE)
try(setwd('C:/Users/Adriana/Dropbox/Musica/Tmp_Kaggle/RedefiningCancer'), silent=TRUE)
rm(list = ls()) # Borra todos los elementos del entorno de R.

# install.packages("stringr")
library(stringr)

s_input_path <- "C:/Servidor/JJTZ-Sync/Kaggle_RedefiningCancer/"
if(!file.exists(paste0(s_input_path, "FileDescriptions.txt")))
  s_input_path <- str_replace(s_input_path, pattern = "JJTZ-Sync", replacement = "JJTZ_Sync")
if(!file.exists(paste0(s_input_path, "FileDescriptions.txt")))
  s_input_path <- "D:/Servidor/JJTZ-Sync/Kaggle_RedefiningCancer/"
if(!file.exists(paste0(s_input_path, "FileDescriptions.txt")))
  s_input_path <- str_replace(s_input_path, pattern = "JJTZ-Sync", replacement = "JJTZ_Sync")
if(!file.exists(paste0(s_input_path, "FileDescriptions.txt")))
  stop(paste0('Fichero <', paste0(s_input_path, "FileDescriptions.txt"), '> NO encontrado. Detenemos el proceso...'))

s_output_path <- "Out/"

# options(echo = FALSE) # ECHO OFF
print('###########################################')
print('# Redefining Cancer Treatment - JJTZ 2017')
# Memorial Sloan Kettering Cancer Center - Belongs to the NIPS 2017 Competition Track (https://nips.cc/Conferences/2017/CompetitionTrack)
# Redefining Cancer Treatment - Predict the effect of Genetic Variants to enable Personalized Medicine
print('###########################################')
# install.packages("data.table")
suppressPackageStartupMessages(library(data.table))
# Process in parallel:
# install.packages("doParallel")
suppressPackageStartupMessages(library(foreach))
library(iterators)
library(parallel)
library(doParallel)
# # Process in parallel: Ejemplo de uso:
# cl <- makeCluster(detectCores(), type='PSOCK') # library(doParallel) [turn parallel processing on]
# registerDoParallel(cl) # library(doParallel) [turn parallel processing on]
# registerDoSEQ() # library(doParallel) [turn parallel processing off and run sequentially again]
# #

# install.packages("bit64")
suppressPackageStartupMessages(library(bit64))
# ##################################################
# ## Funciones:
# ##################################################

source("RedefCancer_jjtz_funciones.R", encoding = "UTF-8")

# ##################################################
# ## Inicio:
# ##################################################
Proyecto <- "Redefining Cancer Treatment"
print(paste0(Sys.time(), ' - ', 'Proyecto = ', Proyecto, ' (Nodo=', Sys.info()["nodename"], ')'))

Proyecto.s <- str_replace_all(Proyecto, "\\(|\\)| |:", "_") # Quitamos espacios, paréntesis, etc.

# Inicializamos variables:
# # NOTA: Dejamos un Core de la CPU "libre" para no "quemar" la máquina:
# cl <- makeCluster(detectCores(), type='PSOCK') # library(doParallel) [turn parallel processing on]
# registerDoParallel(cl) # library(doParallel) [turn parallel processing on]

# memory.limit(size = 16000)

systime_ini <- proc.time()

# --------------------------------------------------------
G_b_DEBUG <- FALSE # Reducimos todo para hacer pruebas más rápido
NUM_BLOQUES <- 3
b_hay_cambios <- TRUE # TRUE para forzar todo
# --------------------------------------------------------
# Cargamos los Levels comunes en Gene y Variation (en trainset y en testset):
# --------------------------------------------------------
if(!mi_load(s_input_path, 'Gene_Variation_Factors.RData'))
{
  jjprint('Creando los Levels comunes en Gene y Variation [Gene_Variation_Factors.RData]...')
  jjprint('Leyendo trainset original...')
  full_trainset <- fread(file.path(s_input_path, "training_variants"))
  jjprint('Leyendo testset original...')
  testset <- fread(file.path(s_input_path, "test_variants"))
  v_Genes <- unique(full_trainset$Gene)[unique(full_trainset$Gene) %in% unique(testset$Gene)]
  v_Genes <- v_Genes[v_Genes %in% names(table(full_trainset$Gene)[table(full_trainset$Gene) > 0.001 * nrow(full_trainset)])]
  v_Genes <- v_Genes[v_Genes %in% names(table(testset$Gene)[table(testset$Gene) > 0.001 * nrow(testset)])]
  v_Genes <- c(v_Genes, '--')
  v_Variations <- unique(full_trainset$Variation)[unique(full_trainset$Variation) %in% unique(testset$Variation)]
  v_Variations <- v_Variations[v_Variations %in% names(table(full_trainset$Variation)[table(full_trainset$Variation) > 0.001 * nrow(full_trainset)])]
  v_Variations <- v_Variations[v_Variations %in% names(table(testset$Variation)[table(testset$Variation) > 0.001 * nrow(testset)])]
  v_Variations <- c(v_Variations, '--')
  save(v_Genes, v_Variations, file = file.path(s_input_path, 'Gene_Variation_Factors.RData'))
  b_hay_cambios <- TRUE
  rm(full_trainset, testset); gc()
}
# --------------------------------------------------------
# Tratamiento de la variable "Variation":
# --------------------------------------------------------
if(!mi_load(s_input_path, 'Variations.RData'))
{
  jjprint('Procesando campo Variation [Variations.RData]...')
  jjprint('Leyendo trainset original...')
  full_trainset <- fread(file.path(s_input_path, "training_variants"))
  jjprint('Leyendo testset original...')
  testset <- fread(file.path(s_input_path, "test_variants"))
  testset[,Class:=NA]
  full_dataset <- rbind(full_trainset, testset)
  full_dataset[, Variation := stringr::str_trim(Variation, "both")]
  full_dataset[, Vari_b_ampl := as.factor(ifelse(grepl("ampl", Variation, ignore.case = T), "1", "0"))]
  # full_dataset[, Vari_b_null := as.factor(ifelse(grepl("null", Variation, ignore.case = T), "1", "0"))]
  full_dataset[, Vari_b_del := as.factor(ifelse(grepl("del", Variation, ignore.case = T), "1", "0"))]
  # full_dataset[, Vari_b_dup := as.factor(ifelse(grepl("dup", Variation, ignore.case = T), "1", "0"))]
  full_dataset[, Vari_b_ins := as.factor(ifelse(grepl("ins", Variation, ignore.case = T), "1", "0"))]
  full_dataset[, Vari_b_Mut := as.factor(ifelse(grepl("Mut", Variation, ignore.case = T), "1", "0"))]
  full_dataset[, Vari_b_Fus := as.factor(ifelse(grepl("Fus", Variation, ignore.case = T), "1", "0"))]
  # full_dataset[, Vari_b_Sil := as.factor(ifelse(grepl("Sil", Variation, ignore.case = T), "1", "0"))]
  # full_dataset[, Vari_b_spli := as.factor(ifelse(grepl("splic", Variation, ignore.case = T), "1", "0"))]
  # full_dataset[, Vari_b_Prom := as.factor(ifelse(grepl("Prom", Variation, ignore.case = T), "1", "0"))]
  full_dataset[, Vari_b_Trun := as.factor(ifelse(grepl("Trun", Variation, ignore.case = T), "1", "0"))]
  full_dataset[, Vari_first := substr(Variation, 1, 1)]
  full_dataset[Vari_first %in% c("1","2","3","4","5","9","B","n","O","p","X","Z"), Vari_first := NA]
  full_dataset[is.na(Vari_first), Vari_first := "--"]
  full_dataset[, Vari_first := as.factor(Vari_first)]
  full_dataset[, Vari_medio := substr(Variation, 2, str_length(Variation)-1)]
  suppressWarnings( # "Warning: NAs introduced by coercion"
    full_dataset[, Vari_medio_num := as.numeric(Vari_medio)]
  )
  full_dataset[, Vari_medio := NULL] # Ya no lo queremos
  if(min(full_dataset[, Vari_medio_num], na.rm = TRUE) != 0) # Si no hay ceros,
    full_dataset[is.na(Vari_medio_num), Vari_medio_num := 0] # pues los NAs a cero!
  # table(is.na(full_dataset$Vari_medio_num), full_dataset$Class)
  # hist(full_dataset$Vari_medio_num)
  full_dataset[, Vari_last := substr(Variation, str_length(Variation), str_length(Variation))]
  full_dataset[Vari_last %in% c("0","1","2","3","4","5","6","7","c","e","g","k","l","m","p","t"), Vari_last := NA]
  full_dataset[is.na(Vari_last), Vari_last := "--"]
  full_dataset[, Vari_last := as.factor(Vari_last)]
  
  # # full_dataset[is.na(Class), Class := 999]
  # # table(full_dataset$Vari_last, full_dataset$Class)
  # sapply(full_dataset[,.(Vari_first, Vari_medio_num, Vari_last)], uniqueN)
  # summary(full_dataset)
  # names(sapply(full_dataset, anyNA))[sapply(full_dataset, anyNA)] # Columnas con algún NA
  
  flds <- c("ID", "Vari_medio_num")
  for(fld in colnames(full_dataset))
    if(is.factor(full_dataset[,(fld),with=F][[1]]))
      flds <- c(flds, fld)
  train_Variations <- full_dataset[!is.na(Class), flds, with=F]
  setkey(train_Variations, ID)
  test_Variations <- full_dataset[is.na(Class), flds, with=F]
  setkey(test_Variations, ID)
  
  save(train_Variations, test_Variations, file = file.path(s_input_path, 'Variations.RData'))
  b_hay_cambios <- TRUE
  rm(full_dataset, full_trainset, testset, flds); gc()
}
# --------------------------------------------------------
# Tratamiento de textos (training "training_text"):
# --------------------------------------------------------
if(b_hay_cambios | !mi_load(s_input_path, "training_text_features.RData"))
{
  if(!mi_load(s_input_path, "training_text.RData"))
  {
    jjprint('Leyendo textos training')
    library(readr)
    full_trainset_txt <- readr::read_lines(file = file.path(s_input_path, "training_text"), skip = 1, n_max = ifelse(G_b_DEBUG, 10L, -1L))
    full_trainset_txt <- str_split(full_trainset_txt, '\\|\\|', n = 2, simplify = TRUE)
    full_trainset_txt <- data.table(ID = as.integer(full_trainset_txt[,1]), texto = full_trainset_txt[,2])
    setkey(full_trainset_txt, texto)
    full_trainset_txt[, train_texto_id := .GRP, by = texto] # Creamos train_texto_id
    save(full_trainset_txt, file = file.path(s_input_path, "training_text.RData"))
  }
  # Quitamos texto original (ocupa mucho!):
  full_trainset_txt[, texto := NULL]; gc()
  setkey(full_trainset_txt, train_texto_id)
  # Preparamos "training_text_features.RData":
  # mi_source("Redef_Cancer_crearTextFeatures.R")
  fichs_LDA <- c('LDAf_1_17_54_54_1-3-0-1_1-1-0-0_1000_1921_LDA_0-54.RData' # V.3
                ,'LDAp_1_17_54_54_1-3-0-1_1-1-0-0_1000_1921_LDA_54-0.RData' # V.4
                ,'LDA_1_17_27_27_1-3-0-1_1-1-0-0_1000_400_LDA_27-27.RData' # V.5 y 6
                ,'LDAf_1_17_72_72_1-3-0-2_1-1-0-0_1000_1921_LDA_0-72.RData' # V.7 y 8
                # ,'LDAp_0.16_17_72_72_1-3-0-2_1-2-0-1_1000_307_LDA_72-0.RData' # V.7 y 8 [before 27/08/2017 22:00]
                ,'LDAp_1_17_72_72_1-3-0-2_1-1-0-0_1000_400_LDA_72-0.RData'   # V.7 y 8 [after  27/08/2017 22:00]
                ,'LDAc_1_17_72_72_1-3-0-2_1-1-0-0_1-5-0-0_100_1921_LDA_0-0-72' # V.9, 10 y 11
                ,'LDAc_1_17_45_45_1-3-0-2_1-1-0-0_1-6-0-0_100_1921_LDA_0-0-45' # V.12
  )
  for(fich_LDA in fichs_LDA)
  {
    stopifnot(mi_load(s_input_path, fich_LDA)) # Cargar Fich_TRN_txt y Fich_TST_txt
    jjprint(paste0('Ok. Fusionando features de textos en full_trainset_txt [',fich_LDA,']...'))
    full_trainset_txt <- merge(full_trainset_txt, Fich_TRN_txt, by = "train_texto_id")
  }
  setkey(full_trainset_txt, ID)
  save(full_trainset_txt, file = file.path(s_input_path, "training_text_features.RData"))
  b_hay_cambios <- TRUE
}
# --------------------------------------------------------
# Tratamiento de textos (test "test_text"):
# --------------------------------------------------------
if(b_hay_cambios | !mi_load(s_input_path, "test_text_features.RData"))
{
  if(!mi_load(s_input_path, "test_text.RData"))
  {
    jjprint('Leyendo textos test')
    library(readr)
    testset_txt <- readr::read_lines(file = file.path(s_input_path, "test_text"), skip = 1, n_max = ifelse(G_b_DEBUG, 10L, -1L))
    testset_txt <- str_split(testset_txt, '\\|\\|', n = 2, simplify = TRUE)
    testset_txt <- data.table(ID = as.integer(testset_txt[,1]), texto = testset_txt[,2])
    setkey(testset_txt, texto)
    testset_txt[, test_texto_id := .GRP, by = texto] # Creamos test_texto_id
    save(testset_txt, file = file.path(s_input_path, "test_text.RData"))
  }
  # Quitamos texto original (ocupa mucho!):
  testset_txt[, texto := NULL]; rm(Fich_TST_txt); gc()
  setkey(testset_txt, test_texto_id)
  # Preparamos "test_text_features.RData":
  # mi_source("Redef_Cancer_crearTextFeatures.R")
  stopifnot(exists("fichs_LDA"))
  for(fich_LDA in fichs_LDA)
  {
    stopifnot(mi_load(s_input_path, fich_LDA)) # Cargar Fich_TRN_txt y Fich_TST_txt
    jjprint(paste0('Ok. Fusionando features de textos en testset_txt [',fich_LDA,']...'))
    testset_txt <- merge(testset_txt, Fich_TST_txt, by = "test_texto_id")
  }
  setkey(testset_txt, ID)
  save(testset_txt, file = file.path(s_input_path, "test_text_features.RData"))
  b_hay_cambios <- TRUE
}
# --------------------------------------------------------
# CARGAMOS EL ÚLTIMO BATCH PARA EMPEZAR CON ALGO:
# --------------------------------------------------------
s_fich_train <- get_batch_train_filename(NUM_BLOQUES)
s_fich_test <- get_batch_test_filename(NUM_BLOQUES)
# Si no existe full_trainset_016, los creamos todos:
if(b_hay_cambios | !file.exists(file.path(s_input_path, s_fich_train)))
{
  if(!G_b_DEBUG)
  {
    jjprint('Leyendo trainset original...')
    full_trainset <- fread(file.path(s_input_path, "training_variants"))
    full_trainset[, Gene_orig := Gene]
    full_trainset[, Variation_orig := Variation]
    full_trainset[, Gene := factor(Gene, levels=v_Genes, labels=v_Genes, ordered=F)]
    full_trainset[, Variation := factor(Variation, levels=v_Variations, labels=v_Variations, ordered=F)]
    full_trainset[is.na(Gene), Gene := '--']
    full_trainset[is.na(Variation), Variation := '--']
    setkey(full_trainset, ID)
    full_trainset <- merge(full_trainset, full_trainset_txt, by = "ID") #, all.x = TRUE)
    rm(full_trainset_txt)
    jjprint('Leyendo train_Variations')
    full_trainset <- merge(full_trainset, train_Variations, by = "ID") #, all.x = TRUE)
    print('Splitting full_trainset...')
    mi_split_train(NUM_BLOQUES)
    b_hay_cambios <- TRUE
    rm(full_trainset); gc()
  }
}
# Si no existe testset_016 (o testset_016_016), los creamos todos:
if(b_hay_cambios | !file.exists(file.path(s_input_path, s_fich_test)) & !file.exists(paste0(s_input_path, get_batch_test_filename(NUM_BLOQUES, NUM_BLOQUES))))
{
  # Si no existe testset_016, los creamos todos:
  load(file = file.path(s_input_path, 'Gene_Variation_Factors.RData'))
  if(!G_b_DEBUG)
  {
    jjprint('Leyendo testset original...')
    testset <- fread(file.path(s_input_path, "test_variants"))
    testset[, Gene_orig := Gene]
    testset[, Variation_orig := Variation]
    testset[, Gene := factor(Gene, levels=v_Genes, labels=v_Genes, ordered=F)]
    testset[, Variation := factor(Variation, levels=v_Variations, labels=v_Variations, ordered=F)]
    setkey(testset, ID)
    testset <- merge(testset, testset_txt, by = "ID") #, all.x = TRUE)
    rm(testset_txt)
    jjprint('Leyendo test_Variations')
    testset <- merge(testset, test_Variations, by = "ID") #, all.x = TRUE)
    print('Splitting testset...')
    testset[, c(paste0("prob_class", 1:9)) := as.numeric(0)]
    mi_split_test(NUM_BLOQUES)
    b_hay_cambios <- TRUE
    rm(testset); gc()
  }
}
if(exists("full_trainset")) rm(full_trainset)
if(exists("testset"))       rm(testset)
gc()
full_trainset <- leer_batch_train(NUM_BLOQUES, "inicio", s_input_path)
if( file.exists(paste0(s_input_path, get_batch_test_filename(NUM_BLOQUES, NUM_BLOQUES))))
  testset <- leer_batch_test(NUM_BLOQUES, "inicio", s_input_path, numSubBatch = NUM_BLOQUES)
if(!file.exists(paste0(s_input_path, get_batch_test_filename(NUM_BLOQUES, NUM_BLOQUES))))
  testset <- leer_batch_test(NUM_BLOQUES, "inicio", s_input_path)

# ------------------------------------------------------------------------------------
if(b_hay_cambios) # Guardamos train_valid_nn__ssss__pppp.RData:
{
  nn <- 1 # numCluster
  ssss <- 1234 # xgb_reduc_seed
  n_porc <- c(0.699, 0.75, 0.85, 0.99, 1) # xgb_train_porc
  # Leemos full_trainset COMPLETO:
  miDescr <- 'Leyendo full_trainset COMPLETO'
  full_trainset <- leer_batch_train(numBatch = 1, miDescr, s_input_path)
  for(numBatch in 2:NUM_BLOQUES)
    full_trainset <- rbindlist(list(full_trainset, leer_batch_train(numBatch, miDescr, s_input_path)))
  jjprint(paste0('Ok. full_trainset COMPLETO: ', nrow(full_trainset), ' registros.'))
  # Dividir muestra en entrenamiento y validación: NOTA: Usamos ID para el corte...
  foreach(mi_porc = n_porc, .packages = c('data.table','stringr')
  )  %do%  { #  %dopar%  {
    fich_name <- paste0("train_valid_", str_pad(nn, 2, "left", "0"), "__", ssss, "__", mi_porc, ".RData") # "train_valid_nn__ssss__pppp.RData"
    jjprint(paste0('preparando [', fich_name, ']...'))
    trainset <- reducir_trainset(mi_set = full_trainset, id_fld = "ID", n_seed = ssss, n_porc = mi_porc)
    validset <- trainset[[2]]
    trainset <- trainset[[1]]
    # Guardamos:
    # if(!file.exists(file.path(s_input_path, fich_name)))
    {
      jjprint(paste0('Guardando [', fich_name, '] en ', s_input_path))
      save(trainset, validset, file = file.path(s_input_path, fich_name))
    }
    return('Ok.')
  }
  for(mi_porc in n_porc){
    fich_name <- paste0("train_valid_", str_pad(nn, 2, "left", "0"), "__", ssss, "__", mi_porc, ".RData") # "train_valid_nn__ssss__pppp.RData"
    load(file = file.path(s_input_path, fich_name))
    # jjprint(paste0(mi_porc, ' - ', ncol(trainset), ' variables.'))
    # print(sapply(trainset, uniqueN))
    # if(!is.null(validset)) if(nrow(validset)>0) print(sapply(validset, uniqueN))
    jjprint(paste0(mi_porc, ' - ', ncol(trainset), ' variables.'))
    jjprint(jjfmt(table(trainset$Class)/nrow(trainset), num_decimals = 2), b_add_datetime = F)
    if(!is.null(validset)) if(nrow(validset)>0)
      jjprint(jjfmt(table(validset$Class)/nrow(validset), num_decimals = 2), b_add_datetime = F)
  }
}

# full_trainset <- fread(file.path(s_input_path, "training_variants"))
# basic_preds_guardar_submit(trainset = full_trainset, k = 2, reg = 4, s_output_path = s_output_path)
# foreach(reg = c(0,3,4,5,6), .inorder=TRUE, .combine = c, .packages = c('data.table'),
#         .export = c('Map_12_diario', 'Map_12', 'mapk_nodup_only_first_12', 'add_tiempo')
# ) %do%
# {
#   return(
#     basic_preds_guardar_submit(trainset = full_trainset, k = 3, reg = reg, s_output_path = s_output_path)
#   ) # ret_val de la función "foreach() %do%" (o foreach( %dopar%))
# }

# # 
# # M12 <- foreach(reg = c(0,3,4,5,6), .inorder=TRUE, .combine = c, .packages = c('data.table'),
# #                .export = c('Map_12_diario', 'Map_12', 'mapk_nodup_only_first_12', 'add_tiempo')
# #                ) %do%
# # {
# #   return(
# #     basic_preds_m12(trainset = trainset, validset = validset, k = 2, reg = reg, b_verbose = 1)
# #     ) # ret_val de la función "foreach() %do%" (o foreach( %dopar%))
# # }
# # print(M12)
# # sort(M12)
# # sort(unlist(sapply(M12, function(x) {return(x$V1[2])})))
# # sort(unlist(sapply(M12, function(x) {return(x$V2[2])})))
# # sort(unlist(sapply(M12, function(x) {return(x$V3[2])})))

minutos_acum <- as.double((proc.time() - systime_ini)['elapsed'])/60
jjprint(paste0(minutos_acum, ' minutos en total.'))
jjprint(paste0('Ok.', ifelse(b_hay_cambios, '', ' (sin cambios)')))

# cleanup:
try(registerDoSEQ(), silent = TRUE) # library(doParallel) [turn parallel processing off and run sequentially again]
try(stopImplicitCluster(), silent = TRUE)
try(stopCluster(cl), silent = TRUE)

