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
s_modelos_path <- paste0(s_input_path, 'modelos/') # <- s_output_path

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
# NOTA: Dejamos un Core de la CPU "libre" para no "quemar" la máquina:
cl <- makeCluster(detectCores(), type='PSOCK') # library(doParallel) [turn parallel processing on]
registerDoParallel(cl) # library(doParallel) [turn parallel processing on]

# memory.limit(size = 16000)

systime_ini <- proc.time()

cmdargs = commandArgs(trailingOnly=TRUE)
if(length(cmdargs) >= 1) print(paste0('cmdargs = [', paste(cmdargs, collapse = ' '),']'))
G_b_NO_LOG <- any(cmdargs == "/q"); cmdargs <- cmdargs[cmdargs != "/q"]
# 1 - b_crearSubmit    # 0 (0, 1)
# 2 - max_depth        # 6 [1 - inf] maximum depth of the tree (demasiado alto => overfitting)
# 3 - colsample_bytree # 1 (0 - 1] (muestreo aleatorio por columnas si es < 1, para intentar evitar overfitting)
# 4 - min_child_weight # 1 (0 - 10]
# 5 - gamma            # 0 [0 - 2] the larger, the more conservative (stable, slow-learner) the algorithm will be.
# 6 - num_parallel_tree# 1 [1 - inf] Experimental. Random trees at each round. Helps reducing over-fitting. Prueba que, con pocos datos, parece mejorar (aunque tarda mucho más)
cmdargs_xgb_max_depth = ifelse(length(cmdargs) >= 2, as.numeric(cmdargs[2]), 6) # 10 [1 - inf] maximum depth of the tree (demasiado alto => overfitting)
cmdargs_xgb_colsample_bytree = ifelse(length(cmdargs) >= 3, as.numeric(cmdargs[3]), 1) # (muestreo aleatorio por columnas si es < 1, para intentar evitar overfitting)
cmdargs_xgb_min_child_weight = ifelse(length(cmdargs) >= 4, as.numeric(cmdargs[4]), 1) # + (numCluster/2 - 1) # 1 [0.1 - 10]
cmdargs_xgb_gamma = ifelse(length(cmdargs) >= 5, as.numeric(cmdargs[5]), 1) # 0 [0 - 2] the larger, the more conservative the algorithm will be.
cmdargs_xgb_num_parallel_tree = ifelse(length(cmdargs) >= 6, as.numeric(cmdargs[6]), 1) # 1 [1 - inf] Experimental. Helps reducing over-fitting. Prueba que, con pocos datos, parece mejorar (aunque tarda mucho más)
xgb_mascara_filename <- paste0(paste('.*',
                             ifelse(length(cmdargs) >= 2, cmdargs_xgb_max_depth, '.*'),
                             '.*',
                             ifelse(length(cmdargs) >= 3, cmdargs_xgb_colsample_bytree, '.*'),
                             '.*',
                             ifelse(length(cmdargs) >= 4, cmdargs_xgb_min_child_weight, '.*'),
                             ifelse(length(cmdargs) >= 5, cmdargs_xgb_gamma, '.*'),
                             ifelse(length(cmdargs) >= 6, cmdargs_xgb_num_parallel_tree, '.*'),
                             '.*',
                             sep = "_"))
# --------------------------------------------------------
G_b_DEBUG <- FALSE # Reducimos todo para hacer pruebas más rápido
G_b_REV <- FALSE # Empezamos por el final (numModelo <- NUM_MODELOS - forNumModelo + 1)
G_maxNumFeatures_OneHotEnc <- 10000 # One-Hot Encode las categóricas con estos levels como máx.
G_CV_nFolds <- 10 # 5 [2 - Inf]
G_b_crearSubmit <- ifelse(length(cmdargs) >= 1, (cmdargs[1]==1), FALSE)
NUM_BLOQUES <- 3
B_CON_CLUSTER <- FALSE # Tiene que haber un campo "Cluster" para ponerlo a TRUE ["Cluster" fue "numAds" en Outbrain]
NUM_MODELOS <- ifelse(B_CON_CLUSTER, 9, 1) # Un modelo por cada "Cluster"...
maxTamFullTrainset <- 100000 # 200.000 para hacerlo "rápido" (y crear las importance_matrix) ó 3.000.000 que es lo máx. que he podido hasta ahora (con muy pocas vars.)
maxImportanceNumVars <- 0 # 100 # Seleccionamos las primeras variables por importancia (de ejecuciones anteriores con misma versión y numAd)
# --------------------------------------------------------
gc()
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# https://www.analyticsvidhya.com/blog/2016/01/xgboost-algorithm-easy-steps/
# install.packages("xgboost")
suppressPackageStartupMessages(library(xgboost))

# suppressPackageStartupMessages(library(lattice)) # Para xgb_prep_datos
# suppressPackageStartupMessages(library(ggplot2)) # Para xgb_prep_datos
# suppressPackageStartupMessages(library(caret))   # Para xgb_prep_datos
# -------------------------------
# Inicializamos:
# -------------------------------
xgb_scale_pos_weight <- 1 # default
# # # mi_sum_pos <- sum(full_trainset$clicked == 1) # 16874593
# # # mi_sum_neg <- sum(full_trainset$clicked == 0) # 70267138
# # mi_sum_pos <- 16874593
# # mi_sum_neg <- 70267138
# xgb_scale_pos_weight <- mi_sum_neg / mi_sum_pos # == 4.164079
# -------------------------------
b_conPond <- 0 # Sin ponderar (V.100)
# Ponderación cancelada: afecta a la métrica de cálculo del error, i.e. mlogloss
# y no parece mejorar en nada el entrenamiento.
# En RedefCancer explícitamente avisan que no se tenga en cuenta la frecuencia de cada clase, así que razón de más...
# A estudiar en el futuro...
# # full_trainset <- fread(file.path(s_input_path, "training_variants"))
# # weights_clases <- table(full_trainset$Class)/nrow(full_trainset)
# # paste0('xgb_weights_clases <- c(', paste0(paste(paste0("'",names(weights_clases),"'"), weights_clases, sep = '='), collapse = ','), ')')
# # xgb_weights_clases <- c('1'=0.171032821439325,'2'=0.136103583258055,'3'=0.0267991568804577,'4'=0.206564287865101,'5'=0.0728696175850647,'6'=0.0828063836193917,'7'=0.286961758506474,'8'=0.00572116832279434,'9'=0.0111412225233363)
# # b_conPond <- 1 # V.200
# # # xgb_weights_clases200 <- 1 / xgb_weights_clases
# # # xgb_weights_clases200 <- xgb_weights_clases200 / sum(xgb_weights_clases200)
# # # b_conPond <- 2 # V.300
# -------------------------------
xgb_num_class = 9 # xgb_objective = "multi:softprob"
xgb_reduc_seed <- 1234
xgb_train_porc <- 1 # 1 # 0.99 # 0.85
xgb_get_predictions <- FALSE # TRUE para hacer stacking de los modelos (Subsemble, SuperLearner)

if(G_b_crearSubmit) {
  print('-------------------------------')
  print(paste0('Primero buscamos algún modelo [', xgb_mascara_filename, '] entrenado pero SIN submit:'))
  print('-------------------------------')
}
xgb_tipo_modelo_XGB = "XGB"    # XGBoost
xgb_tipo_modelo_XGRF = "XGRF"  # Random Forest con XGBoost
xgb_modelos <- vector(mode = "character", length = NUM_MODELOS)
xgb_filenames <- vector(mode = "character")
if(G_b_crearSubmit) # Solo los modelos de los cmdargs
  for(mi_tipo in paste0(c(xgb_tipo_modelo_XGRF, xgb_tipo_modelo_XGB), xgb_mascara_filename))
    xgb_filenames <- c(xgb_filenames, unique(str_replace(dir(path = s_modelos_path, pattern = paste0(mi_tipo, '.*.modelo')), "_[0-9]*\\.[0-9][0-9][0-9]\\.modelo$", "")))
for(mi_tipo in xgb_filenames)
{
  if(length(dir(path = s_output_path, pattern = paste0(mi_tipo, '_submit.*.csv'))) != 0)
  { xgb_filenames <- xgb_filenames[xgb_filenames != mi_tipo]; next }
  if(length(dir(path = s_modelos_path, pattern = paste0(mi_tipo, '_submit.*.csv'))) != 0)
  { xgb_filenames <- xgb_filenames[xgb_filenames != mi_tipo]; next }
  if(length(dir(path = s_output_path, pattern = paste0(mi_tipo, '_submit.*.zip'))) != 0)
  { xgb_filenames <- xgb_filenames[xgb_filenames != mi_tipo]; next }
  if(length(dir(path = s_modelos_path, pattern = paste0(mi_tipo, '_submit.*.zip'))) != 0)
  { xgb_filenames <- xgb_filenames[xgb_filenames != mi_tipo]; next }
  if(length(dir(path = s_output_path, pattern = paste0(mi_tipo, '_submit.*.rar'))) != 0)
  { xgb_filenames <- xgb_filenames[xgb_filenames != mi_tipo]; next }
  if(length(dir(path = s_modelos_path, pattern = paste0(mi_tipo, '_submit.*.rar'))) != 0)
  { xgb_filenames <- xgb_filenames[xgb_filenames != mi_tipo]; next }
}
for(xgb_filename in xgb_filenames)
{
  print('Encontrado un modelo sin submit...')
  print(xgb_filename)
  tmp_xgb_modelos <- str_replace(dir(path = s_modelos_path, pattern = paste0(str_replace(xgb_filename, "_[0-9]*\\.[0-9][0-9][0-9]$", ""), '.*.modelo')), pattern = "\\.modelo", replacement = "")
  n_version <- as.integer(substring(str_extract(xgb_filename, pattern = "v[0-9]+"), first = 2))
  if(n_version > 500)
    next # Es de los de numAds
  if(length(tmp_xgb_modelos) == NUM_MODELOS)
  {
    # Ordenamos modelos por numModelo (aquí NO da igual, porque ya NO los vamos a promediar)
    xgb_modelos[as.integer(substr(tmp_xgb_modelos, nchar(tmp_xgb_modelos)-2, nchar(tmp_xgb_modelos)))] <- tmp_xgb_modelos
    break # Ok. Pasamos directamente a predecir con estos modelos.
  } else {
    print('NOTA: Lo descartamos porque no están todos entrenados')
  }
}
n_versiones <- as.integer(substring(str_extract(xgb_filenames, pattern = "v[0-9]+"), first = 2))
for(n_version in unique(n_versiones)[unique(n_versiones) > 500])
{
  nSubmits <- length(dir(path = s_modelos_path, pattern = paste0('.*', xgb_tipo_modelo_XGRF,'.*_v', n_version, '.*.submit.*')))
  nSubmits <- nSubmits + length(dir(path = s_modelos_path, pattern = paste0('.*', xgb_tipo_modelo_XGB,'.*_v', n_version, '.*.submit.*')))
  if(nSubmits != 0)
  {
    print(paste0('NOTA: Descartamos v', n_version, '. Encontrado(s) ', nSubmits, ' submit(s) con esta versión.'))
    break
  } else { print(paste0('Buscando modelos con versión ', n_version)) }
  tmp_xgb_modelos <- vector(mode = "character", length = NUM_MODELOS)
  xgb_filenames_version <- xgb_filenames[n_versiones == n_version]
  for(xgb_filename in xgb_filenames)
  {
    if(n_version != as.integer(substring(str_extract(xgb_filename, pattern = "v[0-9]+"), first = 2)))
      next
    mis_modelos <- dir(path = s_modelos_path, pattern = paste0(xgb_filename, '.*.modelo'))
    for(mi_modelo in mis_modelos)
    {
      numAds <- 1 + as.integer(substring(str_extract(mi_modelo, pattern = "\\.[0-9][0-9][0-9]\\.modelo"), first = 2, last = 4))
      if(tmp_xgb_modelos[numAds - 1] != ""){
        print(paste0('Nota: Encontrado más de un modelo ', str_pad(numAds-1,3,"left","0"), ' v', n_version, ' (numAds = ', numAds, '). Nos quedamos con el menor...'))
        if(tmp_xgb_modelos[numAds - 1] > str_replace(mi_modelo, pattern = "\\.modelo", replacement = ""))
          tmp_xgb_modelos[numAds - 1] <- str_replace(mi_modelo, pattern = "\\.modelo", replacement = "")
      } else {
        tmp_xgb_modelos[numAds - 1] <- str_replace(mi_modelo, pattern = "\\.modelo", replacement = "")
      }
    }
  }
  if(length(tmp_xgb_modelos[tmp_xgb_modelos != ""]) == NUM_MODELOS)
  {
    # Ordenamos modelos por numModelo (aquí NO da igual, porque ya NO los vamos a promediar)
    xgb_modelos[as.integer(substr(tmp_xgb_modelos, nchar(tmp_xgb_modelos)-2, nchar(tmp_xgb_modelos)))] <- tmp_xgb_modelos
    print(paste0('Hay ', NUM_MODELOS, ' modelos v', n_version, ' (sin submit). Pasamos directamente a predecir con ellos [', xgb_filename, '].'))
    print(xgb_modelos)
    break # Ok. Pasamos directamente a predecir con estos modelos.
  } else {
    print(paste0('NOTA: Descartamos v', n_version, ' porque no están todos entrenados (faltan algunos numAds)'))
    print(substring(str_extract(tmp_xgb_modelos, pattern = "\\.[0-9][0-9][0-9]$"), first = 2, last = 4))
  }
}
s_Fic_log <- NULL # Por si acaso no hemos entrado a entrenar...
if(any(xgb_modelos == "") | anyNA(xgb_modelos))
{
  print('-------------------------------')
  print('Entrenamos:')
  print('-------------------------------')

  for(forNumModelo in 1:NUM_MODELOS) # Un modelo por cada "Cluster"...
  { #PRUEBAS forNumModelo <- 1 #PRUEBAS
    
    if(!G_b_REV) numModelo <- forNumModelo else numModelo <- (NUM_MODELOS - forNumModelo + 1) # Empezamos por el final (numModelo <- NUM_MODELOS - forNumModelo + 1)

    if(!is.na(xgb_modelos[numModelo]) & xgb_modelos[numModelo] != "")
    {
      print(paste0('Warning: Ya hay un modelo ', str_pad(numModelo, 3, "left", "0"), ' (', xgb_modelos[numModelo], '). Pasamos al siguiente...'))
      next # el "continue" de C
    }
    miDescr <- paste0("XGB Training - [cluster ", numModelo, "]")
    numCluster <- numModelo # numModelo + 1 # Era numAdsCluster
    fich_name <- paste0("train_valid_", str_pad(numCluster, 2, "left", "0"), "__", xgb_reduc_seed, "__", xgb_train_porc, ".RData") # "train_valid_nn__ssss__pppp.RData"
    
    if(file.exists(file.path(s_input_path, fich_name)))
    {
      print(paste0(miDescr, ' Batch trainset+validset ALL (cluster == ', numCluster, ')'))
      load(file = file.path(s_input_path, fich_name))
      print(paste0(miDescr, ' Batch trainset+validset ALL (cluster == ', numCluster, ')', ' Ok. ', jjfmt(nrow(trainset)), ' + ', jjfmt(nrow(validset)), ' regs.'))
    } else {
  		fich_name <- paste0("full_trainset_", str_pad(numCluster, 2, "left", "0"), ".RData") # "full_trainset_nn.RData"
  		
  		if(file.exists(file.path(s_input_path, fich_name)))
  		{
  		  print(paste0(miDescr, ' Batch trainset ALL (cluster == ', numCluster, ')'))
  		  load(file = file.path(s_input_path, fich_name))
  		  print(paste0(miDescr, ' Batch trainset ALL (cluster == ', numCluster, ')', ' Ok. ', nrow(full_trainset), ' regs.'))
  		  numBatch <- NUM_BLOQUES # Para dejar claro que hemos cargado todos los ficheros!
  		} else {
  		  numBatch <- 1
  		  fich_name <- paste0("full_trainset_", str_pad(numCluster, 2, "left", "0"), '_', str_pad(numBatch, 3, "left", "0"), ".RData") # "full_trainset_nn_001.RData"
  		  
  		  if(!file.exists(file.path(s_input_path, fich_name)))
  		  {
  		    if(B_CON_CLUSTER)  full_trainset <- leer_batch_train(numBatch, miDescr, s_input_path)[Cluster == numCluster,]#[numAds == numCluster,]
  		    if(!B_CON_CLUSTER) full_trainset <- leer_batch_train(numBatch, miDescr, s_input_path)
  		  } else {
    			# print(paste0(miDescr, ' Batch trainset ', numBatch, ' (Cluster == ', numCluster, ')'))
    			load(file = file.path(s_input_path, fich_name))
  		  }
  		  print(paste0(miDescr, ' Batch trainset ', numBatch, ' (Cluster == ', numCluster, ')', ' Ok. ', nrow(full_trainset), ' regs.'))
  		  nrowsPorBloqueEstim <- nrow(full_trainset)
  		  for(numBatch in 2:NUM_BLOQUES)
  		  {
    			# if(file.exists(paste0(s_input_path, "clicks_full_trainset_2x_debug.csv")))
    			# {
    			#   if(numBatch > 1)
    			#      break # Finished!
    			#   jjprint(paste0('XGB Training (extreme gradient boosting) - PRUEBA CON clicks_full_trainset_2x_debug.csv...'))
    			#   full_trainset <- fread(paste0(s_input_path, "clicks_full_trainset_2x_debug.csv"))
    			# } else {
    				# if(nrow(full_trainset) + nrowsPorBloqueEstim > ifelse(G_b_DEBUG, 50000, maxTamFullTrainset))
    				# {
    				#   print(paste0("NOTA: Modelo ", str_pad(numModelo, 3, "left", "0")," (Cluster == ", numCluster, ") NO está completo! (Solo tiene ", numBatch-1, "/", NUM_BLOQUES, " batches)"))
    				#   break # Paramos antes de pasarnos...
    				# }
  				fich_name <- paste0("full_trainset_", str_pad(numCluster, 2, "left", "0"), '_', str_pad(numBatch, 3, "left", "0"), ".RData")
  	  
  				full_trainset_tmp <- full_trainset
  				if(!file.exists(paste0(s_input_path, fich_name)))
  				{
  				  if(B_CON_CLUSTER)  full_trainset_tmp <- rbindlist(list(full_trainset_tmp, leer_batch_train(numBatch, miDescr, s_input_path)[Cluster == numCluster,]))#[numAds == numCluster,]
  				  if(!B_CON_CLUSTER) full_trainset_tmp <- rbindlist(list(full_trainset_tmp, leer_batch_train(numBatch, miDescr, s_input_path)))
  				} else {
  				  # print(paste0(miDescr, 'Batch trainset ', numBatch, ' (Cluster == ', numCluster, ')'))
  				  load(file = paste0(s_input_path, fich_name))
  				  print(paste0(miDescr, 'Batch trainset ', numBatch, ' (Cluster == ', numCluster, ')', ' Ok. ', nrow(full_trainset), ' regs.'))
  				  full_trainset_tmp <- rbindlist(list(full_trainset_tmp, full_trainset))
  				}
  				full_trainset <- full_trainset_tmp
  				jjprint(paste0('Ok. full_trainset: ', nrow(full_trainset), ' registros.'))
  				if(nrow(full_trainset) > ifelse(G_b_DEBUG, 50000, maxTamFullTrainset))
  				{
  				  if(numBatch < NUM_BLOQUES)
  					print(paste0("NOTA: Modelo ", str_pad(numModelo, 3, "left", "0")," (Cluster == ", numCluster, ") NO está completo! (Solo tiene ", numBatch, "/", NUM_BLOQUES, " batches)"))
  				  break # Paramos antes de pasarnos...
  				}
  				
    			#   if(!G_b_DEBUG)
    			#   {
    			#     # Ampliamos con uno más:
    			#     if(numBatch > 1)
    			#       break # Finished!
    			#     full_trainset <- rbind(full_trainset, leer_batch_train(numBatch + 1, "XGB Training (extreme gradient boosting) 2x Blocks", s_input_path))
    			#   }
    			# }
    			# # sapply(trainset, uniqueN)[sapply(trainset, uniqueN)==1]
    			# # sapply(full_trainset, uniqueN)
    			# # print(100 * round(table(full_trainset$dia)) / nrow(full_trainset), digits = 3)
  		  }
  		  if(exists("full_trainset_tmp")) { rm(full_trainset_tmp); gc() }
  		  if(numBatch == NUM_BLOQUES)
  		  {
    			# Hemos cargado todos los ficheros, así que podemos guardarlo en uno solo:
    			fich_name <- paste0("full_trainset_", str_pad(numCluster, 2, "left", "0"), ".RData") # "full_trainset_nn.RData"
    			save(full_trainset, file = file.path(s_input_path, fich_name))
  		  }
  		}
  		
  		if(G_b_DEBUG)
  		{
  		  print('Hacemos un sample para que vaya más rápido...')
  		  full_trainset <- reducir_trainset(mi_set = full_trainset, n_seed = xgb_reduc_seed, n_porc = 0.1)[[1]]
  		}
  		if(nrow(full_trainset) > maxTamFullTrainset)
  		{
  		  print(paste0('Hacemos un sample (', jjfmt(maxTamFullTrainset), ') para que vaya más rápido...'))
  		  full_trainset <- reducir_trainset(mi_set = full_trainset, n_seed = xgb_reduc_seed, n_porc = (maxTamFullTrainset /  nrow(full_trainset)))[[1]]
  		}
  		# if(!file.exists(paste0(s_input_path, "clicks_full_trainset_2x_debug.csv")))
  		# {
  		#   jjprint(paste0('Guardando full_trainset en ', s_input_path, 'clicks_full_trainset_2x_debug.csv...'))
  		#   write.table(full_trainset, file = paste0(s_input_path, "clicks_full_trainset_2x_debug.csv"), row.names=F, quote=F, sep=",")
  		# }
  		# jjprint(paste0('Ok. full_trainset: ', nrow(full_trainset), ' registros.'))
  		reduc_list <- reducir_trainset(mi_set = full_trainset, id_fld = 'ID', n_seed = xgb_reduc_seed, n_porc = xgb_train_porc)
  		# Reducimos uso de memoria:
  		if(!G_b_DEBUG)  rm(full_trainset); gc()
  		trainset <- reduc_list[[1]]
  		validset <- reduc_list[[2]]
  		setkey(trainset, ID)
  		# Reducimos uso de memoria:
  		if(!G_b_DEBUG)  rm(reduc_list); gc()
  		if(numBatch == NUM_BLOQUES)
  		{
  			# Hemos cargado todos los ficheros, así que podemos guardar trainset y validset directamente:
  			fich_name <- paste0("train_valid_", str_pad(numCluster, 2, "left", "0"), "__", xgb_reduc_seed, "__", xgb_train_porc, ".RData") # "train_valid_nn__ssss__pppp.RData"
  			save(trainset, validset, file = file.path(s_input_path, fich_name))
  		}
  	}
    # print(100 * round(table(trainset$dia)) / nrow(trainset), digits = 3)
    # print(100 * round(table(validset$dia)) / nrow(validset), digits = 3)

    if(nrow(trainset) + nrow(validset) > maxTamFullTrainset)
    {
      print(paste0('Hacemos un sample (', jjfmt(maxTamFullTrainset), ') para que vaya más rápido...'))
      mi_n_porc = maxTamFullTrainset / (nrow(trainset) + nrow(validset))
      trainset <- reducir_trainset(mi_set = trainset, n_seed = xgb_reduc_seed, n_porc = mi_n_porc)[[1]]
      validset <- reducir_trainset(mi_set = validset, n_seed = xgb_reduc_seed, n_porc = mi_n_porc)[[1]]
    }

    jjprint(paste0('Ok. trainset: ', jjfmt(nrow(trainset)), ' registros.  validset: ', jjfmt(nrow(validset)), ' registros.'))
    
    # sapply(trainset, uniqueN)[sapply(trainset, uniqueN)==1]
    
    xgb_preps <- xgb_prep_datos(mi_dt = NULL, b_verbose = 1) # Obtenemos solamente la versión (de 0 a 99!!!)
    
    # ========================  XGB params - START  ==========================
    
    # xgb_mi_version <- 000 + xgb_preps[[1]] # Esto viene de xgb_prep_datos (indica cambios en las "features" seleccionadas)
    xgb_mi_version <- 100 + xgb_preps[[1]] # V.100 Sin usar los pesos (xgb_weights_clases) en el training, en xgb.cv() y xgboost())
    # if(b_conPond == 1) { xgb_mi_version <- xgb_mi_version + 100 } else if(b_conPond == 2) { xgb_mi_version <- xgb_mi_version + 200 } # V.200 / V.300
    if(maxImportanceNumVars != 0) xgb_mi_version <- 400 + xgb_preps[[1]] # V.400 Con selección de variables (maxImportanceNumVars=100 mejores)
    
    # xgb_num_class = 9 # xgb_objective = "multi:softprob"
    xgb_objective = "multi:softprob" # same as softmax, but prediction outputs a vector of ndata * nclass elements, which can be further reshaped to ndata, nclass matrix. The result contains predicted probabilities of each data point belonging to each class.
          # NOTA para multi-class: Class is represented by a number and should be from 0 to num_class - 1.
    # xgb_objective = "binary:logistic" # logistic regression for binary classification, output probability
    # xgb_objective = "rank:pairwise" # set XGBoost to do ranking task by minimizing the pairwise loss
    xgb_eta = 0.01 # 0.3 [0.01 - 0.2] step size of each boosting step (si baja, afina mejor pero habrá que subir nround)
    xgb_max_depth = cmdargs_xgb_max_depth # 10 [1 - inf] maximum depth of the tree (demasiado alto => overfitting)
    xgb_nround = 10000 # max number of iterations
    if(G_b_DEBUG & xgb_nround > 1)  xgb_nround = 100
    xgb_subsample = 1 # (muestreo aleatorio por filas si es < 1, para intentar evitar overfitting)
    xgb_colsample_bytree = cmdargs_xgb_colsample_bytree # (muestreo aleatorio por columnas si es < 1, para intentar evitar overfitting)
    # xgb_eval_metric = "rmse": root mean square error, "mae" = mean absolut error, "logloss": negative log-likelihood
    # xgb_eval_metric = "map"
    # xgb_eval_metric = "error"
    # xgb_eval_metric = "mae"
    # xgb_eval_metric = "logloss"
    # xgb_eval_metric = "auc"
    xgb_eval_metric = "mlogloss" # Para xgb_objective = "multi:softprob"
    # xgb_eval_metric = "merror" # Para xgb_objective = "multi:softprob"
    # xgb_scale_pos_weight = numCluster # Cambio por numAds clustering # mi_sum_neg / mi_sum_pos # == 4.164079
    xgb_min_child_weight = cmdargs_xgb_min_child_weight # + (numCluster/2 - 1) # 1 [0.1 - 10]
    
    xgb_gamma = cmdargs_xgb_gamma # 0 [0 - 2] the larger, the more conservative the algorithm will be.
    xgb_early_stopping_rounds = 300 # stop if performance keeps getting worse consecutively for k rounds.
    xgb_num_parallel_tree = 2000 # Random Forests with XGBoost (si xgb_nround == 1 & xgb_subsample < 1 & xgb_colsample_bytree < 1, claro)
    if(xgb_nround > 1 | xgb_subsample == 1 | xgb_colsample_bytree == 1)
    {
      xgb_s_tipo_modelo = xgb_tipo_modelo_XGB   # "XGB" XGBoost
      xgb_num_parallel_tree = cmdargs_xgb_num_parallel_tree # 1 [1 - inf] Experimental. Helps reducing over-fitting. Prueba que, con pocos datos, parece mejorar (aunque tarda mucho más)
    } else
    {
      xgb_s_tipo_modelo = xgb_tipo_modelo_XGRF # "XGRF" Random Forest con XGBoost
      xgb_nround = 1 # "XGRF" Random Forest con XGBoost
      if(G_b_DEBUG & xgb_num_parallel_tree > 1)  xgb_num_parallel_tree = 100 # "XGRF" Random Forest con XGBoost
    }
    # ===========================  XGB params - END  ==============================
    
    if(G_b_DEBUG)  xgb_mi_version <- xgb_mi_version + 990000
    # Nombre del fichero para guardar estadísticas (de estos modelos) en un único fichero:
    # xgb_filename <- str_replace(xgb_modelos[numModelo], "_[0-9]*\\.[0-9][0-9][0-9]$", "")
    xgb_filename <- paste0(paste(xgb_s_tipo_modelo, G_CV_nFolds,
                                 str_replace(xgb_objective, "^(...).*:(...).*$", "\\1\\2"),
                                 xgb_eta, xgb_max_depth, xgb_subsample, xgb_colsample_bytree, xgb_eval_metric, xgb_early_stopping_rounds, xgb_scale_pos_weight, xgb_min_child_weight, xgb_gamma, xgb_num_parallel_tree,
                                 jjfmt(100*xgb_train_porc, 1, '.', ','), # jjfmt(100*nrow(trainset)/(nrow(trainset)+nrow(validset)), 1, '.', ','),
                                 paste0('v', str_pad(xgb_mi_version, 3, "left", "0")),
                                 sep = "_"))
    # Iniciamos log (xgb_filename".log"):
    if(G_b_NO_LOG == FALSE)
    {
      s_Fic_log <- paste0(s_output_path, xgb_filename, '.log')
      if(file.exists(s_Fic_log))
        jjprint('-------------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE) # Línea separadora para el inicio de este log...
    }
    jjprint(paste0('trainset: ', jjfmt(nrow(trainset)), ' registros.  validset: ', jjfmt(nrow(validset)), ' registros. [', Sys.info()["nodename"], ']'), logfile = s_Fic_log, b_add_datetime = TRUE)
    jjprint(paste0('Preparando datos ', str_pad(numModelo, 3, "left", "0"), '...[', xgb_filename, ']'), logfile = s_Fic_log, b_add_datetime = TRUE)
    # Verificamos que no exista el modelo ya entrenado (fichero):
    if(length(dir(path = s_modelos_path, pattern = paste0(xgb_filename, '.*', '.', str_pad(numModelo, 3, "left", "0"), '.modelo'))) != 0)
    {
      xgb_modelos[numModelo] <- dir(path = s_modelos_path, pattern = paste0(xgb_filename, '.*', '.', str_pad(numModelo, 3, "left", "0"), '.modelo'))[1]
      xgb_modelos[numModelo] <- str_replace(xgb_modelos[numModelo], pattern = "\\.modelo", replacement = "")
      jjprint(paste0('Nota: Encontrado modelo ya entrenado. Lo cargamos ', str_pad(numModelo, 3, "left", "0"), ' (', xgb_modelos[numModelo], ') y pasamos al siguiente...'), logfile = s_Fic_log)
      next # el "continue" de C para el bucle for(forNumModelo...
    }
    
    xgb_preps <- xgb_prep_datos(mi_dt = trainset, b_verbose = 1, maxImportanceNumVars = maxImportanceNumVars, s_modelos_path = s_modelos_path
                                , G_maxNumFeatures_OneHotEnc = G_maxNumFeatures_OneHotEnc
                                , obj_var_name = 'Class', s_Fic_log = s_Fic_log)
    # Seleccionamos variables (si hay otros modelos previos del mismo numAd y misma versión y con importance_matrix):
    # importance_vars <- xgb_prep_datos_busca_imp_matrix(xgb_filename) # vector de variables ordenadas por importancia (las primeras son las más importantes)
    # xgb_preps[[2]] <- xgb_prep_datos_selec_vars(mi_dt = trainset, numVars = 50, importance_vars, b_verbose = 1)
    jjprint(colnames(xgb_preps[[2]]), logfile = s_Fic_log, b_add_datetime = FALSE) # Finalmente, mostramos las columnas elegidas
    
    # sapply(xgb_preps[[2]], uniqueN)[sapply(xgb_preps[[2]], uniqueN)==1]
    X <- data.matrix(xgb_preps[[2]][,-1, with=FALSE])
    y <- xgb_preps[[2]]$Class - 1 # NOTA para multi-class: Class is represented by a number and should be from 0 to num_class - 1.
    # midata <- xgb.DMatrix(data.matrix(xgb_preps[[2]][,-1, with=FALSE]), label = xgb_preps[[2]]$clicked, missing = NA)
    # Ponderación: Pesos (de cada fila) en función de su clase:
    mis_weights <- NULL
    if(b_conPond == 1) {
      mis_weights <- as.numeric(xgb_weights_clases)[(y+1)]    # Frecs de cada clase (en trainset)
    } else if(b_conPond == 2) {
      mis_weights <- as.numeric(xgb_weights_clases200)[(y+1)] # 1 / Frecs de cada clase (en trainset)
    }
    # Reducimos memoria:
    if(!G_b_DEBUG)  rm(xgb_preps); gc()
    s_misDims <- paste0(jjfmt(dim(X)[1]), ' regs/',
                        jjfmt(dim(X)[2]), ' cols')
    # Entrenamos:
    jjprint(paste0('Entrenando (CV) ', str_pad(numModelo, 3, "left", "0"), '...[', s_misDims, '] [', xgb_filename, ']'), logfile = s_Fic_log)
    
    mi_tiempo <- system.time({
    xgb_cv <- xgb.cv(data = X, label = y, missing = NA,
                     nfold = G_CV_nFolds, prediction = xgb_get_predictions, # SOLAMENTE PARA CV
                     nround = xgb_nround,
                     num_class = xgb_num_class,
                     objective = xgb_objective, eta = xgb_eta, max_depth = xgb_max_depth, subsample = xgb_subsample, colsample_bytree = xgb_colsample_bytree,
                     metrics = list(xgb_eval_metric), # eval_metric = xgb_eval_metric,
                     early_stopping_rounds = xgb_early_stopping_rounds,
                     weight = mis_weights, scale_pos_weight = xgb_scale_pos_weight,
                     min_child_weight = xgb_min_child_weight, gamma = xgb_gamma,
                     num_parallel_tree = xgb_num_parallel_tree,
                     nthread = (detectCores()-1),
                     verbose = T,
                     save_period = NULL
                    )
    })
    tiempo_cv <- paste0('Tiempo xgb_cv(): ', mi_tiempo['elapsed']/60, ' minutos')
    jjprint(paste0(tiempo_cv, ' [', s_misDims, ']'), logfile = s_Fic_log)
    if(xgb_get_predictions)
    {
      xgb_preds <- xgb_cv$pred # Predicciones en los N-Folds no usados para entrenar (vector de tamaño nrow(trainset))
      xgb_cv <- xgb_cv$dt      # Este es el mismo data.table que devolvía xgb.cv con prediction=FALSE
    }
    if(xgb_eval_metric == "map")
    {
      xgb_nround <- which.max(xgb_cv$evaluation_log$test_map_mean)
      mi_xgb_cv_train_bestScore <- max(xgb_cv$evaluation_log$train_map_mean)
      mi_xgb_cv_test_bestScore <- max(xgb_cv$evaluation_log$test_map_mean)
    } else if(xgb_eval_metric == "mlogloss") {
      xgb_nround <- which.min(xgb_cv$evaluation_log$test_mlogloss_mean)
      mi_xgb_cv_train_bestScore <- min(xgb_cv$evaluation_log$train_mlogloss_mean)
      mi_xgb_cv_test_bestScore <- min(xgb_cv$evaluation_log$test_mlogloss_mean)
    } else if(xgb_eval_metric == "merror") {
      xgb_nround <- which.min(xgb_cv$evaluation_log$test_merror_mean)
      mi_xgb_cv_train_bestScore <- min(xgb_cv$evaluation_log$train_merror_mean)
      mi_xgb_cv_test_bestScore <- min(xgb_cv$evaluation_log$test_merror_mean)
    } else if(xgb_eval_metric == "error") {
      xgb_nround <- which.min(xgb_cv$evaluation_log$test_error_mean)
      mi_xgb_cv_train_bestScore <- min(xgb_cv$evaluation_log$train_error_mean)
      mi_xgb_cv_test_bestScore <- min(xgb_cv$evaluation_log$test_error_mean)
    } else if(xgb_eval_metric == "auc") {
      xgb_nround <- which.max(xgb_cv$evaluation_log$test_auc_mean)
      mi_xgb_cv_train_bestScore <- max(xgb_cv$evaluation_log$train_auc_mean)
      mi_xgb_cv_test_bestScore <- max(xgb_cv$evaluation_log$test_auc_mean)
    } else {
      save(xgb_cv, file = paste0(s_output_path, "xgb_cv_temp.RData"))
      stop('eval_metric desconocida (Cf. xgb_cv_temp.RData!')
    }
    jjprint(xgb_cv$evaluation_log[1:min((xgb_nround+2),nrow(xgb_cv$evaluation_log))], logfile = s_Fic_log, b_add_datetime = FALSE) # +2 para ver los 2 sgtes, si los hay...
    jjprint(paste0('Best Test eval Nround = ', xgb_nround, ' [', xgb_filename, ']'), logfile = s_Fic_log)
    mi_tiempo <- system.time({
      xgb <- xgboost(data = X, label = y, missing = NA,
                     nround = xgb_nround,
                     num_class = xgb_num_class,
                     objective = xgb_objective, eta = xgb_eta, max_depth = xgb_max_depth, subsample = xgb_subsample, colsample_bytree = xgb_colsample_bytree,
                     eval_metric = xgb_eval_metric, # metrics es para xgb.cv()
                     early_stopping_rounds = xgb_nround, # Para calcular mi_xgb_train_score (i.e. con xgb$evaluation_log)
                     weight = mis_weights, scale_pos_weight = xgb_scale_pos_weight,
                     min_child_weight = xgb_min_child_weight, gamma = xgb_gamma,
                     num_parallel_tree = xgb_num_parallel_tree,
                     nthread = (detectCores()),
                     verbose = 1, # 0,1,2
                     save_period = NULL
                     )
    })
    tiempo_train <- paste0('Tiempo xgb(): ', mi_tiempo['elapsed']/60, ' minutos')
    jjprint(paste0(tiempo_train, ' [', s_misDims, ']'), logfile = s_Fic_log)
    
    # Añadimos nround obtenido en la CV (y el numModelo) al nombre de fichero de cada modelo:
    xgb_modelos[numModelo] <- paste0(paste(xgb_filename,
                                           xgb_nround,
                                        sep = "_")
                                     , '.', str_pad(numModelo, 3, "left", "0")) # lo guardamos para poder hacer ensemble luego
    jjprint(paste0('Predecimos en trainset para medir nuestro score(', xgb_eval_metric, '):'), logfile = s_Fic_log, b_add_datetime = TRUE)
    mi_tiempo <- system.time({
      preds <- predict(xgb, newdata = X, missing = NA, reshape = (xgb_num_class > 2))
      if(xgb_num_class == 2){
        trainset[, prob := preds]
      } else {
        for(i in 1:xgb_num_class)
          trainset[, paste0("prob_class", i) := preds[,i]]
      }
    })
    tiempo_predict <- paste0('Tiempo predict(xgb, trainset): ', mi_tiempo['elapsed'], ' segundos')
    jjprint(paste0(tiempo_predict, ' [', s_misDims, ']'), logfile = s_Fic_log, b_add_datetime = TRUE)
    # # # Tabla de confusión:
    # # table(trainset[, .(clicked, prob>.5)])
    mi_xgb_train_score <- as.numeric(xgb$evaluation_log[xgb_nround, 2])
    jjprint(paste0('xgb_train_score(', xgb_eval_metric,')(xgboost) = ', mi_xgb_train_score, '. Cluster = ', numCluster), logfile = s_Fic_log, b_add_datetime = TRUE)
    mi_xgb_train_score <- mi_mlogloss(tr_valid = trainset, b_verbose = 0)
    jjprint(paste0('xgb_train_score(', xgb_eval_metric,')(mi_', xgb_eval_metric,') = ', mi_xgb_train_score, '. Cluster = ', numCluster), logfile = s_Fic_log, b_add_datetime = TRUE)
    # mi_map12 <- Map_12(tr_valid = trainset[, .(display_id, ad_id, clicked, prob)], b_verbose = 0)
    # jjprint(paste0('map12 = ', mi_map12), logfile = s_Fic_log, b_add_datetime = TRUE)
    jjprint(paste0('Guardando modelo entrenado (', xgb_modelos[numModelo] ,'.modelo', ')'), logfile = s_Fic_log, b_add_datetime = TRUE)
    if(!G_b_DEBUG)
      xgboost::xgb.save(model = xgb, fname = paste0(s_modelos_path, xgb_modelos[numModelo] ,'.modelo'))
    # # Leemos modelo:
    # xgb <- xgb.load(paste0(s_modelos_path, str_replace(xgb_modelos[numModelo], pattern = "\\.modelo", replacement = ""), '.modelo'))
    write.table(x = t(as.data.frame(list(modelo = paste(numModelo, NUM_MODELOS, sep = "/"),
                        dim_x = t(dim(X)),
                        tiempo_cv = tiempo_cv,
                        tiempo_train = tiempo_train,
                        tiempo_predict = tiempo_predict,
                        cols = t(dimnames(X)[[2]]),
                        list(xgb_eta=xgb_eta, xgb_max_depth=xgb_max_depth, xgb_subsample=xgb_subsample, xgb_colsample_bytree=xgb_colsample_bytree,
                             xgb_eval_metric=xgb_eval_metric,xgb_early_stopping_rounds=xgb_early_stopping_rounds,
                             xgb_min_child_weight=xgb_min_child_weight, xgb_gamma=xgb_gamma, xgb_num_parallel_tree=xgb_num_parallel_tree),
                        xgb_train_pond_clases = ifelse(b_conPond == 1, t(xgb_weights_clases), ifelse(b_conPond == 2, t(xgb_weights_clases200), 'No')),
                        xgb_train_best_score = mi_xgb_train_score,
                        reduc_seed = xgb_reduc_seed,
                        reduc_train_porcent = xgb_train_porc,
                        xgb_cv_bestRound = xgb_nround,
                        xgb_cv_train_bestScore = mi_xgb_cv_train_bestScore,
                        xgb_cv_test_bestScore = mi_xgb_cv_test_bestScore))),
              file = paste0(s_output_path, xgb_modelos[numModelo] ,'.txt'), row.names = T, col.names = F, quote = F, sep = " = ")
  
    suppressWarnings( # "Warning: appending column names to file"
      write.table(x = as.data.frame(list(Modelo = xgb_modelos[numModelo])), file = paste0(s_modelos_path, xgb_filename ,'.txt'), row.names = T, col.names = F, quote = F, sep = " = ", append = TRUE)
    )
    
    # Como la importance_matrix tarda bastante en calcularse, no lo hacemos para todos los numModelo:
    # if(numModelo %% 3 == 1) # 1, 4, 7, 10
    # {
    if(Sys.info()["nodename"]=="JJTZAPATA-W10"){
      jjprint('NOTA: No se calcula la feature matrix', logfile = s_Fic_log, b_add_datetime = TRUE)
      suppressWarnings( # "Warning: appending column names to file"
        write.table(x = as.data.frame(list(Nota = 'No se calcula la feature matrix')), file = paste0(s_modelos_path, xgb_filename ,'.txt'), row.names = T, col.names = F, quote = F, sep = " = ", append = TRUE)
      )
    } else {
      jjprint('Compute feature importance matrix...', logfile = s_Fic_log, b_add_datetime = TRUE)
      mi_tiempo <- system.time({
        importance_matrix <- xgb.importance(feature_names = dimnames(X)[[2]], model = xgb)
      })
      tiempo_imp_matrix <- paste0('Tiempo xgb.importance(): ', mi_tiempo['elapsed'], ' segundos')
      jjprint(paste0(tiempo_imp_matrix, ' [', s_misDims, ']'), logfile = s_Fic_log, b_add_datetime = TRUE)
      suppressWarnings(
        # Gain:
        write.table(x = as.data.frame(importance_matrix)[,c(1,2)], file = paste0(s_modelos_path, xgb_filename ,'.txt'), row.names = T, col.names = T, quote = F, sep = "\t", append = TRUE)
      )
      suppressWarnings(
        # Cover:
        write.table(x = as.data.frame(importance_matrix)[,c(1,3)], file = paste0(s_modelos_path, xgb_filename ,'.txt'), row.names = T, col.names = T, quote = F, sep = "\t", append = TRUE)
      )
      suppressWarnings(
        # Frequency:
        write.table(x = as.data.frame(importance_matrix)[,c(1,4)], file = paste0(s_modelos_path, xgb_filename ,'.txt'), row.names = T, col.names = T, quote = F, sep = "\t", append = TRUE)
      )
      # if(G_b_DEBUG) # xgb.plot.importance(), barplot(), chisq.test()...
      # {
      #   jjprint(dimnames(X)[[2]], logfile = NULL, b_add_datetime = FALSE)
      #   # # Nice graph
      #   # # install.packages("Ckmeans.1d.dp")
      #   # library(Ckmeans.1d.dp)
      #   # x11(); xgb.plot.importance(importance_matrix)
      #   
      #   #In case last step does not work for you because of a version issue, you can try following :
      #   x11(); barplot(importance_matrix[1:8]$Gain, names.arg = importance_matrix[1:8]$Feature,
      #                  main = paste0("XGB varImp.Gain - ", str_pad(numModelo, 3, "left", "0")))
      #   # x11(); barplot(importance_matrix[1:8]$Cover, names.arg = importance_matrix[1:8]$Feature,
      #   #                main = paste0("XGB varImp.Cover - ", str_pad(numModelo, 3, "left", "0")))
      #   # x11(); barplot(importance_matrix[1:8]$Frequence, names.arg = importance_matrix[1:8]$Feature,
      #   #                main = paste0("XGB varImp.Frequence - ", str_pad(numModelo, 3, "left", "0")))
      #   
      #   # # install.packages("DiagrammeR")
      #   # library(DiagrammeR)
      #   # x11(); xgb.plot.tree(feature_names = dimnames(X)[[2]], model = xgb, n_first_tree = 3)
      #   
      #   # # To see whether the variable is actually important or not:
      #   # test <- chisq.test(y, as.numeric(trainset$ad_publish_timestamp))
      #   # test <- chisq.test(y, as.numeric(trainset$publish_timestamp))
      #   # test <- chisq.test(y, as.numeric(trainset$numAds))
      #   # jjprint(test, logfile = NULL, b_add_datetime = FALSE)
      # }
    }
    
    if(nrow(validset) != 0)
    {
      jjprint('-------------------------------', logfile = s_Fic_log, b_add_datetime = TRUE)
      jjprint('Ahora predecimos en validset, que son los que NO hemos usado para CV:', logfile = s_Fic_log, b_add_datetime = TRUE)
      jjprint('-------------------------------', logfile = s_Fic_log, b_add_datetime = TRUE)
      jjprint(paste0('trainset: ', jjfmt(nrow(trainset)), ' registros.  validset: ', jjfmt(nrow(validset)), ' registros.'), logfile = s_Fic_log, b_add_datetime = TRUE)
      xgb_preps <- xgb_prep_datos(mi_dt = validset, b_verbose = 1, maxImportanceNumVars = maxImportanceNumVars, s_modelos_path = s_modelos_path, G_maxNumFeatures_OneHotEnc = G_maxNumFeatures_OneHotEnc, obj_var_name = 'Class', s_Fic_log = s_Fic_log)
      # Seleccionamos variables (si hay otros modelos previos del mismo numAd y misma versión y con importance_matrix):
      # importance_vars <- xgb_prep_datos_busca_imp_matrix(xgb_filename) # vector de variables ordenadas por importancia (las primeras son las más importantes)
      # xgb_preps[[2]] <- xgb_prep_datos_selec_vars(mi_dt = trainset, numVars = 50, importance_vars, b_verbose = 1)
      # jjprint(colnames(xgb_preps[[2]]), logfile = NULL, b_add_datetime = FALSE)
      # sapply(xgb_preps[[2]], uniqueN)
      X <- data.matrix(xgb_preps[[2]][,-1, with=FALSE])
      y <- xgb_preps[[2]]$Class - 1 # NOTA para multi-class: Class is represented by a number and should be from 0 to num_class - 1.
      # midata <- xgb.DMatrix(data.matrix(xgb_preps[[2]][,-1, with=FALSE]), label = xgb_preps[[2]]$clicked, missing = NA)
      rm(xgb_preps); gc()
      jjprint(paste0('Predecimos en validset para medir nuestro score (',xgb_eval_metric ,'):'), logfile = s_Fic_log, b_add_datetime = TRUE)
      mi_tiempo <- system.time({
        preds <- predict(xgb, newdata = X, missing = NA, reshape = (xgb_num_class > 2))
        if(xgb_num_class == 2){
          validset[, prob := preds]
        } else {
          for(i in 1:xgb_num_class)
            validset[, paste0("prob_class", i) := preds[,i]]
        }
      })
      tiempo_predict <- paste0('Tiempo predict(xgb, validset): ', mi_tiempo['elapsed'], ' segundos')
      jjprint(tiempo_predict, logfile = s_Fic_log, b_add_datetime = TRUE)
      # # Tabla de confusión:
      # table(validset[, .(clicked, prob>.5)])
      # mi_map12_valid <- Map_12_diario(tr_valid = validset[, .(display_id, ad_id, dia, clicked, prob)], b_verbose = 0,
      #                           dias_list = list(c(1:3), c(3:5), c(5:7), c(7:9), c(9:11), c(11:13), c(1:11), c(12:13), c(1:13)))
      # colnames(mi_map12_valid) <- c('map12val_1-3', 'map12val_3-5', 'map12val_5-7', 'map12val_7-9', 'map12val_9-11', 'map12val_11-13', 'map12val_1-11', 'map12val_12-13', 'map12val_Total')
      # map12val_Total <- mi_map12_valid[2, "map12val_Total"][[1]]
      # p_esp <- 1 / numCluster
      # mejora <- 100 * (map12val_Total - p_esp) * (1 - p_esp)
      # stats2 <- data.frame(numModelo = str_pad(numModelo, 3, "left", "0"), numAds = numCluster, map12val_Total = map12val_Total, mejora = mejora)
      mi_xgb_valid_score <- mi_mlogloss(tr_valid = validset, b_verbose = 1)
    } else {
      jjprint('-------------------------------', logfile = s_Fic_log, b_add_datetime = TRUE)
      jjprint('NO predecimos en validset, porque NO hay validset(hemos usado todo para CV):', logfile = s_Fic_log, b_add_datetime = TRUE)
      jjprint('-------------------------------', logfile = s_Fic_log, b_add_datetime = TRUE)
      mi_xgb_valid_score <- 0
    }
    jjprint(paste0('xgb_train_score(', xgb_eval_metric,')(mi_', xgb_eval_metric,') = ', mi_xgb_train_score, '. Cluster = ', numCluster), logfile = s_Fic_log, b_add_datetime = TRUE)
    jjprint(paste0('xgb_valid_score(', xgb_eval_metric,')(mi_', xgb_eval_metric,') = ', mi_xgb_valid_score, '. Cluster = ', numCluster), logfile = s_Fic_log, b_add_datetime = TRUE)
    jjprint('----------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
    jjprint(data.frame(score_name = xgb_eval_metric, valid_score = mi_xgb_valid_score, train_score = as.numeric(mi_xgb_train_score), CV_test_score = mi_xgb_cv_test_bestScore), logfile = s_Fic_log, b_add_datetime = FALSE) # Volvemos a imprimir mi_xgb_train_score (del Training) al lado.
    # jjprint(cbind(mi_map12_valid[2,], map12Train = mi_xgb_train_score, map12cvTest = mi_xgb_cv_test_bestScore), logfile = s_Fic_log, b_add_datetime = FALSE) # Volvemos a imprimir mi_xgb_train_score (del Training) al lado.
    # jjprint(stats2, logfile = s_Fic_log, b_add_datetime = FALSE)
    jjprint('----------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
    # Añadimos a los ficheros de cada training los scores del validset:
    write.table(x = c(xgb_valid_score = mi_xgb_valid_score),
                file = paste0(s_output_path, xgb_modelos[numModelo] ,'.txt'), row.names = T, col.names = F, quote = F, sep = " = ", append = TRUE)
    # # Añadimos a los ficheros de cada training un resumen con la mejora (porcentaje sobre el máximo posible teórico):
    # write.table(x = t(stats2[,-3]),
    #             file = paste0(s_output_path, xgb_modelos[numModelo] ,'.txt'), row.names = T, col.names = F, quote = F, sep = " = ", append = TRUE)
    
    minutos <- as.double((proc.time() - systime_ini)['elapsed'])/60
    if(minutos < 60) strtmp <- paste0(minutos, ' minutos en total.', ' - ', 'trainset: ', jjfmt(nrow(trainset)), ' registros.  validset: ', jjfmt(nrow(validset)), ' registros.') else  strtmp <- paste0(minutos/60, ' horas en total.', ' - ', 'trainset: ', jjfmt(nrow(trainset)), ' registros.  validset: ', jjfmt(nrow(validset)), ' registros.')
    jjprint(strtmp, logfile = s_Fic_log, b_add_datetime = TRUE)
    minutos_pend <- (as.double((proc.time() - systime_ini)['elapsed'])/60) * (NUM_MODELOS / numModelo  -  1)
    if(minutos_pend < 60) strtmp <- paste0('Faltan aprox. ', minutos_pend, ' minutos.') else strtmp <- paste0('Faltan aprox. ', minutos_pend/60, ' horas.')
    jjprint(strtmp, logfile = s_Fic_log, b_add_datetime = TRUE)
  }
}

# length(xgb_modelos) <- 1 # Si tenemos menos modelos, podemos crear submit con menos...
if(is.null(s_Fic_log))
{
  # Por si acaso no hemos entrado a entrenar...
  if(G_b_NO_LOG == FALSE)
  {
    s_Fic_log <- paste0(s_output_path, xgb_filename, '.log')
    if(file.exists(s_Fic_log))
      jjprint('-------------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE) # Línea separadora para el inicio de este log...
  }
}
jjprint('-------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
jjprint('Predecimos en testset y creamos submit:', logfile = s_Fic_log, b_add_datetime = TRUE)
jjprint('-------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
if(G_b_DEBUG || !G_b_crearSubmit)
{
  if(G_b_DEBUG)
    jjprint('NOTA: G_b_DEBUG == TRUE -> No creamos submit.csv', logfile = s_Fic_log, b_add_datetime = TRUE)
  if(!G_b_crearSubmit)
    jjprint('NOTA: G_b_crearSubmit == FALSE -> No creamos submit.csv', logfile = s_Fic_log, b_add_datetime = TRUE)
} else {
  if(!exists("xgb_modelos"))  xgb_modelos <- vector(mode = "character", length = NUM_MODELOS)
  predict_testset(nombres_modelos = xgb_modelos
                  , filename = xgb_filename, s_input_path = s_input_path, s_output_path = s_output_path, s_modelos_path = s_modelos_path
                  , i_sDescr = "XGB Predicting"
                  , FUN_prep_datos = xgb_prep_datos, prep_datos_b_verbose = 1
                  , obj_var_name = 'Class', num_classes = xgb_num_class, G_maxNumFeatures_OneHotEnc = G_maxNumFeatures_OneHotEnc
                  , FUN_X = data.matrix
                  , FUN_Predict = function(modelo, X){ return(predict(modelo, X, missing = NA, reshape = (xgb_num_class > 2))) }
                  , FUN_loadmodelo = xgb.load
                  , nombreModeloLoaded = ""
                  , NUM_MODELOS = length(xgb_modelos)
                  , b_ReducirFicheros = FALSE # De momento no hace falta reducir más los testset, pero podría hacer falta...
                  , modelos_weights = NULL
  )
}
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# M12 <- foreach(reg = c(0,3,4,5,6), .inorder=TRUE, .combine = c, .packages = c('data.table'),
#                .export = c('Map_12_diario', 'Map_12', 'mapk_nodup_only_first_12', 'add_tiempo')
#                ) %do%
# {
#   return(
#     basic_preds_m12(trainset = trainset, validset = validset, k = 2, reg = reg, b_verbose = 1)
#     ) # ret_val de la función "foreach() %do%" (o foreach( %dopar%))
# }
# # M12    <- basic_preds_m12(trainset = trainset, validset = validset, k = 1)
# # M12[2] <- basic_preds_m12(trainset = trainset, validset = validset, k = 2)
# # M12[3] <- basic_preds_m12(trainset = trainset, validset = validset, k = 3)
# print(M12)
# # sort(M12)
# sort(unlist(sapply(M12, function(x) {return(x$V1[2])})))
# sort(unlist(sapply(M12, function(x) {return(x$V2[2])})))
# sort(unlist(sapply(M12, function(x) {return(x$V3[2])})))

tot_mins <- as.double((proc.time() - systime_ini)['elapsed'])/60
if(tot_mins < 60){
  jjprint(paste0(round(tot_mins,3), ' minutos en total.'), logfile = s_Fic_log, b_add_datetime = TRUE)
} else if(tot_mins < 60 * 24) {
  jjprint(paste0(round(tot_mins/60,3), ' horas en total.'), logfile = s_Fic_log, b_add_datetime = TRUE)
} else {
  jjprint(paste0(round(tot_mins/(60*24),3), ' días en total.'), logfile = s_Fic_log, b_add_datetime = TRUE)
}

jjprint('Ok.', logfile = s_Fic_log, b_add_datetime = TRUE)

# cleanup:
try(registerDoSEQ(), silent = TRUE) # library(doParallel) [turn parallel processing off and run sequentially again]
try(stopImplicitCluster(), silent = TRUE)
try(stopCluster(cl), silent = TRUE)
