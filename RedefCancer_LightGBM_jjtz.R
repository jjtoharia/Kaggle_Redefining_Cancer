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

# # Inicializamos variables:
# # NOTA: Dejamos un Core de la CPU "libre" para no "quemar" la máquina:
# cl <- makeCluster(detectCores(), type='PSOCK') # library(doParallel) [turn parallel processing on]
# registerDoParallel(cl) # library(doParallel) [turn parallel processing on]

# memory.limit(size = 16000)

systime_ini <- proc.time()

cmdargs = commandArgs(trailingOnly=TRUE)
if(length(cmdargs) >= 1) print(paste0('cmdargs = [', paste(cmdargs, collapse = ' '),']'))
G_b_NO_LOG <- any(cmdargs == "/q"); cmdargs <- cmdargs[cmdargs != "/q"]
# 1 - b_crearSubmit           # 0 (0, 1)
# 2 - max_depth = -1          # -1 [1 - inf] maximum depth of the tree (demasiado alto => overfitting) [max_depth = -1 for infinite depth]
# 3 - num_leaves = 511        # 31 [1 - Inf] Typical: 255, usually {15, 31, 63, 127, 255, 511, 1023, 2047, 4095}. Tips: adjust depth accordingly by allowing a slightly higher depth than the theoretical number of leaves.
# 4 - feature_fraction        # 1 (0 - 1] [colsample_bytree en xgboost] (muestreo aleatorio por columnas si es < 1, para intentar evitar overfitting)
# 5 - min_sum_hessian_in_leaf # 10 [0 - Inf] Typical = 1. [min_child_weight en xgboost] # Smaller/Larger is not always better. Leave it alone unless you know what you are doing.
# 6 - min_gain_to_split       # 0 [0 - 2] [gamma en xgboost] the larger, the more conservative (stable, slow-learner) the algorithm will be.
# 7 - min_data_in_leaf        # 20 [1 - 10k] the larger, the more conservative (less over-fitting) the algorithm will be.[desde V.x20]
cmdargs_lgb_max_depth = ifelse(length(cmdargs) >= 2, as.numeric(cmdargs[2]), -1) # -1 # maximum depth of trees (demasiado alto => overfitting) [max_depth = -1 for infinite depth]
cmdargs_lgb_num_leaves = ifelse(length(cmdargs) >= 3, as.numeric(cmdargs[3]), 31) # 31 [1 - Inf] Typical: 255, usually {15, 31, 63, 127, 255, 511, 1023, 2047, 4095}. Tips: adjust depth accordingly by allowing a slightly higher depth than the theoretical number of leaves.
cmdargs_lgb_feature_fraction = ifelse(length(cmdargs) >= 4, as.numeric(cmdargs[4]), 1) # [colsample_bytree en xgboost] (muestreo aleatorio por columnas si es < 1, para intentar evitar overfitting)
cmdargs_lgb_min_sum_hessian_in_leaf = ifelse(length(cmdargs) >= 5, as.numeric(cmdargs[5]), 10) # [min_child_weight en xgboost] Typical = 1. Default = 10 [0 - Inf] # Smaller/Larger is not always better. Leave it alone unless you know what you are doing.
cmdargs_lgb_min_gain_to_split = ifelse(length(cmdargs) >= 6, as.numeric(cmdargs[6]), 0) # [gamma en xgboost] 0 [0 - 2] the larger, the more conservative (less over-fitting) the algorithm will be.
cmdargs_lgb_min_data_in_leaf = ifelse(length(cmdargs) >= 7, as.numeric(cmdargs[7]), 20) # 20 [1 - 10k] the larger, the more conservative (less over-fitting) the algorithm will be.[desde V.x20]
if(cmdargs_lgb_min_data_in_leaf == -1) { # v.114
   lgb_mascara_filename <- paste0(paste('.*',
                                ifelse(length(cmdargs) >= 2, (if(cmdargs_lgb_max_depth==-1) 0 else cmdargs_lgb_max_depth), '.*'),
                                ifelse(length(cmdargs) >= 3, cmdargs_lgb_num_leaves, '.*'),
                                '.*',
                                ifelse(length(cmdargs) >= 4, cmdargs_lgb_feature_fraction, '.*'),
                                '.*',
                                ifelse(length(cmdargs) >= 5, cmdargs_lgb_min_sum_hessian_in_leaf, '.*'),
                                ifelse(length(cmdargs) >= 6, cmdargs_lgb_min_gain_to_split, '.*'),
                                '.*',
                                sep = "_"))
} else { # v.134
   lgb_mascara_filename <- paste0(paste('.*',
                                ifelse(length(cmdargs) >= 2, (if(cmdargs_lgb_max_depth==-1) 0 else cmdargs_lgb_max_depth), '.*'),
                                ifelse(length(cmdargs) >= 3, cmdargs_lgb_num_leaves, '.*'),
                                '.*',
                                ifelse(length(cmdargs) >= 4, cmdargs_lgb_feature_fraction, '.*'),
                                '.*',
                                ifelse(length(cmdargs) >= 5, cmdargs_lgb_min_sum_hessian_in_leaf, '.*'),
                                ifelse(length(cmdargs) >= 6, cmdargs_lgb_min_gain_to_split, '.*'),
                                ifelse(length(cmdargs) >= 7, cmdargs_lgb_min_data_in_leaf, '.*'),
                                '.*',
                                sep = "_"))
}
# --------------------------------------------------------
G_b_DEBUG <- FALSE # Reducimos todo para hacer pruebas más rápido
G_b_REV <- FALSE # Empezamos por el final (numModelo <- NUM_MODELOS - forNumModelo + 1)
G_maxNumFeatures_OneHotEnc <- 10000 # One-Hot Encode las categóricas con estos levels como máx.
G_CV_nFolds <- 20 # 5 [2 - Inf]
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
# # install.packages("devtools")
# library(devtools)
# options(devtools.install.args = "--no-multiarch") # if you have 64-bit R only, you can skip this
# install_github("Microsoft/LightGBM", subdir = "R-package")
# Si lo anterior no funciona (Windows), instalar desde ../Tmp_Kaggle/RedefiningCancer/lightgbm
library(R6) # required by lightgbm
library(lightgbm)
# # Check lightgbm installation:
# data(agaricus.train, package='lightgbm')
# train <- agaricus.train
# dtrain <- lgb.Dataset(train$data, label=train$label)
# params <- list(objective="regression", metric="l2")
# model <- lgb.cv(params, dtrain, 100, nfold=5, min_data=1, learning_rate=0.1, early_stopping_rounds=10
#                 , stratified = TRUE, num_leaves = 127, verbose = 1
#                 , max_depth = 15, num_threads = (detectCores()+0.5)%/%2) # (detectCores()+0.5)%/%2 debería ser el número real de "cores"
G_maxNumCores <- if((detectCores()+0.5)%/%2 > 1) (detectCores()+0.5)%/%2 - 1 else 1 # Dejamos uno libre para no quemar la máquina...
# G_maxNumCores <- max(1, G_maxNumCores - 1)
# -------------------------------
# Inicializamos:
# -------------------------------
lgb_scale_pos_weight <- 1 # default
# # # mi_sum_pos <- sum(full_trainset$clicked == 1) # 16874593
# # # mi_sum_neg <- sum(full_trainset$clicked == 0) # 70267138
# # mi_sum_pos <- 16874593
# # mi_sum_neg <- 70267138
# lgb_scale_pos_weight <- mi_sum_neg / mi_sum_pos # == 4.164079
# -------------------------------
b_conPond <- 0 # Sin ponderar (V.100)
# Ponderación cancelada: afecta a la métrica de cálculo del error, i.e. mlogloss
# y no parece mejorar en nada el entrenamiento.
# En RedefCancer explícitamente avisan que no se tenga en cuenta la frecuencia de cada clase, así que razón de más...
# A estudiar en el futuro...
# # full_trainset <- fread(file.path(s_input_path, "training_variants"))
# # weights_clases <- table(full_trainset$Class)/nrow(full_trainset)
# # paste0('lgb_weights_clases <- c(', paste0(paste(paste0("'",names(weights_clases),"'"), weights_clases, sep = '='), collapse = ','), ')')
# # lgb_weights_clases <- c('1'=0.171032821439325,'2'=0.136103583258055,'3'=0.0267991568804577,'4'=0.206564287865101,'5'=0.0728696175850647,'6'=0.0828063836193917,'7'=0.286961758506474,'8'=0.00572116832279434,'9'=0.0111412225233363)
# # b_conPond <- 1 # V.200
# # # lgb_weights_clases200 <- 1 / lgb_weights_clases
# # # lgb_weights_clases200 <- lgb_weights_clases200 / sum(lgb_weights_clases200)
# # # b_conPond <- 2 # V.300
# -------------------------------
lgb_num_class = 9 # lgb_objective = "multiclass"
lgb_reduc_seed <- 1234
lgb_train_porc <- 1 # 1 # 0.99 # 0.85
# lgb_get_predictions <- FALSE # TRUE para hacer stacking de los modelos (Subsemble, SuperLearner)

if(G_b_crearSubmit) {
  print('-------------------------------')
  print(paste0('Primero buscamos algún modelo [', lgb_mascara_filename, '] entrenado pero SIN submit:'))
  print('-------------------------------')
}
lgb_tipo_modelo_LGB = "LGB"    # LightGBM
lgb_tipo_modelo_LGRF = "LGRF"  # Random Forest con LightGBM
lgb_modelos <- vector(mode = "character", length = NUM_MODELOS)
lgb_filenames <- vector(mode = "character")
if(G_b_crearSubmit) # Solo los modelos de los cmdargs
  for(mi_tipo in paste0(c(lgb_tipo_modelo_LGRF, lgb_tipo_modelo_LGB), lgb_mascara_filename))
    lgb_filenames <- c(lgb_filenames, unique(str_replace(dir(path = s_modelos_path, pattern = paste0(mi_tipo, '.*.modelo')), "_[0-9]*\\.[0-9][0-9][0-9]\\.modelo$", "")))
for(mi_tipo in lgb_filenames)
{
  if(length(dir(path = s_output_path, pattern = paste0(mi_tipo, '_submit.*.csv'))) != 0)
  { lgb_filenames <- lgb_filenames[lgb_filenames != mi_tipo]; next }
  if(length(dir(path = s_modelos_path, pattern = paste0(mi_tipo, '_submit.*.csv'))) != 0)
  { lgb_filenames <- lgb_filenames[lgb_filenames != mi_tipo]; next }
  if(length(dir(path = s_output_path, pattern = paste0(mi_tipo, '_submit.*.zip'))) != 0)
  { lgb_filenames <- lgb_filenames[lgb_filenames != mi_tipo]; next }
  if(length(dir(path = s_modelos_path, pattern = paste0(mi_tipo, '_submit.*.zip'))) != 0)
  { lgb_filenames <- lgb_filenames[lgb_filenames != mi_tipo]; next }
  if(length(dir(path = s_output_path, pattern = paste0(mi_tipo, '_submit.*.rar'))) != 0)
  { lgb_filenames <- lgb_filenames[lgb_filenames != mi_tipo]; next }
  if(length(dir(path = s_modelos_path, pattern = paste0(mi_tipo, '_submit.*.rar'))) != 0)
  { lgb_filenames <- lgb_filenames[lgb_filenames != mi_tipo]; next }
}
for(lgb_filename in lgb_filenames)
{
  print('Encontrado un modelo sin submit...')
  print(lgb_filename)
  tmp_lgb_modelos <- str_replace(dir(path = s_modelos_path, pattern = paste0(str_replace(lgb_filename, "_[0-9]*\\.[0-9][0-9][0-9]$", ""), '.*.modelo')), pattern = "\\.modelo", replacement = "")
  n_version <- as.integer(substring(str_extract(lgb_filename, pattern = "v[0-9]+"), first = 2))
  if(n_version > 500)
    next # Es de los de numAds
  if(length(tmp_lgb_modelos) == NUM_MODELOS)
  {
    # Ordenamos modelos por numModelo (aquí NO da igual, porque ya NO los vamos a promediar)
    lgb_modelos[as.integer(substr(tmp_lgb_modelos, nchar(tmp_lgb_modelos)-2, nchar(tmp_lgb_modelos)))] <- tmp_lgb_modelos
    break # Ok. Pasamos directamente a predecir con estos modelos.
  } else {
    print('NOTA: Lo descartamos porque todavía no están todos los clusters entrenados')
  }
}
n_versiones <- as.integer(substring(str_extract(lgb_filenames, pattern = "v[0-9]+"), first = 2))
for(n_version in unique(n_versiones)[unique(n_versiones) > 500])
{
  nSubmits <- length(dir(path = s_modelos_path, pattern = paste0('.*', lgb_tipo_modelo_LGRF,'.*_v', n_version, '.*.submit.*')))
  nSubmits <- nSubmits + length(dir(path = s_modelos_path, pattern = paste0('.*', lgb_tipo_modelo_LGB,'.*_v', n_version, '.*.submit.*')))
  if(nSubmits != 0)
  {
    print(paste0('NOTA: Descartamos v', n_version, '. Encontrado(s) ', nSubmits, ' submit(s) con esta versión.'))
    break
  } else { print(paste0('Buscando modelos con versión ', n_version)) }
  tmp_lgb_modelos <- vector(mode = "character", length = NUM_MODELOS)
  lgb_filenames_version <- lgb_filenames[n_versiones == n_version]
  for(lgb_filename in lgb_filenames)
  {
    if(n_version != as.integer(substring(str_extract(lgb_filename, pattern = "v[0-9]+"), first = 2)))
      next
    mis_modelos <- dir(path = s_modelos_path, pattern = paste0(lgb_filename, '.*.modelo'))
    for(mi_modelo in mis_modelos)
    {
      numAds <- 1 + as.integer(substring(str_extract(mi_modelo, pattern = "\\.[0-9][0-9][0-9]\\.modelo"), first = 2, last = 4))
      if(tmp_lgb_modelos[numAds - 1] != ""){
        print(paste0('Nota: Encontrado más de un modelo ', str_pad(numAds-1,3,"left","0"), ' v', n_version, ' (numAds = ', numAds, '). Nos quedamos con el menor...'))
        if(tmp_lgb_modelos[numAds - 1] > str_replace(mi_modelo, pattern = "\\.modelo", replacement = ""))
          tmp_lgb_modelos[numAds - 1] <- str_replace(mi_modelo, pattern = "\\.modelo", replacement = "")
      } else {
        tmp_lgb_modelos[numAds - 1] <- str_replace(mi_modelo, pattern = "\\.modelo", replacement = "")
      }
    }
  }
  if(length(tmp_lgb_modelos[tmp_lgb_modelos != ""]) == NUM_MODELOS)
  {
    # Ordenamos modelos por numModelo (aquí NO da igual, porque ya NO los vamos a promediar)
    lgb_modelos[as.integer(substr(tmp_lgb_modelos, nchar(tmp_lgb_modelos)-2, nchar(tmp_lgb_modelos)))] <- tmp_lgb_modelos
    print(paste0('Hay ', NUM_MODELOS, ' modelos v', n_version, ' (sin submit). Pasamos directamente a predecir con ellos [', lgb_filename, '].'))
    print(lgb_modelos)
    break # Ok. Pasamos directamente a predecir con estos modelos.
  } else {
    print(paste0('NOTA: Descartamos v', n_version, ' porque no están todos entrenados (faltan algunos numAds)'))
    print(substring(str_extract(tmp_lgb_modelos, pattern = "\\.[0-9][0-9][0-9]$"), first = 2, last = 4))
  }
}
s_Fic_log <- NULL # Por si acaso no hemos entrado a entrenar...
if(any(lgb_modelos == "") | anyNA(lgb_modelos))
{
  print('-------------------------------')
  print('Entrenamos:')
  print('-------------------------------')

  for(forNumModelo in 1:NUM_MODELOS) # Un modelo por cada "Cluster"...
  { #PRUEBAS forNumModelo <- 1 #PRUEBAS

    if(!G_b_REV) numModelo <- forNumModelo else numModelo <- (NUM_MODELOS - forNumModelo + 1) # Empezamos por el final (numModelo <- NUM_MODELOS - forNumModelo + 1)

    if(!is.na(lgb_modelos[numModelo]) & lgb_modelos[numModelo] != "")
    {
      print(paste0('Warning: Ya hay un modelo ', str_pad(numModelo, 3, "left", "0"), ' (', lgb_modelos[numModelo], '). Pasamos al siguiente...'))
      next # el "continue" de C
    }
    miDescr <- paste0("LGB Training - [cluster ", numModelo, "]")
    numCluster <- numModelo # numModelo + 1 # Era numAdsCluster
    fich_name <- paste0("train_valid_", str_pad(numCluster, 2, "left", "0"), "__", lgb_reduc_seed, "__", lgb_train_porc, ".RData") # "train_valid_nn__ssss__pppp.RData"
    
    if(file.exists(file.path(s_input_path, fich_name)))
    {
      print(paste0(miDescr, ' Batch trainset+validset ALL (cluster == ', numCluster, ')'))
      load(file = file.path(s_input_path, fich_name))
      print(paste0(miDescr, ' Batch trainset+validset ALL (cluster == ', numCluster, ')', ' Ok. ', jjfmt(nrow(trainset)), ' + ', jjfmt(nrow(validset)), ' regs.'))
    } else {
      print(paste0(miDescr, 'Error: fichero [', s_input_path, fich_name, '] no existe (Cf. XGBoost)'))
      stop(paste0(miDescr, 'Fichero [', s_input_path, fich_name, '] no existe (Cf. XGBoost)'))
    }
    if(nrow(trainset) + nrow(validset) > maxTamFullTrainset)
    {
      print(paste0('Hacemos un sample (', jjfmt(maxTamFullTrainset), ') para que vaya más rápido...'))
      mi_n_porc = maxTamFullTrainset / (nrow(trainset) + nrow(validset))
      trainset <- reducir_trainset(mi_set = trainset, n_seed = lgb_reduc_seed, n_porc = mi_n_porc)[[1]]
      validset <- reducir_trainset(mi_set = validset, n_seed = lgb_reduc_seed, n_porc = mi_n_porc)[[1]]
    }

    jjprint(paste0('Ok. trainset: ', jjfmt(nrow(trainset)), ' registros.  validset: ', jjfmt(nrow(validset)), ' registros.'))
    
    # sapply(trainset, uniqueN)[sapply(trainset, uniqueN)==1]
    
    lgb_preps <- lgb_prep_datos(mi_dt = NULL, b_verbose = 1) # Obtenemos solamente la versión (de 0 a 99!!!)
    
    # ========================  LGB (LightGBM) params - START  ==========================
    
    lgb_mi_version <- 000 + lgb_preps[[1]] # Esto viene de lgb_prep_datos (indica cambios en las "features" seleccionadas)
    lgb_mi_version <- 100 + lgb_preps[[1]] # V.100 Sin usar los pesos (lgb_weights_clases) en el training, en lgb.cv() y lightgbm())
    if(maxImportanceNumVars != 0) lgb_mi_version <- 400 + lgb_preps[[1]] # V.400 Con selección de variables (maxImportanceNumVars=100 mejores)
    if(cmdargs_lgb_min_data_in_leaf != -1)
       lgb_mi_version <- lgb_mi_version + 20 # v.x2x # Añadido nuevo param. (lgb_min_data_in_leaf)
    # if(b_conPond == 1) { lgb_mi_version <- lgb_mi_version + 100 } else if(b_conPond == 2) { lgb_mi_version <- lgb_mi_version + 200 } # V.200 / V.300
    
    # lgb_num_class = 9 # lgb_objective = "multiclass"
    lgb_objective = "multiclass" # multi-class classification application, should set num_class as well
          # NOTA para multi-class: Class is represented by a number and should be from 0 to num_class - 1.
    # objective: (application: alias objective, app)
    #   default='regression'
    #   'regression': regression application
    #       'regression_l2': L2 loss, alias=mean_squared_error,mse
    #       'regression_l1': L1 loss, alias=mean_absolute_error,mae
    #       'huber': Huber loss
    #       'fair': Fair loss
    #       'poisson': Poisson regression
    #   'binary': binary classification application
    #   'lambdarank': lambdarank application
    #       The label should be int type in lambdarank tasks, and larger number represent the higher relevance (e.g. 0:bad, 1:fair, 2:good, 3:perfect).
    #       label_gain can be used to set the gain(weight) of int label.
    #   'multiclass': multi-class classification application, should set num_class as well

    lgb_eval_metric = 'multi_logloss' # log loss for multi-class classification
    # metric:
    #   default={'l2' for regression}, {'binary_logloss' for binary classification},{'ndcg' for lambdarank}
    #   'l1': absolute loss, alias='mean_absolute_error', 'mae'
    #   'l2': square loss, alias='mean_squared_error', 'mse'
    #   'l2_root': root square loss, alias='root_mean_squared_error', 'rmse'
    #   'huber': Huber loss
    #   'fair': Fair loss
    #   'poisson': Poisson regression
    #   'ndcg': NDCG
    #   'map': MAP
    #   'auc': AUC
    #   'binary_logloss': log loss
    #   'binary_error': For one sample 0 for correct classification, 1 for error classification.
    #   'multi_logloss': log loss for multi-class classification
    #   'multi_error': error rate for multi-class classification
    
    lgb_learning_rate = 0.001 # 0.3 [0.01 - 0.2] step size of each boosting step (si baja, afina mejor pero habrá que subir nrounds)
      # Learning rate should be tuned according to your training speed and performance tradeoff.
      # It is not a good practice to consider the learning rate as a hyperparameter to tune.
    
    lgb_max_depth = cmdargs_lgb_max_depth # -1 # maximum depth of trees (demasiado alto => overfitting) [max_depth = -1 for infinite depth]
    lgb_num_leaves = cmdargs_lgb_num_leaves # 31 [1 - Inf] Typical: 255, usually {15, 31, 63, 127, 255, 511, 1023, 2047, 4095}. Tips: adjust depth accordingly by allowing a slightly higher depth than the theoretical number of leaves.
      # NOTA: lgb_num_leaves = 2^lgb_max_depth - 1 para parecerse a xgboost.
    
    lgb_nrounds = 20000 # max number of iterations
    if(G_b_DEBUG & lgb_nrounds > 1)  lgb_nrounds = 20
    lgb_sub_row = 1 # alias de bagging_fraction [subsample en xgboost] (muestreo aleatorio por filas si es < 1, para intentar evitar overfitting)
    lgb_feature_fraction = cmdargs_lgb_feature_fraction # [colsample_bytree en xgboost] (muestreo aleatorio por columnas si es < 1, para intentar evitar overfitting)
    
    lgb_min_sum_hessian_in_leaf = cmdargs_lgb_min_sum_hessian_in_leaf # [min_child_weight en xgboost] Typical = 1. Default = 10 [0 - Inf] # Smaller/Larger is not always better. Leave it alone unless you know what you are doing.
    
    lgb_min_gain_to_split = cmdargs_lgb_min_gain_to_split # [gamma en xgboost] 0 [0 - 2] the larger, the more conservative (less over-fitting) the algorithm will be.
    lgb_min_data_in_leaf = cmdargs_lgb_min_data_in_leaf # 20 [1 - 10k] the larger, the more conservative (less over-fitting) the algorithm will be.[desde V.x20]
    lgb_early_stopping_rounds = 300 # stop if performance keeps getting worse consecutively for k rounds.
    lgb_s_tipo_modelo = lgb_tipo_modelo_LGB
    # ===========================  LGB params - END  ==============================
    if(G_b_DEBUG)  lgb_mi_version <- lgb_mi_version + 990000
    # Nombre del fichero para guardar estadísticas (de estos modelos) en un único fichero:
    # lgb_filename <- str_replace(lgb_modelos[numModelo], "_[0-9]*\\.[0-9][0-9][0-9]$", "")
    if(cmdargs_lgb_min_data_in_leaf == -1) { # v.114
       lgb_filename <- paste0(paste(lgb_s_tipo_modelo, G_CV_nFolds,
                                    str_replace(lgb_objective, "^(...).*$", "\\1"),
                                    lgb_learning_rate, (if(lgb_max_depth==-1) 0 else lgb_max_depth), lgb_num_leaves, lgb_sub_row, lgb_feature_fraction, str_replace(lgb_eval_metric, "multi_", "m"), lgb_early_stopping_rounds, lgb_scale_pos_weight, lgb_min_sum_hessian_in_leaf, lgb_min_gain_to_split,
                                    jjfmt(100*lgb_train_porc, 1, '.', ','), # jjfmt(100*nrow(trainset)/(nrow(trainset)+nrow(validset)), 1, '.', ','),
                                    paste0('v', str_pad(lgb_mi_version, 3, "left", "0")),
                                    sep = "_"))
    } else { # v.134
       lgb_filename <- paste0(paste(lgb_s_tipo_modelo, G_CV_nFolds,
                                    str_replace(lgb_objective, "^(...).*$", "\\1"),
                                    lgb_learning_rate, (if(lgb_max_depth==-1) 0 else lgb_max_depth), lgb_num_leaves, lgb_sub_row, lgb_feature_fraction, str_replace(lgb_eval_metric, "multi_", "m"), lgb_early_stopping_rounds, lgb_scale_pos_weight, lgb_min_sum_hessian_in_leaf, lgb_min_gain_to_split, lgb_min_data_in_leaf,
                                    jjfmt(100*lgb_train_porc, 1, '.', ','), # jjfmt(100*nrow(trainset)/(nrow(trainset)+nrow(validset)), 1, '.', ','),
                                    paste0('v', str_pad(lgb_mi_version, 3, "left", "0")),
                                    sep = "_"))
    }
    # Iniciamos log (lgb_filename".log"):
    if(G_b_NO_LOG == FALSE)
    {
      s_Fic_log <- paste0(s_output_path, lgb_filename, '.log')
      if(file.exists(s_Fic_log))
        jjprint('-------------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE) # Línea separadora para el inicio de este log...
    }
    jjprint(paste0('trainset: ', jjfmt(nrow(trainset)), ' registros.  validset: ', jjfmt(nrow(validset)), ' registros. [', Sys.info()["nodename"], ']'), logfile = s_Fic_log, b_add_datetime = TRUE)
    jjprint(paste0('Preparando datos ', str_pad(numModelo, 3, "left", "0"), '...[', lgb_filename, ']'), logfile = s_Fic_log, b_add_datetime = TRUE)
    # Verificamos que no exista el modelo ya entrenado (fichero):
    if(length(dir(path = s_modelos_path, pattern = paste0(lgb_filename, '.*', '.', str_pad(numModelo, 3, "left", "0"), '.modelo'))) != 0)
    {
      lgb_modelos[numModelo] <- dir(path = s_modelos_path, pattern = paste0(lgb_filename, '.*', '.', str_pad(numModelo, 3, "left", "0"), '.modelo'))[1]
      lgb_modelos[numModelo] <- str_replace(lgb_modelos[numModelo], pattern = "\\.modelo", replacement = "")
      jjprint(paste0('Nota: Encontrado modelo ya entrenado. Lo cargamos ', str_pad(numModelo, 3, "left", "0"), ' (', lgb_modelos[numModelo], ') y pasamos al siguiente...'), logfile = s_Fic_log)
      next # el "continue" de C para el bucle for(forNumModelo...
    }
    
    lgb_preps <- lgb_prep_datos(mi_dt = trainset, b_verbose = 1, maxImportanceNumVars = maxImportanceNumVars, s_modelos_path = s_modelos_path
                                , G_maxNumFeatures_OneHotEnc = G_maxNumFeatures_OneHotEnc
                                , obj_var_name = 'Class', s_Fic_log = s_Fic_log)
    # Seleccionamos variables (si hay otros modelos previos del mismo numAd y misma versión y con importance_matrix):
    # importance_vars <- lgb_prep_datos_busca_imp_matrix(lgb_filename) # vector de variables ordenadas por importancia (las primeras son las más importantes)
    # lgb_preps[[2]] <- lgb_prep_datos_selec_vars(mi_dt = trainset, numVars = 50, importance_vars, b_verbose = 1)
    jjprint(colnames(lgb_preps[[2]]), logfile = s_Fic_log, b_add_datetime = FALSE) # Finalmente, mostramos las columnas elegidas
    
    # sapply(lgb_preps[[2]], uniqueN)[sapply(lgb_preps[[2]], uniqueN)==1]
    X <- data.matrix(lgb_preps[[2]][,-1, with=FALSE])
    y <- lgb_preps[[2]]$Class - 1 # NOTA para multi-class: Class is represented by a number and should be from 0 to num_class - 1.
    # midata <- lgb.Dataset(data = X, label = y)
    # Ponderación: Pesos (de cada fila) en función de su clase:
    mis_weights <- NULL
    if(b_conPond == 1) {
      mis_weights <- as.numeric(lgb_weights_clases)[(y+1)]    # Frecs de cada clase (en trainset)
    } else if(b_conPond == 2) {
      mis_weights <- as.numeric(lgb_weights_clases200)[(y+1)] # 1 / Frecs de cada clase (en trainset)
    }
    # Reducimos memoria:
    if(!G_b_DEBUG)  rm(lgb_preps); gc()
    s_misDims <- paste0(jjfmt(dim(X)[1]), ' regs/',
                        jjfmt(dim(X)[2]), ' cols')
    # Entrenamos:
    jjprint(paste0('Entrenando (CV) ', str_pad(numModelo, 3, "left", "0"), '...[', s_misDims, '] [', lgb_filename, ']'), logfile = s_Fic_log)

    mi_tiempo <- system.time({
      # lightgbm::lgb.unloader(wipe = TRUE)
      if(cmdargs_lgb_min_data_in_leaf == -1) { # v.114
         lgb_cv <- lgb.cv(params = list(objective = lgb_objective, metric = lgb_eval_metric, num_class = lgb_num_class),
                          data = X, label = y,
                          nfold = G_CV_nFolds, # prediction = lgb_get_predictions, # SOLAMENTE PARA CV
                          stratified = TRUE,
                          nrounds = lgb_nrounds,
                          scale_pos_weight = lgb_scale_pos_weight,
                          weight = mis_weights, # obj = NULL, eval = NULL,
                          learning_rate = lgb_learning_rate, max_depth = lgb_max_depth, num_leaves = lgb_num_leaves, sub_row = lgb_sub_row, feature_fraction = lgb_feature_fraction,
                          min_sum_hessian_in_leaf = lgb_min_sum_hessian_in_leaf, min_gain_to_split = lgb_min_gain_to_split,
                          record = TRUE, # Boolean, TRUE will record iteration message to booster$record_evals
                          verbose = 1, eval_freq = 1L, showsd = TRUE,
                          # categorical_feature = c('col1', 'col2'),
                          early_stopping_rounds = lgb_early_stopping_rounds,
                          num_threads = G_maxNumCores # detectCores()+0.5)%/%2 debería ser el número real de "cores"
                          )
      } else { # v.134
         lgb_cv <- lgb.cv(params = list(objective = lgb_objective, metric = lgb_eval_metric, num_class = lgb_num_class),
                          data = X, label = y,
                          nfold = G_CV_nFolds, # prediction = lgb_get_predictions, # SOLAMENTE PARA CV
                          stratified = TRUE,
                          nrounds = lgb_nrounds,
                          scale_pos_weight = lgb_scale_pos_weight,
                          weight = mis_weights, # obj = NULL, eval = NULL,
                          learning_rate = lgb_learning_rate, max_depth = lgb_max_depth, num_leaves = lgb_num_leaves, sub_row = lgb_sub_row, feature_fraction = lgb_feature_fraction,
                          min_sum_hessian_in_leaf = lgb_min_sum_hessian_in_leaf, min_gain_to_split = lgb_min_gain_to_split, min_data_in_leaf = lgb_min_data_in_leaf,
                          record = TRUE, # Boolean, TRUE will record iteration message to booster$record_evals
                          verbose = 1, eval_freq = 1L, showsd = TRUE,
                          # categorical_feature = c('col1', 'col2'),
                          early_stopping_rounds = lgb_early_stopping_rounds,
                          num_threads = G_maxNumCores # detectCores()+0.5)%/%2 debería ser el número real de "cores"
                          )
      }                          
    })
    tiempo_cv <- paste0('Tiempo lgb_cv(): ', mi_tiempo['elapsed']/60, ' minutos')
    jjprint(paste0(tiempo_cv, ' [', s_misDims, ']'), logfile = s_Fic_log)
    # if(lgb_get_predictions)
    # {
    #   lgb_preds <- lgb_cv$pred # Predicciones en los N-Folds no usados para entrenar (vector de tamaño nrow(trainset))
    #   lgb_cv <- lgb_cv$dt      # Este es el mismo data.table que devolvía xgb.cv con prediction=FALSE
    # }
    if(lgb_eval_metric == "multi_logloss") {
      # lgb_cv$best_iter y lgb_cv$best_score devuelven -1 si no ha habido early_stop!!!
      lgb_nrounds <- which.min(unlist(lgb_cv$record_evals$valid$multi_logloss[1])) # == lgb_cv$best_iter
      # mi_lgb_cv_train_best_score <- 0
      mi_lgb_cv_test_best_score <- min(unlist(lgb_cv$record_evals$valid$multi_logloss[1])) # == (-1) * lgb_cv$best_score
      jjprint(unlist(lgb_cv$record_evals$valid$multi_logloss[1])[1:min((lgb_nrounds+2),length(unlist(lgb_cv$record_evals$valid$multi_logloss[1])))], logfile = s_Fic_log, b_add_datetime = FALSE) # +2 para ver los 2 sgtes, si los hay...
    } else if(lgb_eval_metric == "multi_error") {
      # lgb_cv$best_iter y lgb_cv$best_score devuelven -1 si no ha habido early_stop!!!
      lgb_nrounds <- which.min(unlist(lgb_cv$record_evals$valid$multi_error[1])) # == lgb_cv$best_iter
      # mi_lgb_cv_train_best_score <- 0
      mi_lgb_cv_test_best_score <- min(unlist(lgb_cv$record_evals$valid$multi_error[1])) # == (-1) * lgb_cv$best_score
      jjprint(unlist(lgb_cv$record_evals$valid$multi_error[1])[1:min((lgb_nrounds+2),length(unlist(lgb_cv$record_evals$valid$multi_error[1])))], logfile = s_Fic_log, b_add_datetime = FALSE) # +2 para ver los 2 sgtes, si los hay...
    # } else if(lgb_eval_metric == "map")
    # {
    #   lgb_nrounds <- which.max(lgb_cv$evaluation_log$test_map_mean)
    #   # mi_lgb_cv_train_best_score <- 0
    #   mi_lgb_cv_test_best_score <- max(lgb_cv$evaluation_log$test_map_mean)
    # } else if(lgb_eval_metric == "error") {
    #   lgb_nrounds <- which.min(lgb_cv$evaluation_log$test_error_mean)
    #   # mi_lgb_cv_train_best_score <- 0
    #   mi_lgb_cv_test_best_score <- min(lgb_cv$evaluation_log$test_error_mean)
    # } else if(lgb_eval_metric == "auc") {
    #   lgb_nrounds <- which.max(lgb_cv$evaluation_log$test_auc_mean)
    #   # mi_lgb_cv_train_best_score <- 0
    #   mi_lgb_cv_test_best_score <- max(lgb_cv$evaluation_log$test_auc_mean)
    } else {                                                      
      try(
        save(lgb_cv, file = paste0(s_output_path, "lgb_cv_temp.RData"))
        , silent = TRUE)
      jjprint(paste0('Best Test eval nrounds = ???? (cf. lgb_cv_temp.RData)[', lgb_filename, ']'), logfile = s_Fic_log)
      stop('eval_metric desconocida (Cf. lgb_cv_temp.RData!')
    }
    jjprint(paste0('Best Test eval nrounds = ', lgb_nrounds, ' [', lgb_filename, ']'), logfile = s_Fic_log)
    # Guardamos lgb_cv para hacer estadísticas de la cross-validation (aparte):
    try(
      save(lgb_cv, file = paste0(s_input_path, lgb_filename, '-lgb_cv.RData'))
      , silent = TRUE)
    lgb_cv_record_evals <- lgb_cv$record_evals
    rm(lgb_cv); gc()
    try(
      save(lgb_cv_record_evals, file = paste0(s_input_path, lgb_filename, '-lgb_cv_recevals.RData'))
      , silent = TRUE)
    rm(lgb_cv_record_evals); gc()
    mi_tiempo <- system.time({
      # lightgbm::lgb.unloader(wipe = TRUE)
      if(cmdargs_lgb_min_data_in_leaf == -1) { # v.114
         lgb <- lightgbm(params = list(objective = lgb_objective, metric = lgb_eval_metric, num_class = lgb_num_class),
                         data = X, label = y,
                         nrounds = lgb_nrounds,
                         scale_pos_weight = lgb_scale_pos_weight,
                         weight = mis_weights, # obj = NULL, eval = NULL,
                         learning_rate = lgb_learning_rate, max_depth = lgb_max_depth, num_leaves = lgb_num_leaves, sub_row = lgb_sub_row, feature_fraction = lgb_feature_fraction,
                         min_sum_hessian_in_leaf = lgb_min_sum_hessian_in_leaf, min_gain_to_split = lgb_min_gain_to_split,
                         record = TRUE, # Boolean, TRUE will record iteration message to booster$record_evals
                         verbose = 1, eval_freq = lgb_nrounds%/%10,
                         # categorical_feature = c('col1', 'col2'),
                         early_stopping_rounds = lgb_early_stopping_rounds,
                         num_threads = G_maxNumCores # detectCores()+0.5)%/%2 debería ser el número real de "cores"
                         , save_name = '' # = paste0(s_modelos_path, lgb_modelos[numModelo] ,'.modelo')
                        )
       } else { # v.134
         lgb <- lightgbm(params = list(objective = lgb_objective, metric = lgb_eval_metric, num_class = lgb_num_class),
                         data = X, label = y,
                         nrounds = lgb_nrounds,
                         scale_pos_weight = lgb_scale_pos_weight,
                         weight = mis_weights, # obj = NULL, eval = NULL,
                         learning_rate = lgb_learning_rate, max_depth = lgb_max_depth, num_leaves = lgb_num_leaves, sub_row = lgb_sub_row, feature_fraction = lgb_feature_fraction,
                         min_sum_hessian_in_leaf = lgb_min_sum_hessian_in_leaf, min_gain_to_split = lgb_min_gain_to_split, min_data_in_leaf = lgb_min_data_in_leaf,
                         record = TRUE, # Boolean, TRUE will record iteration message to booster$record_evals
                         verbose = 1, eval_freq = lgb_nrounds%/%10,
                         # categorical_feature = c('col1', 'col2'),
                         early_stopping_rounds = lgb_early_stopping_rounds,
                         num_threads = G_maxNumCores # detectCores()+0.5)%/%2 debería ser el número real de "cores"
                         , save_name = '' # = paste0(s_modelos_path, lgb_modelos[numModelo] ,'.modelo')
                        )
       }
    })
    tiempo_train <- paste0('Tiempo lightgbm(): ', mi_tiempo['elapsed']/60, ' minutos')
    jjprint(paste0(tiempo_train, ' [', s_misDims, ']'), logfile = s_Fic_log)
    
    # Añadimos nrounds obtenido en la CV (y el numModelo) al nombre de fichero de cada modelo:
    lgb_modelos[numModelo] <- paste0(paste(lgb_filename,
                                           lgb_nrounds,
                                        sep = "_")
                                     , '.', str_pad(numModelo, 3, "left", "0")) # lo guardamos para poder hacer ensemble luego
    jjprint(paste0('Predecimos en trainset para medir nuestro score(', lgb_eval_metric, '):'), logfile = s_Fic_log, b_add_datetime = TRUE)
    mi_tiempo <- system.time({
      preds <- predict(lgb, data = X, reshape = (lgb_num_class > 2))
      if(lgb_num_class == 2){
        trainset[, prob := preds]
      } else {
        for(i in 1:lgb_num_class)
          trainset[, paste0("prob_class", i) := preds[,i]]
      }
    })
    tiempo_predict <- paste0('Tiempo predict(lightgbm, trainset): ', mi_tiempo['elapsed'], ' segundos')
    jjprint(paste0(tiempo_predict, ' [', s_misDims, ']'), logfile = s_Fic_log, b_add_datetime = TRUE)
    # # # Tabla de confusión:
    # # table(trainset[, .(clicked, prob>.5)])
    mi_lgb_train_score <- as.numeric(lgb$record_evals$train$multi_logloss$eval[lgb_nrounds])
    jjprint(paste0('lgb_train_score(', lgb_eval_metric,')(lightgbm) = ', mi_lgb_train_score, '. Cluster = ', numCluster), logfile = s_Fic_log, b_add_datetime = TRUE)
    mi_lgb_train_score <- mi_mlogloss(tr_valid = trainset, b_verbose = 0)
    jjprint(paste0('lgb_train_score(', lgb_eval_metric,')(mi_', lgb_eval_metric,') = ', mi_lgb_train_score, '. Cluster = ', numCluster), logfile = s_Fic_log, b_add_datetime = TRUE)
    # mi_map12 <- Map_12(tr_valid = trainset[, .(display_id, ad_id, clicked, prob)], b_verbose = 0)
    # jjprint(paste0('map12 = ', mi_map12), logfile = s_Fic_log, b_add_datetime = TRUE)
    
    jjprint(paste0('Guardando modelo entrenado (', lgb_modelos[numModelo] ,'.modelo', ')'), logfile = s_Fic_log, b_add_datetime = TRUE)
    if(!G_b_DEBUG)
      lightgbm::lgb.save(booster = lgb, filename = paste0(s_modelos_path, lgb_modelos[numModelo] ,'.modelo'))
    # # Leemos modelo:
    # lgb <- lgb.load(paste0(s_modelos_path, str_replace(lgb_modelos[numModelo], pattern = "\\.modelo", replacement = ""), '.modelo'))
    write.table(x = t(as.data.frame(list(modelo = paste(numModelo, NUM_MODELOS, sep = "/"),
                        dim_x = t(dim(X)),
                        tiempo_cv = tiempo_cv,
                        tiempo_train = tiempo_train,
                        tiempo_predict = tiempo_predict,
                        cols = t(dimnames(X)[[2]]),
                        list(lgb_learning_rate=lgb_learning_rate, lgb_max_depth=lgb_max_depth, lgb_num_leaves=lgb_num_leaves, lgb_sub_row=lgb_sub_row, lgb_feature_fraction=lgb_feature_fraction,
                             lgb_eval_metric=lgb_eval_metric,lgb_early_stopping_rounds=lgb_early_stopping_rounds,
                             lgb_min_sum_hessian_in_leaf=lgb_min_sum_hessian_in_leaf, lgb_min_gain_to_split=lgb_min_gain_to_split, lgb_min_data_in_leaf=lgb_min_data_in_leaf),
                        lgb_train_pond_clases = ifelse(b_conPond == 1, t(lgb_weights_clases), ifelse(b_conPond == 2, t(lgb_weights_clases200), 'No')),
                        # lgb_train_best_score = mi_lgb_train_score,
                        reduc_seed = lgb_reduc_seed,
                        reduc_train_porcent = lgb_train_porc,
                        lgb_cv_bestRound = lgb_nrounds,
                        # lgb_cv_train_best_score = mi_lgb_cv_train_best_score,
                        lgb_cv_test_best_score = mi_lgb_cv_test_best_score,
                        lgb_train_best_score = mi_lgb_train_score))),
              file = paste0(s_output_path, lgb_modelos[numModelo] ,'.txt'), row.names = T, col.names = F, quote = F, sep = " = ")
    
    suppressWarnings( # "Warning: appending column names to file"
      write.table(x = as.data.frame(list(Modelo = lgb_modelos[numModelo])), file = paste0(s_modelos_path, lgb_filename ,'.txt'), row.names = T, col.names = F, quote = F, sep = " = ", append = TRUE)
    )
    
    # Como la importance_matrix tarda bastante en calcularse, no lo hacemos para todos los numModelo:
    # if(numModelo %% 3 == 1) # 1, 4, 7, 10
    # {
    if(Sys.info()["nodename"] %in% c("JJTZAPATA-W10", "JJTZAPATA2-W10")){
      jjprint('NOTA: No se calcula la feature matrix', logfile = s_Fic_log, b_add_datetime = TRUE)
      suppressWarnings( # "Warning: appending column names to file"
        write.table(x = as.data.frame(list(Nota = 'No se calcula la feature matrix')), file = paste0(s_modelos_path, lgb_filename ,'.txt'), row.names = T, col.names = F, quote = F, sep = " = ", append = TRUE)
      )
    } else {
      jjprint('Compute feature importance matrix...', logfile = s_Fic_log, b_add_datetime = TRUE)
      mi_tiempo <- system.time({
        importance_matrix <- lightgbm::lgb.importance(model = lgb, percentage = TRUE)
        # # feature_names = dimnames(X)[[2]]
        rownames(importance_matrix) <- importance_matrix$Feature
      })
      tiempo_imp_matrix <- paste0('Tiempo lgb.importance(): ', mi_tiempo['elapsed'], ' segundos')
      jjprint(paste0(tiempo_imp_matrix, ' [', s_misDims, ']'), logfile = s_Fic_log, b_add_datetime = TRUE)
      suppressWarnings(
        # Gain:
        write.table(x = as.data.frame(importance_matrix)[,c(1,2)], file = paste0(s_modelos_path, lgb_filename ,'.txt'), row.names = T, col.names = T, quote = F, sep = "\t", append = TRUE)
      )
      suppressWarnings(
        # Cover:
        write.table(x = as.data.frame(importance_matrix)[,c(1,3)], file = paste0(s_modelos_path, lgb_filename ,'.txt'), row.names = T, col.names = T, quote = F, sep = "\t", append = TRUE)
      )
      suppressWarnings(
        # Frequency:
        write.table(x = as.data.frame(importance_matrix)[,c(1,4)], file = paste0(s_modelos_path, lgb_filename ,'.txt'), row.names = T, col.names = T, quote = F, sep = "\t", append = TRUE)
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
      lgb_preps <- lgb_prep_datos(mi_dt = validset, b_verbose = 1, maxImportanceNumVars = maxImportanceNumVars, s_modelos_path = s_modelos_path, G_maxNumFeatures_OneHotEnc = G_maxNumFeatures_OneHotEnc, obj_var_name = 'Class', s_Fic_log = s_Fic_log)
      # Seleccionamos variables (si hay otros modelos previos del mismo numAd y misma versión y con importance_matrix):
      # importance_vars <- lgb_prep_datos_busca_imp_matrix(lgb_filename) # vector de variables ordenadas por importancia (las primeras son las más importantes)
      # lgb_preps[[2]] <- lgb_prep_datos_selec_vars(mi_dt = trainset, numVars = 50, importance_vars, b_verbose = 1)
      # jjprint(colnames(lgb_preps[[2]]), logfile = NULL, b_add_datetime = FALSE)
      # sapply(lgb_preps[[2]], uniqueN)
      X <- data.matrix(lgb_preps[[2]][,-1, with=FALSE])
      y <- lgb_preps[[2]]$Class - 1 # NOTA para multi-class: Class is represented by a number and should be from 0 to num_class - 1.
      # midata <- xgb.DMatrix(data.matrix(lgb_preps[[2]][,-1, with=FALSE]), label = lgb_preps[[2]]$clicked, missing = NA)
      rm(lgb_preps); gc()
      jjprint(paste0('Predecimos en validset para medir nuestro score (',lgb_eval_metric ,'):'), logfile = s_Fic_log, b_add_datetime = TRUE)
      mi_tiempo <- system.time({
        preds <- predict(lgb, data = X, reshape = (lgb_num_class > 2))
        if(lgb_num_class == 2){
          validset[, prob := preds]
        } else {
          for(i in 1:lgb_num_class)
            validset[, paste0("prob_class", i) := preds[,i]]
        }
      })
      tiempo_predict <- paste0('Tiempo predict(lightgbm, validset): ', mi_tiempo['elapsed'], ' segundos')
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
      mi_lgb_valid_score <- mi_mlogloss(tr_valid = validset, b_verbose = 1)
    } else {
      jjprint('-------------------------------', logfile = s_Fic_log, b_add_datetime = TRUE)
      jjprint('NO predecimos en validset, porque NO hay validset(hemos usado todo para CV):', logfile = s_Fic_log, b_add_datetime = TRUE)
      jjprint('-------------------------------', logfile = s_Fic_log, b_add_datetime = TRUE)
      mi_lgb_valid_score <- 0
    }
    jjprint(paste0('lgb_train_score(', lgb_eval_metric,')(mi_', lgb_eval_metric,') = ', mi_lgb_train_score, '. Cluster = ', numCluster), logfile = s_Fic_log, b_add_datetime = TRUE)
    jjprint(paste0('lgb_valid_score(', lgb_eval_metric,')(mi_', lgb_eval_metric,') = ', mi_lgb_valid_score, '. Cluster = ', numCluster), logfile = s_Fic_log, b_add_datetime = TRUE)
    jjprint('----------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
    jjprint(data.frame(score_name = lgb_eval_metric, valid_score = mi_lgb_valid_score, train_score = as.numeric(mi_lgb_train_score), CV_test_score = mi_lgb_cv_test_best_score), logfile = s_Fic_log, b_add_datetime = FALSE) # Volvemos a imprimir mi_lgb_train_score (del Training) al lado.
    # jjprint(cbind(mi_map12_valid[2,], map12Train = mi_lgb_train_score, map12cvTest = mi_lgb_cv_test_best_score), logfile = s_Fic_log, b_add_datetime = FALSE) # Volvemos a imprimir mi_lgb_train_score (del Training) al lado.
    # jjprint(stats2, logfile = s_Fic_log, b_add_datetime = FALSE)
    jjprint('----------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
    # Añadimos a los ficheros de cada training los scores del validset:
    write.table(x = c(lgb_valid_score = mi_lgb_valid_score),
                file = paste0(s_output_path, lgb_modelos[numModelo] ,'.txt'), row.names = T, col.names = F, quote = F, sep = " = ", append = TRUE)
    # # Añadimos a los ficheros de cada training un resumen con la mejora (porcentaje sobre el máximo posible teórico):
    # write.table(x = t(stats2[,-3]),
    #             file = paste0(s_output_path, lgb_modelos[numModelo] ,'.txt'), row.names = T, col.names = F, quote = F, sep = " = ", append = TRUE)
    
    minutos <- as.double((proc.time() - systime_ini)['elapsed'])/60
    if(minutos < 60) strtmp <- paste0(minutos, ' minutos en total.', ' - ', 'trainset: ', jjfmt(nrow(trainset)), ' registros.  validset: ', jjfmt(nrow(validset)), ' registros.') else  strtmp <- paste0(minutos/60, ' horas en total.', ' - ', 'trainset: ', jjfmt(nrow(trainset)), ' registros.  validset: ', jjfmt(nrow(validset)), ' registros.')
    jjprint(strtmp, logfile = s_Fic_log, b_add_datetime = TRUE)
    minutos_pend <- (as.double((proc.time() - systime_ini)['elapsed'])/60) * (NUM_MODELOS / numModelo  -  1)
    if(minutos_pend < 60) strtmp <- paste0('Faltan aprox. ', minutos_pend, ' minutos.') else strtmp <- paste0('Faltan aprox. ', minutos_pend/60, ' horas.')
    jjprint(strtmp, logfile = s_Fic_log, b_add_datetime = TRUE)
  }
}

# length(lgb_modelos) <- 1 # Si tenemos menos modelos, podemos crear submit con menos...
if(is.null(s_Fic_log))
{
  # Por si acaso no hemos entrado a entrenar...
  if(G_b_NO_LOG == FALSE)
  {
    s_Fic_log <- paste0(s_output_path, lgb_filename, '.log')
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
  if(!exists("lgb_modelos"))  lgb_modelos <- vector(mode = "character", length = NUM_MODELOS)
  predict_testset(nombres_modelos = lgb_modelos
                  , filename = lgb_filename, s_input_path = s_input_path, s_output_path = s_output_path, s_modelos_path = s_modelos_path
                  , i_sDescr = "LGB (LightGBM) Predicting"
                  , FUN_prep_datos = lgb_prep_datos, prep_datos_b_verbose = 1
                  , obj_var_name = 'Class', num_classes = lgb_num_class, G_maxNumFeatures_OneHotEnc = G_maxNumFeatures_OneHotEnc
                  , FUN_X = data.matrix
                  , FUN_Predict = function(modelo, X){ return(predict(modelo, X, reshape = (lgb_num_class > 2))) }
                  , FUN_loadmodelo = lgb.load
                  , nombreModeloLoaded = ""
                  , NUM_MODELOS = length(lgb_modelos)
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

if(G_b_NO_LOG == FALSE)
{
  mi_log_filename <- paste0('RedefCancer_LightGBM_jjtz_', Sys.info()["nodename"], '.log')
  # for(mi_log_filename in dir(pattern = 'RedefCancer_LightGBM_jjtz_.*.log'))
  {
    mi_log <- readLines(mi_log_filename)
    l_ini <- length(mi_log)
    # Quitamos líneas completas:
    mi_log <- mi_log[mi_log != "[LightGBM] [Info] No further splits with positive gain, best gain: -inf"]
    mi_log <- mi_log[grep(pattern = "\\[LightGBM\\] \\[Info\\] Trained a tree with leaves=[0-9]* and max_depth=[0-9]*", x = mi_log,
         invert = TRUE)]
    # Reducimos tamaño de algunas líneas:
    mis_lins <- grep(pattern = "[LightGBM] [Warning] Stopped training because there are no more leaves that meet the split requirements.", x = mi_log, fixed = T)
    if(length(mis_lins) != 0)
      mi_log[mis_lins] <- "[LightGBM Warn] Stopped training. No more leaves."
    mis_lins <- grep(pattern = '[1] "----------------------------------------------------------------------------"', x = mi_log, fixed = T)
    if(length(mis_lins) != 0)
      mi_log[mis_lins] <- "------------------------------------------------------"
    mis_lins <- grep(pattern = '[1] "-------------------------------"', x = mi_log, fixed = T)
    if(length(mis_lins) != 0)
      mi_log[mis_lins] <- "------------------------"
    # Quitamos líneas entre:
    # regexp=".* - Tiempo lgb_cv(): .*"
    # y
    # regexp=".* - Best Test eval nrounds = .*"
    mis_lins_ini <- grep(pattern = ".* - Tiempo lgb_cv\\(\\): .*", x = mi_log)
    mis_lins_fin <- grep(pattern = ".* - Best Test eval nrounds = .*", x = mi_log) - 1
    if(length(mis_lins_fin) == length(mis_lins_ini) && all(mis_lins_fin >= mis_lins_ini)) {
      mis_lins <- c()
      for(i in seq_along(mis_lins_ini))
        if(mis_lins_fin[i] != mis_lins_ini[i])
          mis_lins <- c(mis_lins, seq.int(from = mis_lins_ini[i] + 1, to = mis_lins_fin[i], by = 1))
      if(length(mis_lins) != 0)
        mi_log <- mi_log[-mis_lins]
    } else {
      print('ERROR al reducir el log[ FALSE == length(mis_lins_fin) == length(mis_lins_ini) && all(mis_lins_fin > mis_lins_ini) ]')
      mi_log <- c(mi_log, 'ERROR al reducir el log[ FALSE == length(mis_lins_fin) == length(mis_lins_ini) && all(mis_lins_fin > mis_lins_ini) ]')
    }
    # Quitamos líneas repetidas:
    N <- length(mi_log)
    mi_log <- mi_log[mi_log[2:N] != mi_log[1:(N-1)]]
    # Si hay cambios, guardamos:
    if(l_ini != length(mi_log))
    {
      print(paste0('Ok. Log [', mi_log_filename, '] reducido.'))
      writeLines(mi_log, con = str_replace(mi_log_filename,
                                           pattern = 'RedefCancer_LightGBM_jjtz_',
                                           replacement = 'RedefCancer_LightGBM_jjtz-')) # Backup, por si acaso falla la siguiente línea...
      writeLines(mi_log, con = mi_log_filename)
    }
  }
}