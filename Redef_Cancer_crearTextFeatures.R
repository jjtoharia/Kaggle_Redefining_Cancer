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

# # install.packages("tm")
# library(tm)
# install.packages("quanteda")
library(quanteda)
# install.packages("topicmodels")
library(topicmodels)

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

memory.limit(size = 32000)

systime_ini <- proc.time()

# --------------------------------------------------------
G_b_DEBUG <- FALSE # Reducimos todo para hacer pruebas más rápido
# NUM_BLOQUES <- 3

b_CV_fit <- FALSE # Clasif. en trainset/validset para determinar k vía CV
b_CrearFeats <- TRUE # Clasif. en full_trainset/testset y guardar!

LDA_prcnt <- 1;   LDA_seed <- 17    # 1921*LDA_prcnt
LDA_k_min <- 45;     LDA_k_max <- 45;  LDA_k_step <- 10
b_palabs <- F;       b_frases <- F;    b_caract <- T
mi_iter.max <- 100 # 1000 # Solo para entrenar
train_repeticiones = 1 # Solo para entrenar (repeticiones del LDA para obtener una perplexity promedio)
mi_batch_size <- 1921
mis_ngrams <- list(
  frases_n_grams = 1:3,  frases_skip_grams = 0:2,  # Hay menos frases, así que genero más n-gramas
  palabs_n_grams = 1,    palabs_skip_grams = 0,  # palabs_n_grams <- 2:3,  palabs_skip_grams <- 0:2
  caract_n_grams = 1:6,  caract_skip_grams = 0)  # Idealmente 1:20... (avg.wrd.len=5.4+3.3; avg.wrd.len.sin.stopwords=6.7+3.3)
# --------------------------------------------------------
# s_Fic_log <- NULL # No escribir log
prefx <- 'LDA'
if(!(b_palabs & b_frases & b_caract))
{
  if(b_frases) prefx <- paste0(prefx, 'f') # 'LDAf'
  if(b_palabs) prefx <- paste0(prefx, 'p') # 'LDAp' (o 'LDAfp')
  if(b_caract) prefx <- paste0(prefx, 'c') # 'LDAc' (o 'LDAfc' o 'LDApc')
}
s_Fichname <- paste(c(prefx, LDA_prcnt, LDA_seed, LDA_k_min, LDA_k_max
                            ,paste0(min(mis_ngrams[["frases_n_grams"]]),'-',max(mis_ngrams[["frases_n_grams"]]),'-',min(mis_ngrams[["frases_skip_grams"]]),'-',max(mis_ngrams[["frases_skip_grams"]]))
                            ,paste0(min(mis_ngrams[["palabs_n_grams"]]),'-',max(mis_ngrams[["palabs_n_grams"]]),'-',min(mis_ngrams[["palabs_skip_grams"]]),'-',max(mis_ngrams[["palabs_skip_grams"]]))
                            ,paste0(min(mis_ngrams[["caract_n_grams"]]),'-',max(mis_ngrams[["caract_n_grams"]]),'-',min(mis_ngrams[["caract_skip_grams"]]),'-',max(mis_ngrams[["caract_skip_grams"]]))
                            ,mi_iter.max, mi_batch_size
                      ), collapse = '_')
s_Fic_log <- paste0(s_Fichname, '.log')
if(file.exists(s_Fic_log))
  jjprint('######################################################################################', logfile = s_Fic_log, b_add_datetime = FALSE) # Línea separadora para el inicio de este log...

# ##################################################
# ## CARGAMOS TEXTOS:
# ##################################################
if(b_CrearFeats)
{
  # Preparamos y leemos test:
  stopifnot(mi_load_uniq_txt('test', s_Fic_log = s_Fic_log))
  testset_txt <- dataset_txt
}
# Preparamos y leemos train:
stopifnot(mi_load_uniq_txt('train', s_Fic_log = s_Fic_log))
if(b_CrearFeats)
  full_trainset_txt <- dataset_txt

if(LDA_prcnt < 1)
{
  jjprint(paste0('Usamos muestra del ', 100 * LDA_prcnt, '% de los regs...'), logfile = s_Fic_log)
  dataset_txt <- reducir_trainset(mi_set = dataset_txt, id_fld = pk,
                                  n_seed = LDA_seed, n_porc = LDA_prcnt)[[1]]
}
if(b_CV_fit)
{
  jjprint(paste0('Trainset_txt/validset_txt (50/50)...'), logfile = s_Fic_log)
  validset_txt <- reducir_trainset(mi_set = dataset_txt, id_fld = pk,
                                   n_seed = LDA_seed, n_porc = 0.5)
  rm(dataset_txt); gc()
  trainset_txt <- validset_txt[[1]]
  validset_txt <- validset_txt[[2]]
  jjprint(paste0('Ok. trainset_txt: ', jjfmt(nrow(trainset_txt)), ' registros.  validset_txt: ', jjfmt(nrow(validset_txt)), ' registros.'), logfile = s_Fic_log)
} else {
  trainset_txt <- dataset_txt
  rm(dataset_txt); gc()
  jjprint(paste0('Ok. trainset_txt: ', jjfmt(nrow(trainset_txt)), ' registros.  validset_txt: (no hay).'), logfile = s_Fic_log)
}
# ##################################################
# ## CLASIFICAMOS (NO SUPERVISADA):
# ##################################################
# Estas funciones hay que definirlas aquí porque algunos valores por defecto son variables:
crearTextFeatures_train <- function(dataset_txt, palabs_k, frases_k, caract_k, mis_ngrams, s_Fic_log) {
  crearTextFeatures(dataset_txt = dataset_txt, b_toLower = FALSE, s_Fic_log = s_Fic_log,
                    b_entrenar = TRUE, batch_size = mi_batch_size, mi_iter.max = mi_iter.max, train_repeticiones = train_repeticiones,
                    palabs_num_features = palabs_k, frases_num_features = frases_k, caract_num_features = caract_k,
                    mis_ngrams = mis_ngrams)
}
crearTextFeatures_fit <- function(txtWhat, dataset_txt, mi_ListaLDA, palabs_k, frases_k, caract_k, mis_ngrams, s_Fic_log) {
  crearTextFeatures(dataset_txt = dataset_txt, b_toLower = FALSE, s_Fic_log = s_Fic_log,
                    b_entrenar = FALSE, batch_size = 200, # Para evitar problemas de memoria, por si acaso, para xxx_fit() conviene bajar el batchsize pq no se quitan columnas al crear el corpus...
                    palabs_num_features = palabs_k, frases_num_features = frases_k, caract_num_features = caract_k,
                    mis_ngrams = mis_ngrams, ListaLDAfit = mi_ListaLDA, txtWhat=txtWhat)
}

full_dfStats <- NULL
for(k in seq.int(LDA_k_min,LDA_k_max,LDA_k_step))
{ # k=LDA_k_min
  palabs_k <- ifelse(b_palabs, k, 0) # Solo para entrenar (pero el 0 vale para clasificar tb)
  frases_k <- ifelse(b_frases, k, 0) # Solo para entrenar (pero el 0 vale para clasificar tb)
  caract_k <- ifelse(b_caract, k, 0) # Solo para entrenar (pero el 0 vale para clasificar tb)
  jjprint('---------------------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
  jjprint(paste0('palabs_k=', palabs_k, '. frases_k=', frases_k, '. caract_k=', caract_k,'. ', 'INICIO'))
  mi_tiempo <- system.time({
    mi_dfStats <- NULL
    mifich <- paste0(s_Fichname, '_00_', palabs_k, '-', frases_k, '-', caract_k, '.RData')
    if(!mi_load(s_input_path, mifich))
    {
      tmpRetVal <- 
        crearTextFeatures_train(trainset_txt, palabs_k, frases_k, caract_k, mis_ngrams, s_Fic_log)
        # crearTextFeatures(dataset_txt = trainset_txt, b_toLower = FALSE, s_Fic_log = s_Fic_log,
        #                   b_entrenar = TRUE, batch_size = mi_batch_size, mi_iter.max = mi_iter.max, train_repeticiones = train_repeticiones,
        #                   palabs_num_features = palabs_k, frases_num_features = frases_k, caract_num_features = caract_k,
        #                   mis_ngrams = mis_ngrams)
      mi_ListaLDA <- tmpRetVal[[1]]
      mi_dfStats <- tmpRetVal[[2]]; rm(tmpRetVal)
      save(mi_ListaLDA, mi_dfStats, file = paste0(s_input_path, mifich))
      # b_CrearFeats <- FALSE # Para evitar problemas de memoria, por si acaso, lo hacemos en dos pasos... (para xxx_fit() conviene bajar el batchsize pq no se quitan columnas al crear el corpus...)
    }
    if(!is.null(mi_dfStats)) {
      jjprint(mi_dfStats, logfile = s_Fic_log, b_add_datetime = FALSE)
      if(is.null(full_dfStats)) {
        full_dfStats <- mi_dfStats
      } else {
        full_dfStats <- rbind(full_dfStats, mi_dfStats)
      }
    } else {
      for(tok_name in c("frases", "palabs", "caract"))
      { # tok_name = "palabs"
        if(!is.null(mi_ListaLDA[[tok_name]]))
        {
          tok_descr <- paste0('[',tok_name,'] ')
          jjprint(paste0(tok_descr, "         k = ", mi_ListaLDA[[tok_name]]@k), logfile = s_Fic_log)
          jjprint(paste0(tok_descr, "     alpha = ", jjfmt(mi_ListaLDA[[tok_name]]@alpha, 4)), logfile = s_Fic_log)
          jjprint(paste0(tok_descr, "       Dim = (", paste(jjfmt(mi_ListaLDA[[tok_name]]@Dim), collapse = ' - '), ')'), logfile = s_Fic_log)
          jjprint(paste0(tok_descr, "__Perplex. = ", jjfmt(perplexity(mi_ListaLDA[[tok_name]]), 2)), logfile = s_Fic_log)
        }
      }
    }
  })
  mi_tiempo_train_LDA <- paste0('Tiempo TOTAL Train LDA(): ', jjfmt(mi_tiempo['elapsed']/60, 3), ' minutos')
  jjprint(paste0('palabs_k=', palabs_k, '. frases_k=', frases_k, '. caract_k=', caract_k,'. ', mi_tiempo_train_LDA), logfile = s_Fic_log)
  
  if(b_CV_fit)
  {
    # mi_tiempo <- system.time({
    #   tmpRetVal <- 
    #     crearTextFeatures_fit(txtWhat='Clasif1 LDA(train)', trainset_txt, mi_ListaLDA, palabs_k, frases_k, caract_k, mis_ngrams, s_Fic_log)
    #     # crearTextFeatures(dataset_txt = trainset_txt, b_toLower = FALSE, s_Fic_log = s_Fic_log,
    #     #                   b_entrenar = FALSE, batch_size = 16,
    #     #                   palabs_num_features = palabs_k, frases_num_features = frases_k, caract_num_features = caract_k,
    #     #                   mis_ngrams = mis_ngrams, ListaLDAfit = mi_ListaLDA)
    # })
    # trainset_txt2 <- tmpRetVal[[1]]
    # mi_dfStats <- tmpRetVal[[2]]; rm(tmpRetVal)
    # jjprint(mi_dfStats, logfile = s_Fic_log, b_add_datetime = FALSE)
    # mi_tiempo_clasif1_LDA <- paste0('Tiempo TOTAL Clasif1 LDA(train): ', jjfmt(mi_tiempo['elapsed']/60, 3), ' minutos')
    # jjprint(paste0('palabs_k=', palabs_k, '. frases_k=', frases_k, '. caract_k=', caract_k,'. ', mi_tiempo_clasif1_LDA, '(', nrow(trainset_txt), ' regs)'), logfile = s_Fic_log)
    mi_tiempo <- system.time({
      tmpRetVal <- 
        crearTextFeatures_fit(txtWhat='Clasif2 LDA(valid)', validset_txt, mi_ListaLDA, palabs_k, frases_k, caract_k, mis_ngrams, s_Fic_log)
      # crearTextFeatures(dataset_txt = validset_txt, b_toLower = FALSE, s_Fic_log = s_Fic_log,
      #                   b_entrenar = FALSE, batch_size = 16,
      #                   palabs_num_features = palabs_k, frases_num_features = frases_k, caract_num_features = caract_k,
      #                   mis_ngrams = mis_ngrams, ListaLDAfit = mi_ListaLDA)
    })
    validset_txt2 <- tmpRetVal[[1]]
    mi_dfStats <- tmpRetVal[[2]]; rm(tmpRetVal)
    jjprint(mi_dfStats, logfile = s_Fic_log, b_add_datetime = FALSE)
    if(is.null(full_dfStats)) {
      full_dfStats <- mi_dfStats
    } else {
      full_dfStats <- rbind(full_dfStats, mi_dfStats)
    }
    mi_tiempo_clasif2_LDA <- paste0('Tiempo TOTAL Clasif2 LDA(valid): ', jjfmt(mi_tiempo['elapsed']/60, 3), ' minutos')
    jjprint(paste0('palabs_k=', palabs_k, '. frases_k=', frases_k, '. caract_k=', caract_k,'. ', mi_tiempo_clasif2_LDA, '(', nrow(validset_txt), ' regs)'), logfile = s_Fic_log)
  }
  
  if(b_CrearFeats)
  {
    mi_tiempo <- system.time({
      tmpRetVal <- 
        crearTextFeatures_fit(txtWhat='Clasif3 LDA( FULL_TRAINSET )', full_trainset_txt, mi_ListaLDA, palabs_k, frases_k, caract_k, mis_ngrams, s_Fic_log)
      # crearTextFeatures(dataset_txt = validset_txt, b_toLower = FALSE, s_Fic_log = s_Fic_log,
      #                   b_entrenar = FALSE, batch_size = 16,
      #                   palabs_num_features = palabs_k, frases_num_features = frases_k, caract_num_features = caract_k,
      #                   mis_ngrams = mis_ngrams, ListaLDAfit = mi_ListaLDA)
    })
    full_trainset_txt2 <- tmpRetVal[[1]]
    mi_dfStats <- tmpRetVal[[2]]; rm(tmpRetVal)
    mi_tiempo_clasif3_LDA <- paste0('Tiempo TOTAL Clasif3 LDA( FULL_TRAINSET ): ', jjfmt(mi_tiempo['elapsed']/60, 3), ' minutos')
    mi_tiempo <- system.time({
      tmpRetVal <- 
        crearTextFeatures_fit(txtWhat='Clasif4 LDA( TESTSET )', testset_txt, mi_ListaLDA, palabs_k, frases_k, caract_k, mis_ngrams, s_Fic_log)
      # crearTextFeatures(dataset_txt = validset_txt, b_toLower = FALSE, s_Fic_log = s_Fic_log,
      #                   b_entrenar = FALSE, batch_size = 16,
      #                   palabs_num_features = palabs_k, frases_num_features = frases_k, caract_num_features = caract_k,
      #                   mis_ngrams = mis_ngrams, ListaLDAfit = mi_ListaLDA)
    })
    testset_txt2 <- tmpRetVal[[1]]
    mi_dfStats <- tmpRetVal[[2]]; rm(tmpRetVal)
    mi_tiempo_clasif4_LDA <- paste0('Tiempo TOTAL Clasif4 LDA( TESTSET ): ', jjfmt(mi_tiempo['elapsed']/60, 3), ' minutos')
  }
  jjprint(paste0('palabs_k=', palabs_k, '. frases_k=', frases_k, '. caract_k=', caract_k,'. ', mi_tiempo_train_LDA,' (', nrow(trainset_txt), ' regs)'), logfile = s_Fic_log)
  if(b_CV_fit)  if(exists('trainset_txt2')) jjprint(paste0('palabs_k=', palabs_k, '. frases_k=', frases_k, '. caract_k=', caract_k,'. ', mi_tiempo_clasif1_LDA, '(', nrow(trainset_txt), ' regs)'), logfile = s_Fic_log)
  if(b_CV_fit)  jjprint(paste0('palabs_k=', palabs_k, '. frases_k=', frases_k, '. caract_k=', caract_k,'. ', mi_tiempo_clasif2_LDA, '(', nrow(validset_txt), ' regs)'), logfile = s_Fic_log)
  if(b_CrearFeats)  jjprint(paste0('palabs_k=', palabs_k, '. frases_k=', frases_k, '. caract_k=', caract_k,'. ', mi_tiempo_clasif3_LDA, '(', nrow(full_trainset_txt), ' regs)'), logfile = s_Fic_log)
  if(b_CrearFeats)  jjprint(paste0('palabs_k=', palabs_k, '. frases_k=', frases_k, '. caract_k=', caract_k,'. ', mi_tiempo_clasif4_LDA, '(', nrow(testset_txt), ' regs)'), logfile = s_Fic_log)
  # if(b_CV_fit)  jjprint(paste0('summary(trainset_txt2[,3:',(2+palabs_k+frases_k+caract_k),']):'), logfile = s_Fic_log)
  # if(b_CV_fit)  jjprint(summary(trainset_txt2[,3:(2+palabs_k+frases_k+caract_k)]), logfile = s_Fic_log, b_add_datetime = FALSE)
  # if(b_CV_fit)  jjprint(paste0('summary(validset_txt2[,3:',(2+palabs_k+frases_k+caract_k),']):'), logfile = s_Fic_log)
  # if(b_CV_fit)  jjprint(summary(validset_txt2[,3:(2+palabs_k+frases_k+caract_k)]), logfile = s_Fic_log, b_add_datetime = FALSE)
  # if(b_CrearFeats)  jjprint(paste0('summary(full_trainset_txt2[,3:',(2+palabs_k+frases_k+caract_k),']):'), logfile = s_Fic_log)
  # if(b_CrearFeats)  jjprint(summary(full_trainset_txt2[,3:(2+palabs_k+frases_k+caract_k)]), logfile = s_Fic_log, b_add_datetime = FALSE)
  # if(b_CrearFeats)  jjprint(paste0('summary(testset_txt2[,3:',(2+palabs_k+frases_k+caract_k),']):'), logfile = s_Fic_log)
  # if(b_CrearFeats)  jjprint(summary(testset_txt2[,3:(2+palabs_k+frases_k+caract_k)]), logfile = s_Fic_log, b_add_datetime = FALSE)
  if(b_CV_fit) {
    if(exists('trainset_txt2')) jjprint(paste0('sapply(trainset_txt2, uniqueN):'), logfile = s_Fic_log)
    if(exists('trainset_txt2')) jjprint(sapply(trainset_txt2, uniqueN), logfile = s_Fic_log, b_add_datetime = FALSE)
    jjprint(paste0('sapply(validset_txt2, uniqueN):'), logfile = s_Fic_log)
    jjprint(sapply(validset_txt2, uniqueN), logfile = s_Fic_log, b_add_datetime = FALSE)
  }
  if(b_CrearFeats) {
    full_trainset_txt2[, texto := NULL]
    testset_txt2[, texto := NULL]
    fich_TRN_TST <- paste0(s_Fichname, '_LDA_', palabs_k, '-', frases_k, '-', caract_k, '.RData')
    Fich_TRN_txt <- full_trainset_txt2
    Fich_TST_txt <- testset_txt2
    save(Fich_TRN_txt, Fich_TST_txt, file = paste0(s_input_path, fich_TRN_TST))
    jjprint(paste0("Ok: stopifnot(mi_load(s_input_path, '",fich_TRN_TST,"')) # Cargar Fich_TRN_txt y Fich_TST_txt"), logfile = s_Fic_log)
    rm(full_trainset_txt2, testset_txt2); gc()
    jjprint(paste0('sapply(Fich_TRN_txt, uniqueN):'), logfile = s_Fic_log)
    jjprint(sapply(Fich_TRN_txt, uniqueN), logfile = s_Fic_log, b_add_datetime = FALSE)
    jjprint(paste0('sapply(Fich_TST_txt, uniqueN):'), logfile = s_Fic_log)
    jjprint(sapply(Fich_TST_txt, uniqueN), logfile = s_Fic_log, b_add_datetime = FALSE)
    jjprint(paste0("Ok: stopifnot(mi_load(s_input_path, '",fich_TRN_TST,"')) # Cargar Fich_TRN_txt y Fich_TST_txt"), logfile = s_Fic_log)
  }
  jjprint(paste0('FIN. palabs_k = ', palabs_k, '. frases_k = ', frases_k, '. caract_k = ', caract_k,
                 '. Total features nuevas = ', palabs_k + frases_k + caract_k), logfile = s_Fic_log)
  jjprint('---------------------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
  if(!is.null(full_dfStats))
  {
    jjprint(full_dfStats, logfile = s_Fic_log, b_add_datetime = FALSE)
    jjprint('---------------------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
    stats_file <- 'LDAc_1_17_9_xxx_1-3-0-2_1-1-0-0_1-x-0-0_1000_1000_1921_stats'
    # Añadimos stats al log de stats:
    jjprint(full_dfStats, logfile = paste0(stats_file, '.log'), b_add_datetime = FALSE)
    # Añadimos stats al csv (en Out/):
    write.table(full_dfStats, file = paste0(s_output_path, stats_file, '.csv'), append = TRUE, quote = TRUE, sep = ",",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = !file.exists(paste0(s_output_path, stats_file, '.csv')),
                fileEncoding = "UTF-8")
    # full_dfStats_from_log <- read.table(file = paste0(s_output_path, '../', stats_file, '.log'), sep = '\t', stringsAsFactors = FALSE, fileEncoding = "UTF-8", header = TRUE)
    # full_dfStats_from_csv <- read.table(file = paste0(s_output_path, stats_file, '.csv'), stringsAsFactors = FALSE, sep = ",", fileEncoding = "UTF-8", header = TRUE)
    
    full_dfStats <- NULL # Reseteamos las stats para el próximo "k".
  }
} # fin de for(k in seq.int(LDA_k_min,LDA_k_max,LDA_k_step))

b_VerStats <- FALSE
if(b_VerStats)
{
  stats_file <- 'LDAf_1_17_9_xxx_1-3-0-2_1-1-0-0_1000_1921_stats'; tok_name <- "frases"
  stats_file <- 'LDAp_1_17_9_xxx_1-3-0-2_1-1-0-0_1000_1921_stats'; tok_name <- "palabs"
  stats_file <- 'LDAc_1_17_9_xxx_1-3-0-2_1-1-0-0_1-x-0-0_1000_1000_1921_stats'; tok_name <- "caract"
  tok_descr <- paste0('[',tok_name,'] ')
  # full_dfStats_from_log <- read.table(file = paste0(s_output_path, '../', stats_file, '.log'), sep = '\t', stringsAsFactors = FALSE, fileEncoding = "UTF-8", header = TRUE)
  full_dfStats_from_csv <- read.table(file = paste0(s_output_path, stats_file, '.csv'), stringsAsFactors = FALSE, sep = ",", fileEncoding = "UTF-8", header = TRUE)
  mi_df <- full_dfStats_from_csv[full_dfStats_from_csv$Perplex_valid > 0 & full_dfStats_from_csv$tok_name == tok_name,]
  mi_df <- mi_df[order(mi_df$model_k),]
  plot(x=mi_df$model_k, y=mi_df$Perplex_valid, type = "l", ylab="Perplex_Valid", xlab = "k (topics)")
  x_vals <- mi_df$model_k[2:(nrow(mi_df))]
  y_vals <- mi_df$Perplex_valid[2:(nrow(mi_df))]
  ydif_vals <- mi_df$Perplex_valid[2:(nrow(mi_df))] - mi_df$Perplex_valid[1:(nrow(mi_df)-1)]
  yRpc_vals <- abs(ydif_vals / (mi_df$model_k[2:(nrow(mi_df))] - mi_df$model_k[1:(nrow(mi_df)-1)]))
  ydifdif_vals <- ydif_vals[2:(length(ydif_vals))] - ydif_vals[1:(length(ydif_vals)-1)]
  plot(x=x_vals, y=ydif_vals, type = "l", ylab="Perplex_Valid Dif.", xlab = "k (topics)")
  plot(x=x_vals, y=yRpc_vals, type = "l", ylab="Perplex_Valid RPC", xlab = "k (topics)")
  plot(x=x_vals[2:length(x_vals)], y=ydifdif_vals, type = "l", ylab="Perplex_Valid Dif.-Dif.", xlab = "k (topics)")
  k_opt = x_vals[which.min(y_vals)]
  print(paste0(tok_descr, 'which.min(y_vals): ', k_opt))
  k_opt = x_vals[which.max(ydif_vals)]
  print(paste0(tok_descr, 'which.max(ydif_vals): ', k_opt))
  
  k_opt = x_vals[2:length(x_vals)][which.min(ydifdif_vals)]
  print(paste0(tok_descr, 'which.min(ydif-dif_vals): ', k_opt))
  
  if(nrow(mi_df) > 9) x9_vals <- mi_df$model_k[9:(nrow(mi_df))]
  if(nrow(mi_df) > 18) x18_vals <- mi_df$model_k[18:(nrow(mi_df))]
  if(nrow(mi_df) > 9) ydif9_vals <- mi_df$Perplex_valid[9:(nrow(mi_df))] - mi_df$Perplex_valid[1:(nrow(mi_df)-8)]
  if(nrow(mi_df) > 18) ydif18_vals <- mi_df$Perplex_valid[18:(nrow(mi_df))] - mi_df$Perplex_valid[1:(nrow(mi_df)-17)]
  if(nrow(mi_df) > 18) ydif9dif_vals <- ydif9_vals[9:(length(ydif9_vals))] - ydif9_vals[1:(length(ydif9_vals)-8)]
  if(nrow(mi_df) > 18) plot(x=x9_vals, y=ydif9_vals, type = "l", ylab="Perplex_Valid Dif-9.", xlab = "k (topics)")
  if(nrow(mi_df) > 18) plot(x=x18_vals, y=ydif18_vals, type = "l", ylab="Perplex_Valid Dif-18.", xlab = "k (topics)")
  if(nrow(mi_df) > 18) plot(x=x9_vals[9:length(x9_vals)], y=ydif9dif_vals, type = "l", ylab="Perplex_Valid Dif-9.-Dif.", xlab = "k (topics)")
  
  if(nrow(mi_df) > 9) k_opt = x9_vals[which.max(ydif9_vals)]
  if(nrow(mi_df) > 9) print(paste0(tok_descr, 'which.max(ydif-9_vals): ', k_opt))
  if(nrow(mi_df) > 18) k_opt = x18_vals[which.max(ydif18_vals)]
  if(nrow(mi_df) > 18) print(paste0(tok_descr, 'which.max(ydif-18_vals): ', k_opt))
  
  if(nrow(mi_df) > 18) k_opt = x9_vals[9:length(x9_vals)][which.min(ydif9dif_vals)]
  if(nrow(mi_df) > 18) print(paste0(tok_descr, 'which.min(ydif-9.dif_vals): ', k_opt))
}

minutos_acum <- as.double((proc.time() - systime_ini)['elapsed'])/60
jjprint(paste0(minutos_acum, ' minutos en total.'), logfile = s_Fic_log)
jjprint(paste0('Ok.', logfile = s_Fic_log))

# cleanup:
try(registerDoSEQ(), silent = TRUE) # library(doParallel) [turn parallel processing off and run sequentially again]
try(stopImplicitCluster(), silent = TRUE)
try(stopCluster(cl), silent = TRUE)
# if(G_b_DEBUG)
# {
#   # library(devtools)
#   # setRepositories()
#   # devtools::install_github("nicolaroberts/hdp", build_vignettes = TRUE)
#  
#   # install.packages("DPpackage")
#   library(DPpackage)
# }