# options(echo = FALSE) # ECHO OFF
###############################################################
#           Redefining Cancer Treatment - JJTZ 2017
###############################################################
# ##################################################
# ## Funciones útiles:
# ##################################################
mi_source <- function(x){
  if(file.exists(x)){ source(x, encoding = "UTF-8")
  } else {
    x <- paste0('../', x)
    if(file.exists(x)){ source(x, encoding = "UTF-8")
    } else {
      x <- paste0('../', x)
      if(file.exists(x)){ source(x, encoding = "UTF-8")
} }}}
mi_source("funciones_utiles.R")

jjprint <- function(x, logfile = NULL, b_add_datetime = TRUE)
{
  if(b_add_datetime) x <- paste0(Sys.time(), ' - ', x)
  print(x)
  if(!is.null(logfile)){
    if(is.data.frame(x)){
      suppressWarnings( # "Warning: appending column names to file"
        write.table(x, file = logfile, row.names = FALSE, append = TRUE, sep = "\t")
      )
    } else {
      write(x, file = logfile, append = TRUE)
    }
  }
}
mi_load <- function(s_input_path, x){
  xn <- sub(pattern = '.RData$', replacement = '', x)
  xf <- paste0(xn, '.RData')
  if(file.exists(file.path(s_input_path, xf))) {
    jjprint(paste0('Leyendo ', xn, '...'))
    load(file = file.path(s_input_path, xf), envir = .GlobalEnv)
    return(TRUE)
  } else {
    warning(paste0('Error leyendo ', xf, '.'))
    return(FALSE)
  }
}
# ##################################################
# ## Funciones:
# ##################################################
guardar_submit <- function(testset, fichero = "submit.csv", b_restore_key = FALSE, b_write_file_append = FALSE, b_verbose = 1, s_output_path = "") # b_verbose=0,1,2
{
  if(!exists("G_b_DEBUG"))  G_b_DEBUG <- FALSE
  # # Cambiamos nombres de campos: (prob_classN -> classN)
  # colnames(testset)[colnames(testset) %in% paste0("prob_class", 1:9)] <- paste0("class", 1:9)
  # ----------------------------------------------
  # PREPARACIÓN DEL SUBMIT:
  # ----------------------------------------------
  if(all(paste0("class", 1:9) %in% colnames(testset))) {
    submitset <- testset[,c("ID", paste0("class", 1:9))]
  } else if(all(paste0("prob_class", 1:9) %in% colnames(testset))) {
    submitset <- testset[,c("ID", paste0("prob_class", 1:9))]
    colnames(submitset) <- c("ID", paste0("class", 1:9))
  }
  stopifnot(all(c("ID", paste0("class", 1:9)) %in% colnames(submitset)))
  # if(b_verbose == 2)  print('Calculando submitset...')
  # mi_tiempo <- system.time({
  # submitset <- fread(file.path(s_input_path, "submissionFile"), colClasses = c("integer", rep("numeric", 9)))
  # for(class in 1:9)
  # {
  #   fld <- paste0("class", class)
  # 
  #   prob <- 1 # class / 10 # PENDIENTE
  # 
  #   submitset[, (fld) := prob]
  # }
  # })
  # if(b_verbose == 2)  print(mi_tiempo['elapsed'])
  # head(submitset)
  
  if(G_b_DEBUG)  fichero <- paste0('Debug_', fichero)

  if(b_verbose == 2)  print(paste0('Guardamos fichero ', paste0(s_output_path, fichero), '...'))
  mi_tiempo <- system.time({
    # Ordenamos por ID:
    setkey(submitset, ID)
    suppressWarnings(
      write.table(submitset, file = paste0(s_output_path, fichero), append = b_write_file_append, row.names=F, col.names=!b_write_file_append, quote=F, sep=",")
    )
  })
  if(b_verbose == 2)  print(mi_tiempo['elapsed'])
  if(b_verbose >= 1)  print(paste0('Ok. Fichero (', paste0(s_output_path, fichero), ') guardado.'))
}

basic_preds_guardar_submit <- function(trainset, validset = NULL, k = 1, reg = 0, i_sFileSuffix = "", s_output_path = "") # k=1,2,3
{
  # -----------------------------------------------------------------------------
  # Crear primera predicción con las frecuencias como prob.:
  # -----------------------------------------------------------------------------
  stopifnot(k %in% 1:3, reg >= 0)

  # Frecuencias de classes en trainset:
  if("prob" %in% names(trainset))
  {
    # Con regularización de las probs estimadas (para no sobre-estimar los que tienen poco):
    setkeyv(trainset, c("Class", "Gene"))
    # probs_Gene <- trainset[, .(prob = (prob)/(.N + reg) ), by = Gene]
    trainset[, probGene := (prob)/(.N + reg), by = c("Class", "Gene")]
    setkeyv(trainset, c("Class", "Variation"))
    # probs_Variation <- trainset[, .(prob = (prob)/(.N + reg) ), by = Variation]
    trainset[, probVariation := (prob)/(.N + reg), by = c("Class", "Variation")]
  } else
  {
    # Con regularización de las probs estimadas (para no sobre-estimar los que tienen poco):
    setkeyv(trainset, c("Class", "Gene"))
    # probs_Gene <- trainset[, .(prob = (prob)/(.N + reg) ), by = Gene]
    trainset[, probGene := (1)/(.N + reg), by = c("Class", "Gene")]
    setkeyv(trainset, c("Class", "Variation"))
    # probs_Variation <- trainset[, .(prob = (prob)/(.N + reg) ), by = Variation]
    trainset[, probVariation := (1)/(.N + reg), by = c("Class", "Variation")]
  }
  # summary(trainset)
  
  # ------------------
  # 2.- Leer testset:
  # ------------------
  testset <- fread(file.path(s_input_path, "test_variants"))
  
  # head(testset)
  # str(testset)
  # summary(testset)
  
  # Proporciones globales (para los NA):
  # Con regularización de las probs estimadas (para no sobre-estimar a los que tienen poco):
  setkey(trainset, Class)
  mean_reg <- function(vv, reg)  return( (sum(vv)) / (length(vv) + reg) )
  for(miclase in 1:9)
  {
    fld <- paste0("prob_class", miclase)
    trainclass <- trainset[Class==miclase,]
    setkey(testset, Gene)
    if(k == 1) prob_global_Gene_en_test <- mean_reg(trainclass$probGene[trainclass$Gene %in% unique(testset$Gene)], reg)
    if(k == 2) prob_global_Gene_no_en_test <- mean_reg(trainclass$probGene[!(trainclass$Gene %in% unique(testset$Gene))], reg)
    if(k == 3) prob_global_Gene <- mean_reg(trainclass$probGene, reg)
    setkey(testset, Variation)
    if(k == 1) prob_global_Variation_en_test <- mean_reg(trainclass$probVariation[trainclass$Variation %in% unique(testset$Variation)], reg)
    if(k == 2) prob_global_Variation_no_en_test <- mean_reg(trainclass$probVariation[!(trainclass$Variation %in% unique(testset$Variation))], reg)
    if(k == 3) prob_global_Variation <- mean_reg(trainclass$probVariation, reg)
    
    if(k == 1) prob_na_Gene <- prob_global_Gene_en_test
    if(k == 2) prob_na_Gene <- prob_global_Gene_no_en_test
    if(k == 3) prob_na_Gene <- prob_global_Gene
    if(k == 1) prob_na_Variation <- prob_global_Variation_en_test
    if(k == 2) prob_na_Variation <- prob_global_Variation_no_en_test
    if(k == 3) prob_na_Variation <- prob_global_Variation
    if(is.na(prob_na_Gene)) prob_na_Gene <- 0
    if(is.na(prob_na_Variation)) prob_na_Variation <- 0
    # prob_na <- mean(c(prob_na_Gene, prob_na_Variation))
    print(paste(fld, "- k", k, "- Gene:", jjformat(mean(trainclass$probGene),2), jjformat(prob_na_Gene,2), "- Variation:", jjformat(mean(trainclass$probVariation),2), jjformat(prob_na_Variation,2)))

    # Quitamos los campos, si existen, para volverlos a meter [merge()]:
    if("probGene" %in% colnames(testset)) testset[, probGene := NULL]
    if("probVariation" %in% colnames(testset)) testset[, probVariation := NULL]
    # Añadir las probs de los Gene:
    setkey(testset, Gene)
    testset <- merge(testset, unique(trainclass[,c("Gene", "probGene")]), all.x = T, by = "Gene")
    setkey(testset, Variation)
    testset <- merge(testset, unique(trainclass[,c("Variation", "probVariation")]), all.x = T, by = "Variation")
    # Predecir las probs de los que no están (i.e. los NAs) con "prob_na":
    # testset[is.na(testset$probGene), probGene := prob_na_Gene]
    # testset[is.na(testset$probVariation), probVariation := prob_na_Variation]
    ## Ligeramente más rápido:
    mi_j <- which(names(testset)=="probGene")
    set(testset, which(is.na(testset$probGene)), mi_j, prob_na_Gene)
    mi_j <- which(names(testset)=="probVariation")
    set(testset, which(is.na(testset$probVariation)), mi_j, prob_na_Variation)
    testset[, (fld) := (probGene + probVariation)/2]
    # print(summary(testset[,fld,with=F]))
    # Quitamos los campos:
    testset[, probGene := NULL]
    testset[, probVariation := NULL]
  }
  print(summary(testset))
  # # Free memory:
  # rm(trainset)
  # gc() # Garbage collector
  
  guardar_submit(testset = testset, fichero = paste0("submitset_k", k, "_Reg", reg, i_sFileSuffix,".csv"), s_output_path = s_output_path)
}

mi_mlogloss <- function(tr_valid, num_classes = 9, b_restore_key = FALSE, b_verbose = 0) # b_verbose=0,1,2
{
  # install.packages("ModelMetrics")
  library(ModelMetrics)
  
  stopifnot(all(c("Class", paste0("prob_class", 1:num_classes)) %in% colnames(tr_valid)))
  if(nrow(tr_valid) == 0)  return(999999)
  if(b_verbose == 2)  print('Calculando mlogloss...')
  mi_tiempo <- system.time({
    eps <- 1e-15
    n <- nrow(tr_valid)
    label_mat <- matrix(0, nrow = n, ncol = num_classes)
    for (i in 1:n)  label_mat[i, tr_valid$Class[i]] <- 1
    pred_mat <- as.matrix(tr_valid[,paste0('prob_class', 1:9),with=F])
    # sum(pred_mat > 1 - eps); sum(pred_mat < eps)
    pred_mat <- pmax(pmin(pred_mat, 1 - eps), eps)
    ret_val <- (-1/n) * sum(label_mat * log(pred_mat)) # (N x 9) * (N x 9) -> (N x 9)
  })
  if(b_verbose == 2)  print(mi_tiempo['elapsed'])
  if(b_verbose >= 1)  print(paste0('mi_mlogloss = ', ret_val))
  return(ret_val)
}
# ----------------------------------------------------------------
# DIVIDE Y VENCERÁS:
# ----------------------------------------------------------------
# IDEA: Split Train in small samples so that we can predict in Testset by several sampled batches (y de paso los comparamos)
mi_split_fich <- function(miNombre='full_trainset', NUM_BLOQUES=16, s_path=s_input_path)
{
  s_fich <- paste0(miNombre, "_", str_pad(NUM_BLOQUES, 3, "left" ,"0"), ".RData")
  if(exists(file.path(s_path, s_fich))) return(0);
  stopifnot(exists(miNombre)) # objeto con nombre miNombre (full_trainset / testset)
  dataset <- get(miNombre)
  setkey(dataset, ID)
  disps <- unique(dataset[,.(ID)], by = "ID") # data.table sin duplicados por clave (y de una columna)
  index <- 1:nrow(disps)
  # Primero ordenamos aleatoriamente:
  if(miNombre != 'testset') # NOTA: testset se deja con el orden original para predecir y crear el submitset en orden.
  {
    index <- sample(index, nrow(disps))
    disps <- disps[index,]
  }
  # Ok.
  rows_per_block <- 1 + as.integer(nrow(disps) / NUM_BLOQUES)
  for(i in 1:NUM_BLOQUES)
  {
    s_fich <- paste0(miNombre, "_", str_pad(i, 3, "left" ,"0"), ".RData")
    index_from <- 1 + (i-1) * rows_per_block
    index_to <- i * rows_per_block
    if(index_to > nrow(disps)) index_to <- nrow(disps)
    disps_set <- disps[index_from:index_to,]
    dataset_batch <- merge(dataset, disps_set, by = 'ID')
    save(dataset_batch, file = file.path(s_path, s_fich)) # full_trainset_001 a full_trainset_016
    rm(dataset_batch)
    gc()
    print(paste0(s_fich, " Ok."))
  }
}
mi_split_train <- function(NUM_BLOQUES=16, s_path=s_input_path)
{ return(mi_split_fich(miNombre='full_trainset', NUM_BLOQUES=NUM_BLOQUES, s_path=s_path)) }
mi_split_test <- function(NUM_BLOQUES=16, s_path=s_input_path)
{ return(mi_split_fich(miNombre='testset', NUM_BLOQUES=NUM_BLOQUES, s_path=s_path)) }

reducir_trainset <- function(mi_set, id_fld="display_id", n_seed = 1, n_porc = 0.1) # 10% por defecto
{
  # En este caso, sencillamente reducimos trainset y devolvemos ambos [lista de dos data.table (trainset y validset)]:
  if(n_porc == 1) return(list(mi_set, data.table()))
  stopifnot(id_fld %in% colnames(mi_set))
  # Reducimos todo para hacer pruebas más rápido:
  set.seed(n_seed)
  setkeyv(mi_set, id_fld)
  ids <- unique(mi_set, by = id_fld)[,id_fld,with=F] # data.table sin duplicados por clave (y de una columna)
  index <- sample(1:nrow(ids), trunc(n_porc * nrow(ids))) # Random sample of ids
  
  return(list(merge(mi_set, ids[index,],  by = id_fld),  # trainset
              merge(mi_set, ids[-index,], by = id_fld))) # validset
  # Ok.
}
# ----------------------------------------------------------------
get_batch_train_filename <- function(numBatch)
{
  return(paste0("full_trainset_", str_pad(numBatch, 3, "left" ,"0"), ".RData"))
}
get_batch_test_filename  <- function(numBatch, numSubBatch = 0)
{
  if(numSubBatch == 0){
    return(paste0("testset_",       str_pad(numBatch, 3, "left" ,"0"), ".RData"))
  } else {
    return(paste0("testset_", str_pad(numBatch, 3, "left" ,"0"), "_", str_pad(numSubBatch, 3, "left" ,"0"), ".RData"))
  }
}
leer_batch_train <- function(numBatch, s_descr = "", s_input_path, i_bVerbose = TRUE)
{
  if(i_bVerbose)  print(paste0(ifelse(s_descr != "", paste0(s_descr, ' - '), ""), 'Batch trainset ', numBatch))
  s_fich_train <- get_batch_train_filename(numBatch)
  stopifnot(file.exists(file.path(s_input_path, s_fich_train)))
  if(exists("full_trainset")) {rm(full_trainset, inherits = TRUE); gc()}
  load(file.path(s_input_path, s_fich_train), .GlobalEnv) # Cuidado que se puede llamar dataset_batch!!!
  if(exists("dataset_batch")) {full_trainset <- dataset_batch; rm(dataset_batch, inherits = TRUE); gc()}
  if(i_bVerbose)  print(paste0(ifelse(s_descr != "", paste0(s_descr, ' - '), ""), 'Batch trainset ', numBatch, ' Ok. ', nrow(full_trainset), ' regs.'))
  return(full_trainset)
}
leer_batch_test <- function(numBatch, s_descr = "", s_input_path, i_bVerbose = TRUE, numSubBatch = 0)
{
  if(i_bVerbose)  print(paste0(ifelse(s_descr != "", paste0(s_descr, ' - '), ""), 'Batch testset ', numBatch, ifelse(numSubBatch == 0, '', paste0('_', numSubBatch))))
  s_fich_test <- get_batch_test_filename(numBatch, numSubBatch)
  stopifnot(file.exists(file.path(s_input_path, s_fich_test)))
  if(exists("testset")) {rm(testset, inherits = TRUE); gc()}
  load(file.path(s_input_path, s_fich_test), .GlobalEnv) # Cuidado que se puede llamar dataset_batch!!!
  if(exists("dataset_batch")) {testset <- dataset_batch; rm(dataset_batch, inherits = TRUE); gc()}
  if(i_bVerbose)  print(paste0(ifelse(s_descr != "", paste0(s_descr, ' - '), ""), 'Batch testset ', numBatch, ifelse(numSubBatch == 0, '', paste0('_', numSubBatch)), ' Ok. ', nrow(testset), ' regs.'))
  return(testset)
}
# ----------------------------------------------------------------
crearTextFeatures <- function(dataset_txt = NULL, # LDA (no supervisado) con textos
                              b_toLower = FALSE,
                              batch_size = 500,
                              ListaLDAfit = NULL,
                              b_entrenar = TRUE, mi_iter.max = 5000, train_repeticiones = 5, # = parallel::detectCores() - 1,
                              palabs_num_features = 27, # = 0 para no usar palabs
                              frases_num_features = 54, # = 0 para no usar frases
                              caract_num_features = 72, # = 0 para no usar chars
                              mis_ngrams = list(
                                frases_n_grams = 1:3,  frases_skip_grams = 0:1,  # Hay menos frases, así que genero más n-gramas
                                palabs_n_grams = 2:3,  palabs_skip_grams = 0:2,
                                caract_n_grams = 1:5,  caract_skip_grams = 0 # Idealmente 1:20... (avg.wrd.len=5.4+3.3; avg.wrd.len.sin.stopwords=6.7+3.3)
                              ),
                              txtWhat = '', # Descriptivo para logs (en el caso b_entrenar = FALSE)
                              s_Fic_log = 'crearTextFeatures.log') # La función añade las features al data.table!
{
  if(b_entrenar && is.null(ListaLDAfit))
    ListaLDAfit <- list(frases=NULL, palabs=NULL, caract=NULL)
  stopifnot(palabs_num_features + frases_num_features + caract_num_features != 0) # O bien palabs. o bien frases, o bien chars, o bien dos de ellos o los tres
  if(txtWhat != '')
    txtWhat <- paste0('[', txtWhat, '] ')
  jjprint('---------------------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE) # Línea separadora para el inicio de este log...
  jjprint(paste0(txtWhat, 'crearTextFeatures(b_entrenar:', b_entrenar, ') - Inicio'), logfile = s_Fic_log)
  jjprint('---------------------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
  # # G_b_DEBUG <- TRUE # Para pruebas...
  # if(G_b_DEBUG){
  #   s_Fic_log <- NULL
  #   b_toLower = FALSE; batch_size = 5000; ListaLDAfit = list(frases=NULL, palabs=NULL)
  #   b_entrenar = TRUE; mi_iter.max = 5; palabs_num_features = 27; frases_num_features = 54
  #   
  #   if(!exists("dataset_txt")) mi_load_uniq_txt('train', s_Fic_log = s_Fic_log) # mi_load(s_input_path, "training_text.RData") # Para pruebas...
  #   if(!exists("dataset_txt")) dataset_txt <- full_trainset_txt # Para pruebas...
  #   if(is.null(dataset_txt)) dataset_txt <- full_trainset_txt # Para pruebas...
  # }
  
  # Reducimos a los unique usando la clave train_texto_id/test_texto_id:
  pk <- ifelse("train_texto_id" %in% colnames(dataset_txt), "train_texto_id",
               ifelse("test_texto_id" %in% colnames(dataset_txt), "test_texto_id",
                      ""))
  stopifnot(pk %in% c("train_texto_id", "test_texto_id"))
  # setkeyv(dataset_txt, pk)
  # nrow(unique(dataset_txt[, c(pk, "texto"), with=F]))
  jjprint('Extraemos textos...', logfile = s_Fic_log)
  
  if(G_b_DEBUG){
    jjprint('NOTA: Solo los primeros N docs...', logfile = s_Fic_log)
    mis_textos <- unique(dataset_txt[1:40, c(pk, "texto"), with=F])
  } else {
    mis_textos <- unique(dataset_txt[, c(pk, "texto"), with=F])
  }
  
  if(b_toLower) {
    jjprint('---------------------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
    mi_tiempo <- system.time({
      jjprint('char_tolower()...', logfile = s_Fic_log)
      mis_textos[, texto := char_tolower(texto, keep_acronyms=TRUE)]
    })
    mi_tiempo_tolower <- paste0('Tiempo tolower(): ', jjfmt(mi_tiempo['elapsed']/60, 3), ' minutos')
    jjprint(paste0(mi_tiempo_tolower), logfile = s_Fic_log)
    jjprint('---------------------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
  }
  n_tot <- nrow(mis_textos)
  byseq <- batch_size
  if(n_tot %% byseq == 1) # Misterio de LDA(): NO FUNCIONA CON UNA ÚNICA FILA
  {
    sTmp <- paste0('ERROR! Misterios de la cripta: n_tot(', n_tot, ') MOD byseq(', byseq, ') == 1. [LDA(): NO FUNCIONA CON UNA ÚNICA FILA!]')
    jjprint(sTmp, logfile = s_Fic_log)
    stop(sTmp)
  }
  jjprint(paste0(txtWhat, 'Tokenizing ', n_tot, ' docs (batch_size=', byseq,')...'), logfile = s_Fic_log)
  Vector_num_features <- c()
  mis_features <- list(palabs=list(0), frases=list(0), caract=list(0)) # mis_features[['palabs']][[length(mis_features[['palabs']])+1]] <- obj # añade obj a la lista
  dfStats <- NULL
  for(i1 in seq(1,n_tot,by=byseq)) { # Uno por uno (batch_size) por problemas de memoria...
    # i1 = 1 # PRUEBAS
    i2 <- min(c(i1+byseq-1 , n_tot))
    jjprint(paste0(txtWhat, 'Procesando... ', '[', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
    mis_textos_tmp <- mis_textos[i1:i2,]
    corp <- corpus(mis_textos_tmp, text_field = "texto", docid_field = pk)
    # summary(corp, 5) # Ver solo los 10 primeros docs
    # # full_corp <- corpus(dataset_txt[, c("ID", "texto"), with=F], text_field = "texto", docid_field = "ID")
    # # summary(full_corp, 5)
    # # # options(width = 80)
    # # # kwic(corp, "cancer")
    ListaDims_antes <- list(frases=c(0,0), palabs=c(0,0), caract=c(0,0))
    ListaDims_desp <- list(frases=c(0,0), palabs=c(0,0), caract=c(0,0))
    ListaPerplex_train <- list(frases=0, palabs=0, caract=0)
    ListaPerplex_valid <- list(frases=0, palabs=0, caract=0)
    if(frases_num_features != 0)
    {
      jjprint('-----------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
      numfeat <- ifelse(b_entrenar, frases_num_features, ListaLDAfit[['frases']]@k)
      Vector_num_features <- c(Vector_num_features, frases = numfeat) # frases_num_features <- 54
      # rm(dataset_txt, mis_textos, frases); gc(); save.image(paste0(s_input_path, 'borrar.RData'))
      # load(paste0(s_input_path, 'borrar.RData')); gc()
      n.grams <- mis_ngrams[["frases_n_grams"]] # 1:3 # 1
      n.skip <- mis_ngrams[["frases_skip_grams"]] # 0:1 # 0
      frases <- tokens(corp, remove_punct = TRUE, remove_symbols = TRUE
                       , what = "sentence"
                       , remove_numbers = FALSE  # En frases nunca quitamos números
                       , remove_separators = FALSE  # idem
                       , ngrams = n.grams
                       , skip = n.skip  # cercanía == desde 0 hasta N frases más allá (i.e. 3 n-gramas por cada frase y cada n-grama)
                       , verbose = TRUE
      )
      {
        jjprint(paste0("frases: ANTES [n.grams=(", paste(n.grams, collapse = ','), ')] docs: ', nrow(corp$documents), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        jjprint(paste0('Número total de tokens: ', jjfmt(as.integer(sum(sapply(frases, length)))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Número total de tokens
        ListaDims_antes[['frases']] <- c(nrow(corp$documents), as.integer(sum(sapply(frases, length))))

        # jjprint(paste0('frases: Quitando stopwords...[', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        # frases <- tokens_remove(frases, features = stopwords("english"), verbose=T) # En frases no hace nada (CREO!)
        jjprint(paste0('frases: Quitando algunas frases < 20 chars o con < 6 palabras más...[', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        # Error: cannot allocate vector of size 171.9 Gb: frases <- tokens_remove(frases, features = c("fig\\.", "FIG\\.", "[1-9]\\. fig\\."), valuetype = c("regex"), case_insensitive = TRUE, verbose=TRUE)
        p=as.list(unique(unlist(sapply(frases, function(x) { return(x[str_length(x) < 20 | str_count(x, " ") < 5]) }))))
        if(length(p) != 0)
        {
          jjprint(paste0('frases: Quitamos ', jjfmt(length(p)), ' unique tokens...[', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          # jjprint(paste0('frases: tokens antes: ', jjfmt(as.integer(sum(sapply(frases, length))))), logfile = s_Fic_log)
          frases <- tokens_remove(frases, features = p, valuetype = c("fixed"), verbose=TRUE)
          # jjprint(paste0('frases: tokens desp.: ', jjfmt(as.integer(sum(sapply(frases, length))))), logfile = s_Fic_log)
        }
        jjprint(paste0("frases: DESP. [n.grams=(", paste(n.grams, collapse = ','), ')] docs: ', nrow(corp$documents), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        jjprint(paste0('Número total de tokens: ', jjfmt(as.integer(sum(sapply(frases, length)))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Número total de tokens
        # jjprint(paste0('Doc con más tokens: ', jjfmt(as.integer(which.max(sapply(frases, length)))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Doc con más tokens
        jjprint(paste0('Máx. tokens: ', jjfmt(as.integer(sapply(frases, length)[which.max(sapply(frases, length))])), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Cuántos tokens tiene
        # jjprint(paste0('Doc con menos tokens (pero más de uno): ', jjfmt(as.integer(which.min(sapply(frases, length)[sapply(frases, length)>1]))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Doc con menos tokens (pero más de uno)
        jjprint(paste0('Mín. tokens (>1): ', jjfmt(as.integer(sapply(frases, length)[sapply(frases, length)>1][which.min(sapply(frases, length)[sapply(frases, length)>1])])), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Cuántos tokens tiene
        # jjprint(frases[sapply(frases, length)>1][which.min(sapply(frases, length)[sapply(frases, length)>1])], logfile = s_Fic_log, b_add_datetime = FALSE) # Los tokens
      }
    } else {
      Vector_num_features <- c(Vector_num_features, frases = 0)
    }
    if(palabs_num_features != 0)
    {
      jjprint('-----------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
      numfeat <- ifelse(b_entrenar, palabs_num_features, ListaLDAfit[['palabs']]@k)
      # rm(dataset_txt, mis_textos, palabs); gc(); save.image(paste0(s_input_path, 'borrar.RData'))
      # load(paste0(s_input_path, 'borrar.RData')); gc()
      Vector_num_features <- c(Vector_num_features, palabs = numfeat) # palabs_num_features <- 27
      n.grams <- mis_ngrams[["palabs_n_grams"]] # 2:3 # n.grams <- 1
      n.skip <- mis_ngrams[["palabs_skip_grams"]] # 0:2 # 0
      palabs <- tokens(corp, remove_punct = TRUE, remove_symbols = TRUE
                       , what = "word"  # , what = "fasterword", remove_url = TRUE
                       , remove_numbers = all(n.grams==1)  # Quitamos números salvo para n-grams
                       , remove_separators = TRUE
                       , ngrams = 1 # Cf. más abajo (n.grams)
                       , skip = 0 # Cf. más abajo (n.skip) # cercanía == desde 0 hasta 2 palabras más allá (i.e. 3 n-gramas por cada palabra y cada n-grama)
                       , verbose = TRUE
      )
      {
        jjprint(paste0("palabs: ANTES [n.grams=(", paste(n.grams, collapse = ','), ')] docs: ', nrow(corp$documents), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        jjprint(paste0('Número total de tokens: ', jjfmt(as.integer(sum(sapply(palabs, length)))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Número total de tokens
        ListaDims_antes[['palabs']] <- c(nrow(corp$documents), as.integer(sum(sapply(palabs, length))))
        
        jjprint(paste0('palabs: Quitando stopwords...[', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        palabs <- tokens_remove(palabs, features = stopwords("english"), verbose=T)
        jjprint(paste0("palabs: ANTES [n.grams=(", paste(n.grams, collapse = ','), ')] docs: ', nrow(corp$documents), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        jjprint(paste0('Número total de tokens: ', jjfmt(as.integer(sum(sapply(palabs, length)))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Número total de tokens
        
        if(any(n.grams!=1)) # | any(n.skip!=0))
        {
          jjprint(paste0('palabs: Finalmente, creamos n-grams...[', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          if(b_entrenar)
          {
            tmp_dfm <- dfm(x = palabs, tolower = FALSE, verbose = TRUE)
            jjprint(paste0('palabs: Ok. (dfm) = ', jjfmt(ncol(tmp_dfm)), ' tokens.', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
            jjprint(paste0('palabs: Reducimos matriz antes de n.grams (dfm): ',
                           'Quitamos tokens con freqs <= log(numdocs) y docfreqs <= log(numdocs) y docfreqs >= 0.9*numdocs',
                           '[freqs <= ', jjfmt(log(nrow(tmp_dfm)),1), '] y ',
                           '[docfreqs <= ', jjfmt(log(nrow(tmp_dfm)),1), '] y ',
                           '[docfreqs >= ', jjfmt(0.9 * nrow(tmp_dfm),1), ']',
                           '...', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
            if(log(nrow(tmp_dfm)) > 1) {
              tmp_dfm <- dfm_trim(x = tmp_dfm
                                  , min_count   = log(nrow(tmp_dfm))
                                  , min_docfreq = log(nrow(tmp_dfm))
                                  , max_docfreq = 0.9 * nrow(tmp_dfm)
                                  , verbose = T)
            } else {
              tmp_dfm <- dfm_trim(x = tmp_dfm
                                  , max_docfreq = 0.9 * nrow(tmp_dfm)
                                  , verbose = T)
            }
            jjprint(paste0('palabs: Ok. (dfm) = ', jjfmt(ncol(tmp_dfm)), ' tokens.', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
            
            palabs <- tokens_select(x = palabs, features = featnames(tmp_dfm))
            rm(tmp_dfm); gc();
            # save.image(paste0(s_input_path, 'borrar2.RData'))
            # load(paste0(s_input_path, 'borrar2.RData')); gc()
          }
          palabs <- tokens_ngrams(x = palabs, n = n.grams, skip = n.skip)
        }
        # jjprint(paste0('palabs: Quitando algunos tokens más...[', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        # tokens_remove(palabs, features = c("fig\\.", "FIG\\.", "[1-9]\\. fig\\."), valuetype = c("regex"), case_insensitive = TRUE, verbose=TRUE)
        # # p=unique(unlist(sapply(palabs, function(x) { return(x[str_length(x) < 20 | str_count(x, " ") < 5]) })))
        # # palabs <- tokens_remove(palabs, features = p, valuetype = "fixed", verbose=T)
        
        jjprint(paste0("palabs: DESP. [n.grams=(", paste(n.grams, collapse = ','), ')] docs: ', nrow(corp$documents), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        jjprint(paste0('Número total de tokens: ', jjfmt(as.integer(sum(sapply(palabs, length)))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Número total de tokens
        # jjprint(paste0('Doc con más tokens: ', jjfmt(as.integer(which.max(sapply(palabs, length)))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Doc con más tokens
        jjprint(paste0('Máx. tokens: ', jjfmt(as.integer(sapply(palabs, length)[which.max(sapply(palabs, length))])), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Cuántos tokens tiene
        # jjprint(paste0('Doc con menos tokens (pero más de uno): ', jjfmt(as.integer(which.min(sapply(palabs, length)[sapply(palabs, length)>1]))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Doc con menos tokens (pero más de uno)
        jjprint(paste0('Mín. tokens (>1): ', jjfmt(as.integer(sapply(palabs, length)[sapply(palabs, length)>1][which.min(sapply(palabs, length)[sapply(palabs, length)>1])])), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Cuántos tokens tiene
        # jjprint(palabs[sapply(palabs, length)>1][which.min(sapply(palabs, length)[sapply(palabs, length)>1])], logfile = s_Fic_log, b_add_datetime = FALSE) # Los tokens
      }
    } else {
      Vector_num_features <- c(Vector_num_features, palabs = 0)
    }
    
    if(caract_num_features != 0)
    {
      jjprint('-----------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
      numfeat <- ifelse(b_entrenar, caract_num_features, ListaLDAfit[['caract']]@k)
      # rm(dataset_txt, mis_textos, caract); gc(); save.image(paste0(s_input_path, 'borrar.RData'))
      # load(paste0(s_input_path, 'borrar.RData')); gc()
      Vector_num_features <- c(Vector_num_features, caract = numfeat) # caract_num_features <- 27
      n.grams <- mis_ngrams[["caract_n_grams"]] # 1:20 # n.grams <- 1
      n.skip <- mis_ngrams[["caract_skip_grams"]] # n.skip <- 0 # Los caracteres tienen un orden preciso, no tiene sentido buscar ngramas "cercanos"...
      caract <- tokens(corp, remove_punct = TRUE, remove_symbols = TRUE
                       , what = "character"
                       , remove_numbers = all(n.grams==1)  # Quitamos números salvo para n-grams
                       , remove_separators = TRUE
                       , ngrams = 1 # Cf. más abajo (n.grams)
                       , skip = 0 # Cf. más abajo (n.skip) # cercanía == desde 0 hasta 2 palabras más allá (i.e. 3 n-gramas por cada palabra y cada n-grama)
                       , verbose = TRUE
      )
      {
        jjprint(paste0("caract: ANTES [n.grams=(", paste(n.grams, collapse = ','), ')] docs: ', nrow(corp$documents), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        jjprint(paste0('Número total de tokens: ', jjfmt(as.integer(sum(sapply(caract, length)))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Número total de tokens
        ListaDims_antes[['caract']] <- c(nrow(corp$documents), as.integer(sum(sapply(caract, length))))

        jjprint(paste0('caract: Quitando stopwords...[', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        caract <- tokens_remove(caract, features = stopwords("english"), verbose=T)
        jjprint(paste0("caract: ANTES [n.grams=(", paste(n.grams, collapse = ','), ')] docs: ', nrow(corp$documents), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        jjprint(paste0('Número total de tokens: ', jjfmt(as.integer(sum(sapply(caract, length)))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Número total de tokens

        if(any(n.grams!=1)) # | any(n.skip!=0))
        {
          jjprint(paste0('caract: Finalmente, creamos n-grams...[', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          if(b_entrenar)
          {
            tmp_dfm <- dfm(x = caract, tolower = FALSE, verbose = TRUE)
            jjprint(paste0('caract: Ok. (dfm) = ', jjfmt(ncol(tmp_dfm)), ' tokens.', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
            jjprint(paste0('caract: Reducimos matriz antes de n.grams (dfm): ',
                           'Quitamos tokens con freqs <= log(numdocs) y docfreqs <= log(numdocs) y docfreqs >= 0.9*numdocs',
                           '[freqs <= ', jjfmt(log(nrow(tmp_dfm)),1), '] y ',
                           '[docfreqs <= ', jjfmt(log(nrow(tmp_dfm)),1), '] y ',
                           '[docfreqs >= ', jjfmt(0.9 * nrow(tmp_dfm),1), ']',
                           '...', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
            if(log(nrow(tmp_dfm)) > 1) {
              tmp_dfm <- dfm_trim(x = tmp_dfm
                                  , min_count   = log(nrow(tmp_dfm))
                                  , min_docfreq = log(nrow(tmp_dfm))
                                  , max_docfreq = 0.9 * nrow(tmp_dfm)
                                  , verbose = T)
            } else {
              tmp_dfm <- dfm_trim(x = tmp_dfm
                                  , max_docfreq = 0.9 * nrow(tmp_dfm)
                                  , verbose = T)
            }
            jjprint(paste0('caract: Ok. (dfm) = ', jjfmt(ncol(tmp_dfm)), ' tokens.', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)

            caract <- tokens_select(x = caract, features = featnames(tmp_dfm))
            rm(tmp_dfm); gc();
            # save.image(paste0(s_input_path, 'borrar2.RData'))
            # load(paste0(s_input_path, 'borrar2.RData')); gc()
          }
          caract <- tokens_ngrams(x = caract, n = n.grams, skip = n.skip)
        }
        # jjprint(paste0('caract: Quitando algunos tokens más...[', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        # tokens_remove(caract, features = c("fig\\.", "FIG\\.", "[1-9]\\. fig\\."), valuetype = c("regex"), case_insensitive = TRUE, verbose=TRUE)
        # # p=unique(unlist(sapply(caract, function(x) { return(x[str_length(x) < 20 | str_count(x, " ") < 5]) })))
        # # caract <- tokens_remove(caract, features = p, valuetype = "fixed", verbose=T)

        jjprint(paste0("caract: DESP. [n.grams=(", paste(n.grams, collapse = ','), ')] docs: ', nrow(corp$documents), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        jjprint(paste0('Número total de tokens: ', jjfmt(as.integer(sum(sapply(caract, length)))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Número total de tokens
        # jjprint(paste0('Doc con más tokens: ', jjfmt(as.integer(which.max(sapply(caract, length)))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Doc con más tokens
        jjprint(paste0('Máx. tokens: ', jjfmt(as.integer(sapply(caract, length)[which.max(sapply(caract, length))])), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Cuántos tokens tiene
        # jjprint(paste0('Doc con menos tokens (pero más de uno): ', jjfmt(as.integer(which.min(sapply(caract, length)[sapply(caract, length)>1]))), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Doc con menos tokens (pero más de uno)
        jjprint(paste0('Mín. tokens (>1): ', jjfmt(as.integer(sapply(caract, length)[sapply(caract, length)>1][which.min(sapply(caract, length)[sapply(caract, length)>1])])), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log) # Cuántos tokens tiene
        # jjprint(caract[sapply(caract, length)>1][which.min(sapply(caract, length)[sapply(caract, length)>1])], logfile = s_Fic_log, b_add_datetime = FALSE) # Los tokens
      }
    } else {
      Vector_num_features <- c(Vector_num_features, caract = 0)
    }
    
    # PRUEBAS }
    # hist(sapply(frases, length))
    # hist(sapply(palabs, length))
    # hist(sapply(caract, length))
    
    if(G_b_DEBUG){
      pp <- c(text='hola me llamo pepito de los 50 palotes',
              'adiós te llamas joselito de los 25 árboles')
      # palabs <- tokens(pp, remove_punct = TRUE, remove_symbols = TRUE
      #                  , remove_numbers = all(n.grams==1) # Quitamos números salvo para n-grams
      #                  , remove_separators = TRUE
      #                  , what = "word"
      #                  # , what = "fasterword", remove_url = TRUE
      #                  , ngrams = n.grams # ngrams = 2:3
      #                  , skip = 0:2
      #                  , verbose = TRUE)
    }
    if(G_b_DEBUG){
      s_Fic_log <- NULL
      rm(caract, tmp_dfm, pp); gc()
      stopifnot(mi_load_uniq_txt('train', s_Fic_log = NULL))
      pp <- corpus(dataset_txt, text_field = "texto", docid_field = "train_texto_id")
      rm(dataset_txt); gc()
      caract <- tokens(pp, remove_punct = TRUE, remove_symbols = TRUE # Quitamos todo porque dan problemas (bad coding of quanteda::dfm(), tokens(), etc...))
                       , remove_separators = TRUE, remove_numbers = FALSE # Pero no quitamos números
                       , what = "character"
                       , ngrams = 1 # n.grams después!!!
                       , skip = 0
                       , verbose = TRUE)
      # c(sum(str_length(pp)), length(caract$text), sum(str_length(caract$text)))
      jjfmt(c(sum(str_length(pp)), sum(sapply(caract, length))))
      tmp_dfm <- dfm(x = caract, tolower = FALSE, verbose = TRUE)
      tmp_dfm <- dfm_trim(x = tmp_dfm
                          , min_count   = log(nrow(tmp_dfm))
                          , min_docfreq = log(nrow(tmp_dfm))
                          , max_docfreq = 0.9 * nrow(tmp_dfm)
                          , verbose = T)
      caract <- tokens_select(x = caract, features = featnames(tmp_dfm))
      jjfmt(c(sum(str_length(pp)), sum(sapply(caract, length))))
      rm(tmp_dfm); gc();
      n.grams = 1:20; n.skip = 0
      caract <- tokens_ngrams(x = caract, n = n.grams, skip = n.skip)
      jjfmt(c(sum(str_length(pp)), sum(sapply(caract, length))))
    }
    for(tok_name in c("frases", "palabs", "caract"))
    { # tok_name = "frases"
      num_features <- Vector_num_features[tok_name]
      if(num_features != 0)
      {
        tok_descr <- paste0(txtWhat, '[',tok_name,'] ')
        jjprint('-----------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
        
        mis_tokens <- get(tok_name) # frases / palabs / caract
        myLDAfit <- ListaLDAfit[[tok_name]] # ListaLDAfit[["frases"]] / ListaLDAfit[["palabs"]] / ListaLDAfit[["caract"]]
        
        jjprint(paste0(tok_descr, 'Creando matriz (dfm)...', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        quantdfm <- dfm(mis_tokens, tolower = FALSE, verbose = TRUE)
        jjprint(paste0(tok_descr, paste0(jjfmt(dim(quantdfm)), c(" regs", " cols"), collapse = "/")), logfile = s_Fic_log)
        if(b_entrenar)
        {
          rm(mis_tokens); gc()
          # feats_ok <- featnames(quantdfm)[colSums(quantdfm) > log(nrow(quantdfm)/num_features)]
          # quantdfm_reduc <- dfm_select(x = quantdfm, features = feats_ok, valuetype = "fixed", selection = "keep", min_nchar = 5, max_nchar = NULL, verbose = TRUE)
          
          jjprint(paste0(tok_descr, 'Reducimos matriz después de n.grams (dfm): ',
                         'Quitamos tokens con freqs <= log(log(numdocs)) y docfreqs <= log(log(numdocs)) y docfreqs >= 0.9*numdocs',
                         '[freqs <= ', jjfmt(log(log(nrow(quantdfm))),1), '] y ',
                         '[docfreqs <= ', jjfmt(log(log(nrow(quantdfm))),1), '] y ',
                         '[docfreqs >= ', jjfmt(0.9 * nrow(quantdfm),1), ']',
                         '...', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          if(log(log(nrow(quantdfm))) > 1) {
            quantdfm <- dfm_trim(x = quantdfm
                                , min_count   = log(log(nrow(quantdfm)))
                                , min_docfreq = log(log(nrow(quantdfm)))
                                , max_docfreq = 0.9 * nrow(quantdfm)
                                , verbose = T)
          } else {
            quantdfm <- dfm_trim(x = quantdfm
                                 , max_docfreq = 0.9 * nrow(quantdfm)
                                 , verbose = T)
          }
          jjprint(paste0(tok_descr, 'Ok. (dfm) = ', jjfmt(ncol(quantdfm)), ' tokens.', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          
          # feats_ok_idx <- which(as.integer(colSums(quantdfm)) > log(nrow(quantdfm)/num_features))
          # jjprint(paste0(tok_descr, 'Reducimos matriz (dfm): Quitamos ', jjfmt(ncol(quantdfm) - length(feats_ok_idx))
          #                , ' tokens con freqs <= log(numdocs)/k [freqs <= ', jjfmt(log(nrow(quantdfm)/num_features),1), ']...', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          # if(length(feats_ok_idx) != 0 & length(feats_ok_idx) != ncol(quantdfm)) {
          #   quantdfm_reduc <- quantdfm[, feats_ok_idx]
          #   quantdfm <- quantdfm_reduc; rm(quantdfm_reduc)
          # }
          # jjprint(paste0(tok_descr, 'Ok. (dfm) = ', jjfmt(ncol(quantdfm)), ' tokens.', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          
          # feats_ok_idx <- which(as.integer(colSums(quantdfm)) < 0.60 * nrow(quantdfm))
          # jjprint(paste0(tok_descr, 'Reducimos matriz (dfm): Quitamos ', jjfmt(ncol(quantdfm) - length(feats_ok_idx))
          #                , ' tokens con freqs >= 60 % numdocs [freqs >= ', jjfmt(0.60 * nrow(quantdfm),1), ']...', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          # if(length(feats_ok_idx) != 0 & length(feats_ok_idx) != ncol(quantdfm)) {
          #   quantdfm_reduc <- quantdfm[, feats_ok_idx]
          #   quantdfm <- quantdfm_reduc; rm(quantdfm_reduc)
          # }
          # rm(feats_ok_idx); gc()
          # jjprint(paste0(tok_descr, 'Ok. (dfm) = ', jjfmt(ncol(quantdfm)), ' tokens.', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          # # jjprint(paste0(tok_descr, paste0(jjfmt(dim(quantdfm)), c(" regs", " cols"), collapse = "/")), logfile = s_Fic_log)
        }
        
        # rm(frases, palabs); gc(); save.image(paste0(s_input_path, 'borrar2.RData'))
        # load(paste0(s_input_path, 'borrar2.RData')); gc()
        jjprint(paste0(tok_descr, 'convirtiendo a matriz para topicmodels (dtm)...', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
        quantdtm <- convert(quantdfm, to = "topicmodels", docvars = data.frame(doc_id=rownames(quantdfm)))
        jjprint(paste0(tok_descr, paste0(jjfmt(dim(quantdtm)), c(" regs", " cols"), collapse = "/")), logfile = s_Fic_log)
        ListaDims_desp[[tok_name]] <- dim(quantdtm) # ListaDims_desp[['frases']] / ListaDims_desp[['palabs']] / ListaDims_desp[['caract']]
        
        if(b_entrenar)
        {
          
          # rm(quantdfm); gc(); save.image(paste0(s_input_path, 'borrar3.RData'))
          # load(paste0(s_input_path, 'borrar3.RData')); gc(); mi_repe <- 1
          # quantdtm <- quantdtm[sample.int(nrow(quantdtm), nrow(quantdtm)/2), ]
          # rm(corp, mis_tokens); gc(); save.image(paste0(s_input_path, 'borrar4.RData'))
          # load(paste0(s_input_path, 'borrar4.RData')); gc()
          # # quantdtm <- quantdtm[sample.int(nrow(quantdtm), nrow(quantdtm)/2), ]

          # Entrenar:
          NumRepes <- train_repeticiones # = parallel::detectCores() - 1
          jjprint('*********************************************************', logfile = s_Fic_log, b_add_datetime = FALSE)
          jjprint(paste0(tok_descr, 'Entrenando el clasificador (no-supervisado) Latent Dirichlet Allocation (LDA) para ', num_features,' clases...', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          jjprint(paste0(tok_descr, ' - Dim(quantdtm) = (',dim(quantdtm)[1],',',dim(quantdtm)[2],')'), logfile = s_Fic_log)
          jjprint(paste0(tok_descr, ' - NumRepes = ',NumRepes,'; iter.max = ', mi_iter.max), logfile = s_Fic_log)
          # # Usamos myLDAfit anterior como modelo inicial si ya existía:
          # if(!exists("myLDAfit")) myLDAfit_ant <- NULL else myLDAfit_ant <- myLDAfit
          myLDAfit_ant <- myLDAfit # Viene en la llamada
          retL <- foreach(mi_repe = seq_len(NumRepes) # = seq_len(parallel::detectCores() - 1)
                          , .inorder=FALSE, .multicombine = T, .combine = list
                          , .packages = c('topicmodels')
                          , .export = c('jjfmt', 's_input_path')
          ) %do%  #  %dopar%
          {
            mi_tiempo <- system.time({
              myLDAfit <- LDA(x = quantdtm, k = num_features, method = "VEM",
                              model = myLDAfit_ant,
                              control = list(seed = mi_repe * (as.integer(Sys.time()) %% 10000),
                                             save = 1, prefix = paste0(s_input_path, Sys.info()["nodename"], "_LDA_", mi_repe, "_"), # Para %dopar%
                                             # save = 0, # Para %do%
                                             nstart = 1, # Number of repeated random starts
                                             # alpha = 0.02, estimate.alpha = FALSE,
                                             estimate.alpha = TRUE,
                                             em=list(iter.max=mi_iter.max, tol=1e-5),
                                             initialize = ifelse(is.null(myLDAfit_ant), "random", "model"),
                                             verbose = as.integer((500 - 1 + mi_iter.max) / 500) # Print cada N iters
                              ))
            })
            tiempo_lda <- paste0('Tiempo LDA(): ', jjfmt(mi_tiempo['elapsed']/60, 3), ' minutos')
            return(list(myLDAfit, tiempo_lda))
          }
          for(mi_repe in seq_len(NumRepes))
          {
            tok_descr <- paste0(tok_descr, '[', mi_repe, '/', length(retL), '] ')
            if(NumRepes!=1) {
              myLDAfit <- retL[[mi_repe]][[1]]; tiempo_lda <- retL[[mi_repe]][[2]]
            } else {
              myLDAfit <- retL[[1]]; tiempo_lda <- retL[[2]]
            }
            jjprint(paste0(tok_descr, tiempo_lda, ' [', num_features, ' clases]', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
            if(tok_name != 'frases'){
              jjprint(paste0(tok_descr, 'Terms: most likely terms:', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
              jjprint(terms(myLDAfit, k = 3), logfile = s_Fic_log, b_add_datetime = F)
            }
            # jjprint(paste0(tok_descr, 'Topics: most likely topics:', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
            # # jjprint(topics(x = myLDAfit, k = min(15, num_features)), logfile = s_Fic_log, b_add_datetime = F)
            # jjprint(sapply(data.table(posterior(object = myLDAfit)$topics), summary), logfile = s_Fic_log, b_add_datetime = F)
            ListaPerplex_train[[tok_name]] <- ListaPerplex_train[[tok_name]] + perplexity(myLDAfit)
            jjprint(paste0(tok_descr, 'Perplexity: ', jjfmt(perplexity(myLDAfit), 3), ' - alpha: ', jjfmt(myLDAfit@alpha, 3), ' - Topics(k): ', jjfmt(myLDAfit@k), ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
            tok_descr <- paste0(txtWhat, '[',tok_name,'] ')
          }
          jjprint('*********************************************************', logfile = s_Fic_log, b_add_datetime = FALSE)
          ListaPerplex_train[[tok_name]] <- ListaPerplex_train[[tok_name]] / NumRepes # ¡Promedio!
          ListaLDAfit[[tok_name]] <- myLDAfit # ListaLDAfit[["frases"]] / ListaLDAfit[["palabs"]] / ListaLDAfit[["caract"]]
        } else {
          # Clasificar:
          num_features <- myLDAfit@k # El modelo ya entrenado tiene num_topics (k) ya definido
          jjprint(paste0(tok_descr, 'Clasificamos con el modelo LDA entrenado previamente:', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          lda_inf <- posterior(myLDAfit, newdata = quantdtm)
          ListaPerplex_train[[tok_name]] <- perplexity(myLDAfit)
          ListaPerplex_valid[[tok_name]] <- perplexity(myLDAfit, newdata = quantdtm)
          jjprint(paste0(tok_descr, "__Perplex. = ", jjfmt(perplexity(myLDAfit, newdata = quantdtm), 2), '/', jjfmt(perplexity(myLDAfit), 2)), logfile = s_Fic_log)
          jjprint(paste0(tok_descr, 'Preparamos DT...', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          mis_features_tmp <- data.table(lda_inf$topics, keep.rownames = TRUE) # probs de las num_features features (columnas) para cada doc (filas)
          mis_features_tmp[, rn := as.integer(rn)]
          # NOTA: mis_features_tmp$rn es el doc_id (train_texto_id/test_texto_id)
          colnames(mis_features_tmp) <- c(pk, paste0('LDA_', tok_name, num_features, '_', 1:num_features))
          # mis_features_tmp[, (pk) := mis_textos_tmp[, (pk), with=FALSE]] # colnames(mis_features_tmp)[1] <- pk
          jjprint(paste0('Ok. DT de ', nrow(mis_features_tmp), ' regs. y ', ncol(mis_features_tmp), ' cols.'), logfile = s_Fic_log)
          stopifnot(ncol(mis_features_tmp) == num_features + 1)
          jjprint(paste0(tok_descr, 'Insertamos DT en lista [mis_features]...', ' [', jjfmt(100 * i2/n_tot, 2), '%]'), logfile = s_Fic_log)
          mis_features[[tok_name]][[length(mis_features[[tok_name]])+1]] <- mis_features_tmp
        } # Fin de if(b_entrenar)
        tmpStats <- data.frame(
          batchPrcnt = i2/n_tot, # de 0 a 1
          tok_name = tok_name, # 'frases' / 'palabs' / 'caract'
          rows_ini = ListaDims_antes[[tok_name]][1], # Rows/Tokens iniciales (ListaDims_antes[['frases']] / ListaDims_antes[['palabs']])
          toks_ini = ListaDims_antes[[tok_name]][2], # Rows/Tokens iniciales (ListaDims_antes[['frases']] / ListaDims_antes[['palabs']])
          rows_fin = ListaDims_desp[[tok_name]][1],   # Rows/Tokens finales (para entrenar LDA) (ListaDims_desp[['frases']] / ListaDims_desp[['palabs']])
          toks_fin = ListaDims_desp[[tok_name]][2],   # Rows/Tokens finales (para entrenar LDA) (ListaDims_desp[['frases']] / ListaDims_desp[['palabs']])
          model_alpha = ListaLDAfit[[tok_name]]@alpha,
          model_k = ListaLDAfit[[tok_name]]@k,
          Perplex_train = ListaPerplex_train[[tok_name]], # Perplexity del modelo (i.e. con trainset) (ListaPerplex_train[['frases']] / ListaPerplex_train[['palabs']])
          Perplex_valid = ListaPerplex_valid[[tok_name]]  # Perplexity del modelo con validset (ListaPerplex_valid[['frases']] / ListaPerplex_valid[['palabs']])
        )
        if(is.null(dfStats)) {
          dfStats <- tmpStats
        } else {
          dfStats <- rbind(dfStats, tmpStats)
        }
      } # Fin de if(num_features != 0)
    } # Fin de for(tok_name)
  } # Fin de for(i)
  
  if(b_entrenar) {
    for(tok_name in c("frases", "palabs", "caract"))
    { # tok_name = "palabs"
      if(!is.null(ListaLDAfit[[tok_name]]))
      {
        tok_descr <- paste0(txtWhat, '[',tok_name,'] ')
        jjprint('-----------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
        jjprint(paste0(tok_descr, 'Ok. Modelo LDA entrenado [', Vector_num_features[tok_name], ' topics].'), logfile = s_Fic_log)
      }
      # save(ListaLDAfit, file = paste0(s_input_path, 'borrar4.RData'))
      # load(paste0(s_input_path, 'borrar4.RData')); gc();
    }
    # retVal <- ListaLDAfit
  } else {
    for(tok_name in c("frases", "palabs", "caract"))
    { # tok_name = "palabs"
      if(!is.null(ListaLDAfit[[tok_name]]))
      {
        tok_descr <- paste0(txtWhat, '[',tok_name,'] ')
        jjprint('-----------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
        myLDAfit <- ListaLDAfit[[tok_name]] # ListaLDAfit[["frases"]] / ListaLDAfit[["palabs"]]
        mis_features[[tok_name]] <- mis_features[[tok_name]][sapply(mis_features[[tok_name]], typeof) == "list"] # Nos quedamos con los elementos que son listas (i.e. quitamos el primero que es un número) [NOTA: Un data.frame es una lista!]
        jjprint(paste0(tok_descr, 'Creando DT final: (rbindlist de los ', length(mis_features[[tok_name]]), ' DTs de [mis_features])...'), logfile = s_Fic_log)
        mis_features_def <- rbindlist(mis_features[[tok_name]], use.names = TRUE)
        jjprint(paste0(tok_descr, 'Incluimos features en dataset_txt: (merge en los ', nrow(dataset_txt), ' docs)...'), logfile = s_Fic_log)
        setkeyv(mis_features_def, pk) # num_features columnas (probabilidad de cada topic) para cada texto_id
        # Solo falta asignarlas a nuestros ID iniciales (con merge pero devolviendo la merged DT)
        setkeyv(dataset_txt, pk) # num_features columnas (probabilidad de cada topic) para cada texto_id
        stopifnot(Vector_num_features[tok_name] == ncol(mis_features_def) - 1)
        dataset_txt <- merge(dataset_txt, mis_features_def, by = pk, all.x = TRUE)
        jjprint(paste0(tok_descr, 'Ok. Features creadas [', ncol(mis_features_def) - 1, ' topics].'), logfile = s_Fic_log)
      }
    }
    # retVal <- dataset_txt # ListaLDAfit / dataset_txt
  }
  jjprint('---------------------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
  jjprint(paste0(txtWhat, 'crearTextFeatures(b_entrenar:', b_entrenar, ') - Final'), logfile = s_Fic_log)
  jjprint('---------------------------------------------------------------------------------------', logfile = s_Fic_log, b_add_datetime = FALSE)
  if(b_entrenar)
  {
    retVal_train <- list( # para b_entrenar == TRUE
      ListaLDAfit = ListaLDAfit,
      dfStats = dfStats
    )
    retVal <- retVal_train
  } else {
    retVal_fit <- list( # para b_entrenar == FALSE
      dataset_txt = dataset_txt,
      dfStats = dfStats
    )
    retVal <- retVal_fit
  }
  return(retVal)
}
# ----------------------------------------------------------------
xgb_leerVarsModelo <- function(fichModelo, s_modelos_path)
{
  fich <- str_replace(str_replace(fichModelo, pattern = "\\.modelo$", replacement = ""), pattern = "\\_[0-9]*\\.[0-9][0-9][0-9]$", replacement = "\\.txt")
  stopifnot(file.exists(paste0(s_modelos_path, fich)))
  vars <- read.delim(file = paste0(s_modelos_path, fich), header = FALSE, col.names = c("varNum", "VarName", "featGain"), quote = "", dec = ".",
                     colClasses = c("character", "character", "numeric"), skip = 2)
  ultF <- nrow(vars)
  if(anyNA(vars$featGain))
  {
    ultF <- min(which(is.na(vars$featGain))) - 1
    if(vars$varNum[(ultF+1)] == "Feature")
      vars <- vars[1:ultF,]
  }
  stopifnot(all(vars$varNum[1:ultF] == c(1:ultF)))
  # Ya tenemos TODAS las vars del fichero ordenadas por Gain (desc.).
  return(vars)
}
xgb_prep_datos <- function(mi_dt, b_verbose = 2, maxImportanceNumVars = 0, fichero_Modelo = NULL, s_modelos_path = "", G_maxNumFeatures_OneHotEnc = 300, obj_var_name = 'clicked', s_Fic_log = NULL)
{
  # if(!is.null(mi_dt))
  # {
  #   # # mi_dt_nombre <- deparse(substitute(mi_dt)) # Para los eval(parse())
  #   # if("publish_time" %in% colnames(mi_dt))
  #   # {
  #   #   if(!("publish_timestamp" %in% colnames(mi_dt)))
  #   #   {
  #   #     setkey(mi_dt, publish_time)
  #   #     mi_dt[publish_time == "", publish_time := NA]
  #   #     setkey(mi_dt, publish_time)
  #   #     mi_dt[, publish_timestamp := 1000 * as.numeric(as.POSIXct(publish_time)) - as.numeric(1465876799998), by = publish_time] # (ms since 1970-01-01 - 1465876799998)
  #   #     setkey(mi_dt, NULL)
  #   # } }
  #   # mi_dt[, dia := as.integer(dia)]
  #   # # sapply(mi_dt, uniqueN)
  #   # # str(mi_dt)
  # } # if(!is.null(mi_dt))
  
  xgb_mi_version <- 0
  xgb_mi_version <- 1 # Añadidas vars 'Vari_medio_num' y 'Vari_b_...'
  xgb_mi_version <- 2 # Quitamos categ.NA para las categs (añadida como Level previamente en las variables con NAs)
  xgb_mi_version <- 3 # Añadimos numéricas LDA_frases54xxx (Latent Dirichlet Allocation -LDA- con los textos)
  xgb_mi_version <- 4 # Añadimos numéricas LDA_palabs54xxx (Latent Dirichlet Allocation -LDA- con los textos)
  # xgb_mi_version <- 5 # Numéricas LDA_frases27xxx y LDA_palabs27xxx (SIN LDA_*54)
  # xgb_mi_version <- 6 # Numéricas LDA_frases27xxx y LDA_palabs27xxx (CON LDA_*54)
  xgb_mi_version <- 7   # Numéricas LDA_frases72xxx y LDA_palabs72xxx (SIN LDA_*54, SIN LDA_*27)
  xgb_mi_version <- 8 # Numéricas LDA_frases72xxx y LDA_palabs72xxx (CON LDA_*54, SIN LDA_*27)
  xgb_mi_version <- 9 # Numéricas LDA_caract72xxx  (SIN LDA_*54, SIN LDA_*27)
  xgb_mi_version <- 10# Numéricas LDA_*72xxx (caract+palabs+frases) (SIN LDA_*54, SIN LDA_*27)
  xgb_mi_version <- 11# Numéricas LDA_caract72xxx  (CON todo!)
  # xgb_mi_version <- 12# Numéricas LDA_caract45xxx  (SIN LDA_*54, SIN LDA_*27, SIN LDA_*72)
  # xgb_mi_version <- 13# Numéricas LDA_*72xxx, caract45 (caract45+palabs72+frases72) (SIN LDA_*54, SIN LDA_*27)
  xgb_mi_version <- 14# Numéricas LDA_caract45xxx  (CON todo!)
  
  if(is.null(mi_dt))
    return(list(xgb_mi_version,NULL)) # Sólo devolvemos la versión
  
  # one-hot-encoding categorical features:
  if(exists("categs")) { rm(categs) }
  categs = c('Gene')
  categs = c(categs, 'Variation') # V.0
  categs = c(categs, "Vari_b_ampl", "Vari_b_del", "Vari_b_ins", "Vari_b_Mut") # V.1
  categs = c(categs, "Vari_b_Fus", "Vari_b_Trun", "Vari_first", "Vari_last") # V.1
  
  if(exists("numerics")) { rm(numerics) }
  numerics = c('Vari_medio_num') # V.1
  if(xgb_mi_version == 3)
  {
    # LDA_frases54*:
    LDA_vars <- colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases54_')]
  } else if(xgb_mi_version == 4)
  {
    # LDA_*54*:
    LDA_vars <- colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases54_')]
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs54_')])
  } else if(xgb_mi_version == 5)
  {
    # LDA_*27*:
    LDA_vars <- colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases27_')]
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs27_')])
  } else if(xgb_mi_version == 6)
  {
    # LDA_*27* & LDA_*54*:
    LDA_vars <- colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases27_')]
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs27_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases54_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs54_')])
  } else if(xgb_mi_version == 7)
  {
    # LDA_*72*:
    LDA_vars <- colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases72_')]
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs72_')])
  } else if(xgb_mi_version == 8)
  {
    # LDA_*72* & LDA_*54*:
    LDA_vars <- colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases72_')]
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs72_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases54_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs54_')])
  }  else if(xgb_mi_version == 9) # Añadimos "character" LDA (ademas de palabs y frases)
  {
    # LDA_caract*72*:
    LDA_vars <- colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_caract72_')]
  } else if(xgb_mi_version == 10)
  {
    # LDA*72* (caract, palabs, frases):
    LDA_vars <- colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_caract72_')]
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases72_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs72_')])
  }else if(xgb_mi_version == 11)
  {
    # LDA_caract*72* & LDA_*54* & LDA_*27*:
    LDA_vars <- colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_caract72_')]
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases72_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs72_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases54_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs54_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases27_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs27_')])
  }  else if(xgb_mi_version == 12) # (caract45)
  {
    # LDA_caract*45*:
    LDA_vars <- colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_caract45_')]
  }else if(xgb_mi_version == 13)
  {
    # LDA_caract*45* & LDA_*72*: # (caract45+palabs72+frases72) (SIN LDA_*54, SIN LDA_*27)
    LDA_vars <- colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_caract45_')]
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases72_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs72_')])
  }else if(xgb_mi_version == 14)
  {
    # LDA_caract*72*+*45* & LDA_*54* & LDA_*27*:
    LDA_vars <- colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_caract45_')]
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_caract72_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases72_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs72_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases54_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs54_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_frases27_')])
    LDA_vars <- c(LDA_vars, colnames(mi_dt)[startsWith(colnames(mi_dt), prefix = 'LDA_palabs27_')])
  } else { stopifnot(xgb_mi_version < 15) }
  
  numerics = c(numerics, LDA_vars)
  if(b_verbose >= 1)
    for(prfx in unique(str_extract(LDA_vars, "LDA\\_.*\\_")))
    {
      n <- length(LDA_vars[startsWith(LDA_vars, prfx)])
      jjprint(paste0('Añadidas ', n,' variables ', prfx ,'xx ...'), logfile = s_Fic_log, b_add_datetime = TRUE)
    }
  
  if(exists("categs"))
  {
    if(b_verbose >= 1)  jjprint('Procesando categs...', logfile = s_Fic_log, b_add_datetime = TRUE)
    mi_dt2 <- mi_dt[, c("ID", categs), with=FALSE] # Copiamos categóricas en otro data.table
    categsv8 <- vector(mode = "character")
    # One-Hot-Encoding de las categs con "pocos" levels:
    for(miCampo in categs) {
      if(!is.factor(mi_dt[,miCampo,with=FALSE][[1]])) {
        misLvls <- unlist(unique(mi_dt[,miCampo,with=FALSE]))
      } else {
        misLvls <- levels(mi_dt[,miCampo,with=FALSE][[1]])
      }
      numLvls <- length(misLvls)
      if(numLvls < G_maxNumFeatures_OneHotEnc){
        flds <- paste0(miCampo, 1:numLvls)
        for(lvl in 1:numLvls) {
          mi_dt2[, flds[lvl] := as.numeric(ifelse(get(miCampo) == misLvls[lvl], 1, 0))]
          # v.2 mi_dt2[is.na(get(miCampo)), flds[lvl] := 0] # Por si acaso, ya que vamos a añadir el campo xxNA
        }
        stopifnot(sum(mi_dt2[!is.na(get(miCampo)),flds,with=F]) == nrow(mi_dt2[!is.na(get(miCampo)),]))
        categsv8 <- c(categsv8, flds)
        # v.2 mi_dt2[, paste0(miCampo, 'NA') := as.numeric(ifelse(is.na(get(miCampo)), 1, 0))] # Otra más para los NAs
        # v.2 categsv8 <- c(categsv8, paste0(miCampo, 'NA'))
        # v.2 if(b_verbose >= 1)  jjprint(paste0('Variable categ. ', miCampo, ' Ok (', numLvls, ' niveles + NA)'), logfile = s_Fic_log, b_add_datetime = TRUE)
        if(b_verbose >= 1)  jjprint(paste0('Variable categ. ', miCampo, ' Ok (', numLvls, ' niveles)'), logfile = s_Fic_log, b_add_datetime = TRUE) # v.2
      } else {
        if(b_verbose >= 1)  jjprint(paste0('Variable categ. ', miCampo, ' tiene más de ', G_maxNumFeatures_OneHotEnc, ' niveles (', numLvls, '!)'), logfile = s_Fic_log, b_add_datetime = TRUE)
        jjprint(paste0('Variable categ. ', miCampo, ': PENDIENTE ver qué hacemos con ella...'), logfile = s_Fic_log, b_add_datetime = TRUE)
      }
    }
    # # Quitamos variables con un único valor (constantes):
    # varsUnValor <- names(sapply(mi_dt2[,categsv8,with=F], uniqueN)[sapply(mi_dt2[,categsv8,with=F], uniqueN) == 1])
    # if(length(varsUnValor) != 0)
    # {
    #   if(b_verbose >= 1)  jjprint(paste0('Warning: Quitamos variables categ. con un único valor: [', paste(varsUnValor, collapse = ','), ']'), logfile = s_Fic_log, b_add_datetime = TRUE)
    #   categsv8 <- categsv8[!(categsv8 %in% varsUnValor)]
    # }
    
    # if("platform" %in% categs)
    # {
    #   mi_dt2[, platform1 := ifelse(platform == 1, 1, 0)]
    #   mi_dt2[, platform2 := ifelse(platform == 2, 1, 0)]
    #   mi_dt2[, platform3 := ifelse(platform == 3, 1, 0)]
    #   categsv8 <- c(categsv8, "platform1", "platform2", "platform3")
    #   # mi_dt2[, platformNA := ifelse(is.na(platform), 1, 0)] # Otra más para los NAs
    #   # categsv8 <- c(categsv8, "platformNA")
    # }
    # # # sparse_matrix <- model.matrix(~ 0 + get(categs), data = campaign)
    # # # # # install.packages("caret")
    # # # # library(lattice)
    # # # # library(ggplot2)
    # # # # library(caret)
    # # formula <- paste0("~ ", paste(categs, collapse = "+"))
    # # dummies <- dummyVars(formula = formula, data = mi_dt)
  }
  
  importanceVars <- vector("character", 0)
  if(!is.null(fichero_Modelo))
  {
    if(b_verbose >= 1)  jjprint('Seleccionando variables del modelo...', logfile = s_Fic_log, b_add_datetime = TRUE)
    importanceVars <- xgb_leerVarsModelo(fichero_Modelo, s_modelos_path)$VarName
    if(b_verbose >= 1)  jjprint(paste0('Seleccionadas las ', length(importanceVars),' variables del modelo...'), logfile = s_Fic_log, b_add_datetime = TRUE)
  } else if(maxImportanceNumVars != 0)
  {
    if(b_verbose >= 1)  jjprint('Seleccionando variables de los siguientes ficheros:', logfile = s_Fic_log, b_add_datetime = TRUE)
    # Creamos la lista de las maxImportanceNumVars variables más importantes:
    fichs <- dir(path = s_modelos_path, pattern = paste0("XGB_.*v[0-9]", str_pad(xgb_mi_version, 2, "left","0"), "\\.txt$"))
    stopifnot(length(fichs) != 0)
    jjprint(fichs, logfile = s_Fic_log, b_add_datetime = FALSE)
    mis_vars <- data.table(mivar = character(0), mipos = integer(0))
    for(fich in fichs)
    {
      vars <- xgb_leerVarsModelo(fich, s_modelos_path)
      if(nrow(vars) < maxImportanceNumVars)
        next # Este fichero tiene menos variables (o no tiene ninguna), así que no lo usamos...
      # Nos quedaremos, después, con las maxImportanceNumVars más frecuentes en todos los ficheros que encontremos:
      mis_vars <- rbindlist(list(mis_vars, list(mivar = vars$VarName, mipos = as.integer(vars$varNum))), use.names = T, fill = F)
    }
    setkey(mis_vars, mivar)
    pp <- mis_vars[, .(mx=max(mipos),mn=min(mipos),av=mean(mipos),sd=sd(mipos),cnt=.N), by=mivar]
    # plot(x = 1:nrow(pp), y = pp$mx)
    # plot(x = 1:nrow(pp), y = pp$av)
    # plot(x = 1:nrow(pp), y = pp$mn)
    setorderv(pp, "mx") # ordenamos por máximo (i.e. la "peor" posición de cada variable)
    # Y nos quedamos con las mejores:
    importanceVars <- pp$mivar[1:maxImportanceNumVars]
    if(b_verbose >= 1)  jjprint(paste0('Nos quedamos con las ', maxImportanceNumVars,' variables más importantes (XGB Feature Gain)...'), logfile = s_Fic_log, b_add_datetime = TRUE)
  } # if(maxImportanceNumVars != 0)
  
  # Filtramos variables por importancia:
  if(length(importanceVars) != 0)
  {
    if(exists("numerics"))
    {
      numerics <- numerics[numerics %in% importanceVars]
      if(length(numerics) == 0)
      {
        if(b_verbose >= 1)  jjprint("NOTA: Quitamos numéricas porque no ha 'sobrevivido' ninguna", logfile = s_Fic_log, b_add_datetime = TRUE)
        rm(numerics) # Quitamos numéricas porque no ha 'sobrevivido' ninguna...
      }
    }
    if(exists("categs"))
    {
      categsv8 <- categsv8[categsv8 %in% importanceVars]
      if(length(categsv8) == 0)
      {
        if(b_verbose >= 1)  jjprint("NOTA: Quitamos categóricas porque no ha 'sobrevivido' ninguna", logfile = s_Fic_log, b_add_datetime = TRUE)
        rm(categs, categsv8) # Quitamos categóricas porque no ha 'sobrevivido' ninguna...
      }
    }
  } # if(length(importancevars) != 0)
  
  # # V.10 (y 3) - Y finalmente. añadimos el contenido de las variables numerics_clk_1 & categs_clk_1 a mi_dt:
  # if(length(numerics_clk_1) != 0)
  #   numerics_clk_1 <- numerics_clk_1[paste0('clk1_', numerics_clk_1) %in% numerics]
  # if(length(numerics_clk_1) != 0)
  # {
  #   if(b_verbose >= 1)  jjprint(paste0('Creando variables numerics_clk_1 (', length(numerics_clk_1), ' vars.)...'), logfile = s_Fic_log, b_add_datetime = TRUE)
  #   miconv_ini <- "mi_dt[,c("
  #   miconv_fin <- "), by='display_id']"
  #   for(i in seq.int(1,length(numerics_clk_1), 25))
  #   {
  #     j <- min(i+24, length(numerics_clk_1))
  #     if(b_verbose >= 1)  jjprint(paste0('Creando variables numerics_clk_1[', i, ':', j, ']...'), logfile = s_Fic_log, b_add_datetime = TRUE)
  #     miconv_mid_1 <- paste0("'clk1_", numerics_clk_1[i:j], "'", collapse = ",")
  #     miconv_mid_2 <- paste0(numerics_clk_1[i:j], "[1]", collapse = ",")
  #     if(b_verbose >= 2)  jjprint(miconv_mid_1, logfile = s_Fic_log, b_add_datetime = FALSE)
  #     miconv <- paste0(miconv_ini, miconv_mid_1, ") := list(", miconv_mid_2, miconv_fin)
  #     eval(parse(text = miconv))
  #   }
  #   # for(j in numerics_clk_1)
  #   # { mi_dt[, c(paste0('clk1_', j)) := get(j)[1], by = "display_id"]
  #   #   if(b_verbose >= 2)  jjprint(paste0('Creando variable ', j, ' (', which(numerics_clk_1 == j),'/',length(numerics_clk_1),')...'), logfile = s_Fic_log, b_add_datetime = TRUE)
  #   # }
  # }
  # 
  # if(exists("categs"))
  # {
  #   categsv8 = c(categsv8, 'numAds') # V.11 (numAds como última columna)
  # } else {
  #   numerics = c(numerics, 'numAds') # V.11 (numAds como última columna)
  # }
  
  if(obj_var_name %in% colnames(mi_dt))
  {
    obj_var.col <- as.numeric(mi_dt[,(obj_var_name),with=FALSE][[1]])
  } else 
  {
    obj_var.col <- as.numeric(rep.int(0, nrow(mi_dt)))
  }
  if(exists("categs") && exists("numerics"))
  {
    if(b_verbose >= 1)  jjprint('Preparando data.table [numerics + categs]...', logfile = s_Fic_log, b_add_datetime = TRUE)
    return(list(xgb_mi_version,
                cbind(setNames(data.table(obj_var.col), obj_var_name),
                      mi_dt[, numerics, with=FALSE],
                      mi_dt2[, categsv8, with=FALSE] # predict(dummies, newdata = mi_dt)
                )
    )    )
  } else if(exists("categs"))
  {
    if(b_verbose >= 1)  jjprint('Preparando data.table [Sólo categs]...', logfile = s_Fic_log, b_add_datetime = TRUE)
    return(list(xgb_mi_version,
                cbind(setNames(data.table(obj_var.col), obj_var_name),
                      mi_dt2[, categsv8, with=FALSE] # predict(dummies, newdata = mi_dt)
                )
    )    )
  } else if(exists("numerics"))
  {
    if(b_verbose >= 1)  jjprint('Preparando data.table [Sólo numerics]...', logfile = s_Fic_log, b_add_datetime = TRUE)
    return(list(xgb_mi_version,
                cbind(setNames(data.table(obj_var.col), obj_var_name),
                      mi_dt[, numerics, with=FALSE]
                )
    )    )
  } else
  {
    if(b_verbose >= 1)  jjprint('Preparando data.table [ERROR: Ni numerics Ni categs]...', logfile = s_Fic_log, b_add_datetime = TRUE)
    return(list(xgb_mi_version,
                cbind(setNames(data.table(obj_var.col), obj_var_name),
                      data.table(Error="No hay variables numerics ni categs")
                )
    )    )
  }
}
# ----------------------------------------------------------------
lgb_prep_datos <- function(mi_dt, b_verbose = 2, maxImportanceNumVars = 0, fichero_Modelo = NULL, s_modelos_path = "", G_maxNumFeatures_OneHotEnc = 1000, obj_var_name = 'Class', s_Fic_log = NULL)
{
  return(xgb_prep_datos(mi_dt = mi_dt, b_verbose = b_verbose, maxImportanceNumVars = maxImportanceNumVars, fichero_Modelo = fichero_Modelo, s_modelos_path = s_modelos_path, G_maxNumFeatures_OneHotEnc = G_maxNumFeatures_OneHotEnc, obj_var_name = obj_var_name, s_Fic_log = s_Fic_log))
}
# ----------------------------------------------------------------
nbay_prep_datos <- function(mi_dt, b_verbose = 2)
{
  # mi_dt_nombre <- deparse(substitute(mi_dt)) # Para los eval(parse())
  setkey(mi_dt, NULL)
  # mi_dt[, hora := as.integer(hora)]
  # mi_dt[, dia := as.integer(dia)]
  # sapply(mi_dt, uniqueN)
  # str(mi_dt)
  
  # Para Naive Bayes ponemos todo como categóricas:
  categs = c('platform', 'pais')
  categs = c(categs, 'uuid', 'geo_location', 'geo_loc.country')
  categs = c(categs, 'ad_campaign_id', 'ad_advertiser_id')
  categs = c(categs, 'source_id', 'publisher_id')
  categs = c(categs, 'ad_source_id', 'ad_publisher_id')
  categs = c(categs, 'publish_timestamp', 'ad_publish_timestamp')
  
  numerics = c('dia', 'hora', 'timestamp')
  # numerics = c(numerics, 'publish_timestamp', 'ad_publish_timestamp') # NAs ???
  numerics = c(numerics, 'topics_prob', 'ad_topics_prob')
  numerics = c(numerics, 'topic_prob_1', 'topic_prob_2', 'topic_prob_3', 'topic_prob_4', 'topic_prob_5', 'topic_prob_6', 'topic_prob_7', 'topic_prob_8')
  numerics = c(numerics, 'entities_prob', 'ad_entities_prob')
  numerics = c(numerics, 'entity_prob_1', 'entity_prob_2', 'entity_prob_3', 'entity_prob_4', 'entity_prob_5', 'entity_prob_6', 'entity_prob_7', 'entity_prob_8')
  numerics = c(numerics, 'categories_prob', 'ad_categories_prob')
  numerics = c(numerics, 'category_prob_1', 'category_prob_2')
  
  nbay_mi_version <- 5 # Numéricas normalizadas
  nbay_mi_version <- 6 # publish_timestamp y ad_publish_timestamp como categóricas
  
  if(b_verbose >= 1)  print('Convertimos numéricas a numeric (double) y normalizamos: (WHY?)')
  for(miconv in paste0("mi_dt[, ", numerics, " := scale(as.numeric(", numerics, "), center=T, scale=T)]"))
  {
    if(b_verbose == 2)  print(miconv)
    eval(parse(text = miconv))
  }
  if(b_verbose >= 1)  print('Convertimos categóricas a Factor:')
  for(miconv in paste0("mi_dt[, ", categs, " := as.factor(", categs, ")]"))
  {
    if(b_verbose == 2)  print(miconv)
    eval(parse(text = miconv)) # Convertimos a Factor
  }
  for(miconv in paste0("if(anyNA(mi_dt$", categs, "))  mi_dt$", categs, " <- addNA(mi_dt$", categs, ")"))
  {
    if(b_verbose == 2)  print(miconv)
    eval(parse(text = miconv)) # Insertamos xxx_NA, un level para los NA, si hay
  }
  if(b_verbose >= 1)  print('Preparando data.table...')
  if("clicked" %in% colnames(mi_dt))
  {
    clicked.col <- as.integer(mi_dt$clicked)
  } else 
  {
    clicked.col <- as.integer(rep.int(0, nrow(mi_dt)))
  }
  clicked.col <- as.factor(clicked.col)
  
  return(list(nbay_mi_version,
              cbind(clicked = clicked.col,
                    mi_dt[, numerics, with=FALSE],
                    mi_dt[, categs, with=FALSE])
  )    )
}
# ----------------------------------------------------------------
randfor_prep_datos <- function(mi_dt, b_verbose = 2)
{
  return(nbay_prep_datos(mi_dt = mi_dt, b_verbose = b_verbose))
}
# ----------------------------------------------------------------
predict_testset <- function(nombres_modelos
                            , filename, s_input_path, s_output_path, s_modelos_path
                            , i_sDescr = "" # "XGB Predicting (extreme gradient boosting)"
                            , obj_var_name = 'clicked', num_classes = 2, G_maxNumFeatures_OneHotEnc = 300
                            , FUN_prep_datos, prep_datos_b_verbose = 0
                            , FUN_X
                            , FUN_Predict
                            , FUN_loadmodelo, nombreModeloLoaded
                            # , i_sDescr = "LGB (LightGBM) Predicting" # LightGBM
                            # , FUN_prep_datos = lgb_prep_datos # LightGBM
                            # , FUN_X = data.matrix # FUN_X = function(x){ return(data.matrix(x)) } # LightGBM
                            # , FUN_Predict = function(modelo, X){ return(predict(modelo, X)) } # LightGBM
                            # , FUN_loadmodelo = lgb.load, nombreModeloLoaded = "" # LightGBM
                            # , i_sDescr = "NBAY Predicting (Naive Bayes)" # NaiveBayes
                            # , FUN_prep_datos = nbay_prep_datos # NaiveBayes
                            # , FUN_X = identity # FUN_X = function(x){ return(data.matrix(x)) } # NaiveBayes
                            # , FUN_Predict = function(modelo, X){ return(predict(modelo, newdata = X, type = 'prob')$posterior[,"1"]) } # NaiveBayes
                            # , FUN_loadmodelo = load, nombreModeloLoaded = "nbay" # NaiveBayes
                            # , i_sDescr = "XGB Predicting (extreme gradient boosting)" # XGBoost
                            # , FUN_prep_datos = xgb_prep_datos # XGBoost
                            # , FUN_X = data.matrix # FUN_X = function(x){ return(data.matrix(x)) } # XGBoost
                            # , FUN_Predict = function(modelo, X){ return(xgboost::predict(modelo, X, missing = NA)) } # XGBoost
                            # , FUN_loadmodelo = xgb.load, nombreModeloLoaded = "" # XGBoost
                            # , FUN_prep_datos = randfor_prep_datos # RandomForest
                            # , FUN_X = data.matrix # FUN_X = function(x){ return(data.matrix(x)) } # RandomForest
                            # , FUN_Predict = function(modelo, X){ return(predict(modelo, newdata = X, type = 'prob')[,"1"]) } # RandomForest
                            # , FUN_loadmodelo = load, nombreModeloLoaded = "randfor" # RandomForest
                            , NUM_MODELOS = NUM_BLOQUES
                            , b_ReducirFicheros = FALSE # - COMENTADO - De momento no hace falta, pero podría hacer falta...
                            , modelos_weights = NULL # Para mezclar modelos sin ton ni son (bueno, con estos pesos)...
)
{
  # Leemos modelo(s):
  if(length(nombres_modelos[nombres_modelos != ""]) == 0)
  {
    tmp_nombres_modelos <- str_replace(dir(path = s_modelos_path, pattern = paste0(str_replace(filename, "_[0-9]*\\.[0-9][0-9][0-9]$", ""), '.*.modelo')), pattern = "\\.modelo", replacement = "")
    # Ordenamos modelos por numModelo (aquí NO da siempre igual, porque ya NO los vamos a promediar siempre)
    nombres_modelos[as.integer(substr(tmp_nombres_modelos, nchar(tmp_nombres_modelos)-2, nchar(tmp_nombres_modelos)))] <- tmp_nombres_modelos
  }
  stopifnot(length(nombres_modelos[nombres_modelos != ""]) == NUM_MODELOS)
  b_con_media_ponderada <- !is.null(modelos_weights)
  if(b_con_media_ponderada)
  {
    stopifnot(length(modelos_weights[!is.na(modelos_weights)]) == NUM_MODELOS)
    if(substr(filename, 1, 2) != "w_")
      filename <- paste0("w_", filename)
  }
  # s_Fichero_submit <- paste0(unique(str_replace(string = nombres_modelos, pattern = "_[0-9]*\\.[0-9][0-9][0-9]$", replace = "")), '_submit.csv')
  s_Fichero_submit <- paste0(str_replace(filename, "_[0-9]*\\.[0-9][0-9][0-9]$", ""), '_submit.csv')
  s_Fic_submit_log <- paste0(str_replace(filename, "_[0-9]*\\.[0-9][0-9][0-9]$", ""), '_submit.log')
  s_Fic_submit_log <- paste0(s_output_path, s_Fic_submit_log)
  if(file.exists(s_Fic_submit_log))
    jjprint('---------------------------------------------------------------------------------------', logfile = s_Fic_submit_log, b_add_datetime = FALSE) # Línea separadora para el inicio de este log...
  jjprint(paste0('Preparando submit [', s_Fichero_submit, ']...'), logfile = s_Fic_submit_log)

  if(b_ReducirFicheros)
  {
    # systime_ini2 <- proc.time()
    # # Creamos una división en cada bloque de otros 16 bloques (para acelerar y evitar problemas de memoria):
    # if(!file.exists(file.path(s_input_path, get_batch_test_filename(NUM_BLOQUES, NUM_BLOQUES))))
    # {
    #   for(numBatch in 1:NUM_BLOQUES)
    #   {
    #     jjprint(paste0('Creando sub_bloques de testset...'), logfile = s_Fic_submit_log)
    #     testset_big <- leer_batch_test(numBatch, i_sDescr, s_input_path)
    #     setkey(testset_big, display_id)
    #     disps <- unique(testset_big[,.(display_id)], by = "display_id") # data.table sin duplicados por clave (y de una columna)
    #     index <- 1:nrow(disps)
    #     # NOTA: testset_big se deja con el orden original para predecir y crear el submitset en orden.
    #     rows_per_block <- 1 + as.integer(nrow(disps) / NUM_BLOQUES)
    #     for(numSubBatch in 1:NUM_BLOQUES)
    #     {
    #       s_fich <- get_batch_test_filename(numBatch, numSubBatch)
    #       index_from <- 1 + (numSubBatch-1) * rows_per_block
    #       index_to <- numSubBatch * rows_per_block
    #       if(index_to > nrow(disps)) index_to <- nrow(disps)
    #       disps_set <- disps[index_from:index_to,]
    #       testset <- merge(testset_big, disps_set, by = 'display_id')
    #       save(testset, file = file.path(s_input_path, s_fich)) # testset_XXX_001 a testset_XXX_016
    #       rm(testset); gc()
    #       jjprint(paste0(s_fich, " Ok."), logfile = s_Fic_submit_log)
    #       # print(paste0(as.double((proc.time() - systime_ini2)['elapsed'])/60, ' minutos.'))
    #       minutos_pend <- (as.double((proc.time() - systime_ini2)['elapsed'])/60) * ( (NUM_BLOQUES^2) / ((numBatch-1) * NUM_BLOQUES + numSubBatch)  -  1)
    #       if(minutos_pend < 60) print(paste0('Faltan aprox. ',minutos_pend, ' minutos.')) else print(paste0('Faltan aprox. ',minutos_pend/60, ' horas.'))
    #     }
    #     rm(testset_big, disps, disps_set); gc()
    #   }
    # }
  }
  
  nombres_modelos <- str_replace(nombres_modelos, pattern = "\\.modelo$", replacement = "")
  modelos <- list() # Inicializamos lista de modelos y cargamos los modelos:
  for(i in 1:NUM_MODELOS)
  {
    str_tmp <- paste0('Leyendo ', str_pad(i, 3, "left" ,"0"), ' (modelo ', i, ' de ', NUM_MODELOS, '): [', nombres_modelos[i],']...')
    ix_mod_ya_cargado <- ifelse(i == 1, integer(0), match(nombres_modelos[i], nombres_modelos[1:(i-1)]))
    if(!is.na(ix_mod_ya_cargado))
    { # El modelo ya está cargado:
      str_tmp <- paste0(str_tmp, ' Ok. (Ya estaba cargado)')
      jjprint(str_tmp, logfile = s_Fic_submit_log)
      ix_mod_ya_cargado <- min(ix_mod_ya_cargado)
      modelos[[i]] <- ix_mod_ya_cargado # Luego hay que verificar si el modelo es o no es un integer
    } else {
      # Leemos fichero del modelo:
      jjprint(str_tmp, logfile = s_Fic_submit_log)
      stopifnot(file.exists(paste0(s_modelos_path, nombres_modelos[i] ,'.modelo')))
      modelos[[i]] <- FUN_loadmodelo(paste0(s_modelos_path, nombres_modelos[i] ,'.modelo'))
      if(nombreModeloLoaded != "")
      {
        modelos[[i]] <- get(nombreModeloLoaded); rm(list = c(nombreModeloLoaded), inherits = TRUE); gc()
      } # deparse(nombreModeloLoaded) deparse(get(nombreModeloLoaded))
    }
    jjprint(str_tmp, logfile = s_Fic_submit_log)
  }
  
  n_versiones <- as.integer(substring(stringr::str_extract(nombres_modelos, pattern = "v[0-9]+"), first = 2))
  n_version <- n_versiones[1]
  jjprint(paste0('Versión ', n_version, '. Con_media_ponderada = ', b_con_media_ponderada), logfile = s_Fic_submit_log)
  if(!b_con_media_ponderada)
    stopifnot(all(n_versiones == n_version)) # La versión debe coincidir si NO estamos haciendo medias ponderadas
  systime_ini2 <- proc.time()
  for(numBatch in 1:NUM_BLOQUES)
  {
    for(numSubBatch in 1:NUM_BLOQUES)
    {
      if(!file.exists(paste0(s_input_path, get_batch_test_filename(NUM_BLOQUES, NUM_BLOQUES))))
      {
        # No hay subdivisión de testsets. Lo hacemos con los "grandes":
        if(numSubBatch != 1)  break # Sólo el primero!
        testset <- leer_batch_test(numBatch, i_sDescr, s_input_path)
      } else {
        testset <- leer_batch_test(numBatch, numSubBatch = numSubBatch, s_descr = i_sDescr, s_input_path = s_input_path)
      }
      gc()
      dt_all <- FUN_prep_datos(mi_dt = testset, b_verbose = prep_datos_b_verbose, s_modelos_path = s_modelos_path, G_maxNumFeatures_OneHotEnc = G_maxNumFeatures_OneHotEnc, obj_var_name = obj_var_name)[[2]]
      # print(colnames(dt_all))
      jjprint(paste0('Preparando matrix y predicting...', ' Versión ', n_version, '. Con_media_ponderada = ', b_con_media_ponderada), logfile = s_Fic_submit_log)
      X <- FUN_X(dt_all[,-1,with=FALSE])
      # y <- data.matrix(dt_all[,1,with=FALSE])
      y_pred <- vector(mode = "numeric", length = NUM_MODELOS)
      for(i in 1:NUM_MODELOS)
      {
        jjprint(paste0('Predicting ', str_pad(numBatch, 3, "left" ,"0"), '_', str_pad(numSubBatch, 3, "left" ,"0"), ' (modelo ', i, ' de ', NUM_MODELOS, '): ', nombres_modelos[i],'...', ' Versión ', n_version, '. Con_media_ponderada = ', b_con_media_ponderada), logfile = s_Fic_submit_log)
        # modelos[[i]] <- xgb.load(paste0(s_modelos_path, nombres_modelos[i] ,'.modelo'))
        
        if(n_version > 300) # Hay selección de variables (que puede ser diferente para cada modelo)
        {
          X <- FUN_X(FUN_prep_datos(mi_dt = testset, b_verbose = prep_datos_b_verbose, fichero_Modelo = nombres_modelos[i], s_modelos_path = s_modelos_path, G_maxNumFeatures_OneHotEnc = G_maxNumFeatures_OneHotEnc, obj_var_name = obj_var_name)[[2]][,-1, with=FALSE])
        }
        
        if(b_con_media_ponderada)
        {
          y_pred[i] <- list(modelos_weights[i] * 
                             FUN_Predict( if(is.integer(modelos[[i]])) modelos[[(modelos[[i]])]] else modelos[[i]], X) )
        } else {
          y_pred[i] <- list( FUN_Predict( if(is.integer(modelos[[i]])) modelos[[(modelos[[i]])]] else modelos[[i]], X) )
        }
      }
      if(NUM_MODELOS > 1)
      {
        if(num_classes == 2) {
          mean_y_pred <- rowSums(rbindlist(list(y_pred[1])))
          if(b_con_media_ponderada) {
            mean_y_pred <-  mean_y_pred / sum(modelos_weights)
          } else {
            mean_y_pred <- mean_y_pred / NUM_MODELOS
          }
        } else {
          for(j in 1:num_classes)
          {
            mean_y_pred[j] <- rowSums(rbindlist(list(y_pred[,j])))
            if(b_con_media_ponderada) {
              mean_y_pred[j] <-  mean_y_pred[j] / sum(modelos_weights)
            } else {
              mean_y_pred[j] <- mean_y_pred[j] / NUM_MODELOS
            }
          }
        }
      } else {
        mean_y_pred <- y_pred[[1]]
      }
      # Finalmente, ponemos las probs en el testset:
      if(num_classes == 2)   testset[, prob := mean_y_pred]
      if(num_classes != 2)   for(i in 1:num_classes) testset[, paste0('prob_class', i) := mean_y_pred[,i]]
      # Verificamos que no hay NAs (en las probs):
      if(num_classes == 2) stopifnot(!anyNA(testset$prob))
      if(num_classes != 2) stopifnot(!anyNA(testset[,paste0('prob_class', 1:num_classes),with=F]))
      # Guardamos...
      if(!file.exists(paste0(s_input_path, get_batch_test_filename(NUM_BLOQUES, NUM_BLOQUES))))
      {
        # No hay subdivisión de testsets. Lo hacemos con los "grandes":
        jjprint(paste0('Guardando ', s_Fichero_submit, ' (', (numBatch) , '/', NUM_BLOQUES, ')...'), logfile = s_Fic_submit_log)
      } else {
        jjprint(paste0('Guardando ', s_Fichero_submit, ' (', (numBatch-1) * NUM_BLOQUES + numSubBatch, '/', NUM_BLOQUES^2, ')...'), logfile = s_Fic_submit_log)
      }
      guardar_submit(testset = testset, fichero = s_Fichero_submit, s_output_path = s_output_path,
                     b_write_file_append = (numBatch>1 | numSubBatch>1)
                    )
      if(numSubBatch != 1)
      {
        minutos <- as.double((proc.time() - systime_ini2)['elapsed'])/60
        if(minutos < 60) str_tmp <- paste0(minutos, ' minutos.')  else  str_tmp <- paste0(minutos/60, ' horas.')
        jjprint(str_tmp, logfile = s_Fic_submit_log)
        minutos_pend <- (as.double((proc.time() - systime_ini2)['elapsed'])/60) * ( (NUM_BLOQUES^2) / ((numBatch-1) * NUM_BLOQUES + numSubBatch)  -  1)
        if(minutos_pend < 60) print(paste0('Faltan aprox. ', minutos_pend, ' minutos.')) else print(paste0('Faltan aprox. ', minutos_pend/60, ' horas.'))
      }
      if(exists("mean_y_pred")) rm(mean_y_pred)
      rm(testset, dt_all, X, y_pred); gc()
    }
    minutos <- as.double((proc.time() - systime_ini2)['elapsed'])/60
    if(minutos < 60) str_tmp <- paste0('(', (numBatch) , '/', NUM_BLOQUES, ')', minutos, ' minutos.')  else  str_tmp <- paste0(Sys.time(), ' - (', (numBatch) , '/', NUM_BLOQUES, ')', minutos/60, ' horas.')
    jjprint(str_tmp, logfile = s_Fic_submit_log)
    minutos_pend <- (as.double((proc.time() - systime_ini2)['elapsed'])/60) * ( NUM_BLOQUES / numBatch  -  1)
    if(minutos_pend < 60) print(paste0('Faltan aprox. ', minutos_pend, ' minutos.')) else print(paste0('Faltan aprox. ', minutos_pend/60, ' horas.'))
  }
}
# ----------------------------------------------------------------
mi_load_uniq_txt <- function(tipo = 'train', s_Fic_log = 'crearTextFeatures_mi_load_uniq_txt.log') # "tolower_uniq_train_txt.RData" / "tolower_uniq_test_txt.RData"
{
  stopifnot(tipo %in% c('train', 'test'))
  fich <- paste0('tolower_uniq_', tipo, '_txt.RData')
  if(!mi_load(s_input_path, fich))
  {
    if(tipo == 'train') {
      # Leemos full_trainset_txt:
      stopifnot(mi_load(s_input_path, "training_text.RData"))
      dataset_txt <- full_trainset_txt; rm(full_trainset_txt); gc()
      pk <- "train_texto_id"
    } else {
      # Leemos testset_txt:
      stopifnot(mi_load(s_input_path, "test_text.RData"))
      dataset_txt <- testset_txt; rm(testset_txt); gc()
      pk <- "test_texto_id"
    }
    jjprint(paste0('UNIQUE de textos...'), logfile = s_Fic_log)
    dataset_txt <- unique(dataset_txt[, c(pk, "texto"), with=F])
    mi_tiempo <- system.time({
      jjprint(paste0('char_tolower() ', nrow(dataset_txt), ' regs...'), logfile = s_Fic_log)
      dataset_txt[, texto := char_tolower(texto, keep_acronyms=TRUE)]
    })
    mi_tiempo_tolower <- paste0('Tiempo tolower(): ', jjfmt(mi_tiempo['elapsed']/60, 3), ' minutos (', nrow(dataset_txt), ' regs)')
    jjprint(paste0(mi_tiempo_tolower), logfile = s_Fic_log)
    setkeyv(dataset_txt, pk)
    mi_tiempo <- system.time({
      save(dataset_txt, pk, file = file.path(s_input_path, fich))
    })
    mi_tiempo_tolower <- paste0('Tiempo tolower(): ', jjfmt(mi_tiempo['elapsed']/60, 3), ' minutos (', nrow(dataset_txt), ' regs)')
    jjprint(paste0('Ok. [', fich, '] guardado.', mi_tiempo_tolower), logfile = s_Fic_log)
  }
  return(TRUE)
} # "tolower_uniq_train_txt.RData" / "tolower_uniq_test_txt.RData"
# ----------------------------------------------------------------
