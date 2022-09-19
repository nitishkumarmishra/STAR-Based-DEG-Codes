








##########  packages

require(Biostrings)
require(plyr)
require(dplyr)
require(stringr)
require(limma)
#require(cqn)
require(edgeR)
require(methods)
require(utils)
require(stats)
require(reshape2)
require(ggplot2)
require(Hmisc)
require(Matrix)  # rank
require(gplots)
require(cowplot)
require(scales) # muted
#require(GGally) # scatmat
require(data.table)

require(AnnotationDbi)




myParams.limmaFuncs.grey  <-  "grey54"

myParams.limmaFuncs.de.colors  <-  c(down = "firebrick3" , up = "chartreuse4" ,
				sig = "dodgerblue3" , notSig = "grey60" ,
				downUp =  "firebrick3" ,  upDown = "chartreuse4" ,
				downDown =  "firebrick3" ,  upUp = "chartreuse4" ,
				XdownLdownMore = "firebrick3" , XupLupMore = "chartreuse4" ,
				both = "dodgerblue3" , opposite = "brown" ,
				"+" = "orange" , "-" = "darkorchid2" )
myParams.limmaFuncs.de.colors.na  <- "grey90"

myParams.limmaFuncs.gg.scale.color  <-  scale_colour_manual( values = myParams.limmaFuncs.de.colors , na.value = myParams.limmaFuncs.de.colors.na , drop=FALSE )
myParams.limmaFuncs.gg.scale.fill  <-  scale_fill_manual( values = myParams.limmaFuncs.de.colors , na.value = myParams.limmaFuncs.de.colors.na , drop=FALSE )




























######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
##################   Differential Expression  ########################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################



get_project_info  <-  function( project = NULL , individual = NULL , runid = NULL , strategy = NULL , verbose = FALSE ) {



        dbdir  <-  "./data/my_databases"
        poss.names  <-  list.files( path = dbdir , pattern = "*.txt" )
        poss.names  <-  poss.names[!(poss.names %in% c("sample.db.txt","sample.db.notes.txt","public.db.txt","bad.1000G.txt",
                                                                "general.database.txt","organism.databases.txt","bad.hipsci.txt" , "inhouse.txt","ecoliOperonExpr.txt" ,
                                                                "protoVariant.database.txt" , "OLD.PRJEB9586.txt"))]


        db  <-  paste( dbdir , poss.names , sep = "/" )

        if (verbose) { show(db) }

        indf  <-  data.frame()
        for (dbx  in  db) {
                if (verbose) { verb("\t\tdbx = [%s]\n", dbx ) }

                invardb  <-  read.table( file = dbx , header = TRUE , sep = "\t" ,  row.names = NULL , stringsAsFactors = FALSE , quote = "" , fill = TRUE , comment.char = "" )

                invardb$project  <-  as.character(invardb$project)
                invardb$individual  <-  as.character(invardb$individual)
                invardb$runid  <-  as.character(invardb$runid)

                ### find
                sdf  <-  data.frame()
                if (!is.null(runid)) {
                        sdf  <-  invardb[ invardb$runid %in% runid ,,drop=FALSE]
                }
                if (!is.null(individual)) {
                        sdf  <-  invardb[ invardb$individual %in% individual ,,drop=FALSE]
                }
                if (!is.null(project)) {
                        sdf  <-  invardb[ invardb$project %in% project ,,drop=FALSE]
                } # project

                if (verbose) { verb("\t\t\t%d rows\n", nrow(sdf) ) }

                if (nrow(sdf) > 0) {
                        indf  <-  rbind.fill( indf , sdf )
                } # sdf
        } # dbx



        if (!is.null(strategy)) {
                stopifnot("strategy"  %in%  colnames(indf))
                indf  <-  indf[ indf$strategy %in% strategy  &  !is.na(indf$strategy) ,, drop=FALSE]
        } # strategy



        ### check
        if (nrow(indf) == 0) {
                verb("\n\n\nERROR!  unable to find entries!!!\n" )
                show(project)
                show(runid)
                stop()
        } # check

        return(indf)

} # get_project_info






process_project_input  <-  function( project = NULL , runid = NULL , individual = NULL , strategy = NULL , verbose = FALSE ) {

        if ( (is.character(project)  ||  is.character(runid))  &&  !is.data.frame(project)  && !is.data.table(project)  ) {
                project  <-  get_project_info( project = project , individual = individual , runid = runid , strategy = strategy , verbose = verbose )
        }

        stopifnot(!any(is.na(project$condition)))
        stopifnot(!any(is.na(project$runid)))
        stopifnot(!any(is.na(project$individual)))
        stopifnot(!any(is.na(project$strategy)))
        stopifnot(!any(is.na(project$layout)))
        stopifnot(!any(is.na(project$species)))

        if ("order" %in% colnames(project)) {
                stopifnot(!any(is.na(project$order)))

                if (is.data.table(project)) {
                        setkey(project , order , runid)
                } else {
                        project  <-  project[ order(project$order , project$runid) , ,drop=FALSE]
                }

                ### each condition must have a distinct order number, with a single order number per condition.
                if ("condition" %in% colnames(project)) {
                        uco  <-  unique(project[,c("condition","order"),drop=FALSE])
                        #stopifnot(!any(duplicated(uco$condition)))
                        stopifnot(!any(duplicated(uco$order)))
                } # condition
        } # order

        rownames(project)  <-  project$runid

        return(project)

} # process_project_input




parse_rlist  <-  function(      op = NULL ,
                                rlist = NULL ,
                                feats = NULL )  {

        general.cols <-  c(
                "rundate",
                "sname",
                "condition",
                "protocol",
                "drug",
                "stranded",
                "biorep",
                "mrna.org",
                "trna.org",
                "rrna.org",
                "read.length.low",
                "read.length.high",
                "strategy",
                "adaptor",
                "layout",
                "trim",
                "abc.dir",
                "order",
                "rrna.fa",
                "frompublic",
                "notes",
                "subunit")



        if (op == "load") {

                in.data  <-  read.table( file = rlist , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = ""  )
                colnames(in.data)  <-  general.cols[1:ncol(in.data)]
                in.data$sname  <-  as.character(in.data$sname)
                in.data$rundate  <-  as.character(in.data$rundate)
                return(in.data)

        } else if (op == "sort") {
                return( rlist[ order(rlist$order , rlist$protocol , rlist$rundate , rlist$sname) ,] )

        }




} # parse_rlist







get_DNA_level_for_strategy  <-  function( strategy )  {

        vec  <-  ifelse( strategy  %in% c("DNAseq","ChIPseq") , "DNA" , "RNA" )
        return(vec)

} # get_DNA_level_for_strategy





de_analysis_type_to_directory_and_tag  <-  function( type , project ) { 

        if (type == "mrna") {
                bdir  <-  "out/limma-voom.mrna"
                btag  <-  "limma-voom.mrna"
        } else if (type == "rrna-var") {
                bdir  <-  "out/limma-voom.rrna-var"
                btag  <-  "limma-voom.rrna-var"
        } else {
                verb("\n\n\nERROR!  unrecognized type=[%s]!!!\n", type)
                stop()
        } # type

        res  <-  list(dir = bdir , tag = btag )
        return(res)

} # de_analysis_type_to_directory_and_tag






read_gtf  <-  function( file , header = FALSE ) {
	message('[read_gtf][file] ', file)
        in.gtf  <-  fread( cmd = sprintf("zcat  -f  %s", file) , sep = "\t" , header = FALSE ,  stringsAsFactors = FALSE , quote = "" , colClasses = c( rep("character",3) , rep("numeric",2) , rep("character",4) ) )

        def.colnames  <-  c("seqname" , "source" , "feature" , "start" , "end" , "score" , "strand" , "frame" , "attribute")
        colnames(in.gtf)  <-  def.colnames

        in.gtf  <-  as.data.frame(in.gtf)

        return(in.gtf)

} # read_gtf




## Start modification by Nitish
load_transcriptome_gtf  <-  function( species , verbose = FALSE ) {
	message('[load_transcriptome_gtf][species] ', species)
        ### transcriptome
	if (species == "mm10") {
          fname  <-  sprintf("data/mouse/mm10/Mus_musculus.GRCm38.97.rdna_rn18s.gtf.gz")
        }

	else if (species == "mm38") {
          fname  <-  sprintf("data/mouse/mm10/Mus_musculus.GRCm38.97.rdna_rn18s.gtf.gz")
        }

	else if (species == "mm39") {
          fname  <-  sprintf("data/mouse/ensembl/Mus_musculus.GRCm39.104.rdna_rn18s.gtf.gz")
        }
	else if (species == "GENCODEmm39") {
          fname  <-  sprintf("data/mouse/gencode/gencode.vM27.annotation.rdna_rn18s.gtf.gz")
        }
 
	else {
                verb("\n\n\nERROR!  unrecognized species=[%s]!!!\n", species)
                stop()
        }
        if (verbose) { verb("\t\t\tloading transcriptome gtf: [%s]\n", fname ) }

	message('[load_transcriptome_gtf][fname] ', fname)


        indf  <-  read_gtf( file = fname )
	
	print(head(indf))

        return(indf)

} # load_transcriptome_gtf






get_gtf_attribute_field  <-  function( gtf , field , num = FALSE ) {
        # always returns as a string
        # num = TRUE -> will convert to numeric (via as.numeric)

        library(stringr)

        fmt  <-  paste( field , ' "[^"]+' , sep="" )
        exf  <-  str_extract( gtf$attribute , fmt )

        fmt  <-  paste( field , ' "' , sep="")
        exf  <-  str_replace( exf , fmt , '' )

        if (num) { exf  <-  as.numeric(exf) }
        return(exf)

} # get_gtf_attribute_field





prepare_transcriptome_gtf  <-  function( gtf ) {

	print(str(gtf))	

        if (!("transcript_id" %in% colnames(gtf))) {
		message('[transcript_id] start')
                gtf$transcript_id  <-  gtf$seqname
		message('[transcript_id] done')
        }
        if (!("gene_name" %in% colnames(gtf))) {
		message('[gene_name] start')
                gtf$gene_name  <-  get_gtf_attribute_field( gtf = gtf , field = "gene_name" )
		message('[gene_name] end')
        }
        if (!("gene_id" %in% colnames(gtf))) {
		message('[gene_id] start')
                gtf$gene_id  <-  get_gtf_attribute_field( gtf = gtf , field = "gene_id" )
		message('[gene_id] end')
        }

        return(gtf)

} # prepare_transcriptome_gtf



convert_ensembl_gene_id_to_gene_name_and_aggregate_matrix  <-  function( mat , species = NULL , fun = sum , gtf = NULL , verbose = FALSE ) {

        if (is.null(gtf)) {
                if (verbose) { verb("\t\tloading gtf.\n") }
                ## Modification by Nitish
		#gtf  <-  load_genome_gtf( species = species )
		gtf  <-  load_transcriptome_gtf( species = species)

                if (verbose) { verb("\t\tpreparing gtf.\n") }
                gtf  <-  prepare_transcriptome_gtf( gtf = gtf )
        } # gtf


        #### table
        if (verbose) { verb("\t\tconvert table.\n") }

        namedf  <-  unique(gtf[,c("gene_name","gene_id"),drop=FALSE])

        if (verbose) { show(head(namedf)) }


        ### melt
        if (verbose) { verb("\t\tmelt.\n") }

        mdf  <-  reshape2::melt(mat)
        colnames(mdf)  <-  c("gene_id" , "col" , "value" )

        if (verbose) { show(head(mdf)) }

        mdf  <-  merge( mdf , namedf , by = "gene_id" , all.x = TRUE )

        if (verbose) { verb("\t\tmerge.d\n") ; show(head(mdf)) }

        # begin of modification by H. Kim
        #mdf$gene_name[is.na(mdf$gene_name)]  <-  mdf$gene_id
        f_no_conversion <- is.na(mdf$gene_name)
        if (any(f_no_conversion)) {
          # ENS1,ENS2 --> sym1,sym2
          mdf[f_no_conversion, "gene_name"] <- sapply(as.character(mdf[f_no_conversion, "gene_id"]), function(x) {
                items <- strsplit(x, ',')[[1]]
                idx <- fmatch(items, gtf$gene_id)
                if (all(!is.na(idx))) {
                  paste(gtf[idx,"gene_name"], collapse=',')
                } else {
                  NA
                }
          })
          f_no_conversion <- is.na(mdf$gene_name)
          mdf[f_no_conversion, "gene_name"]  <-  as.character(mdf[f_no_conversion, "gene_id"])
        }
        # end of modification
	## Modification by Nitish
	mdf <- mdf[!is.na(mdf$gene_name),]
	#mdf <- mdf[!grep("^NA$", mdf$gene_name),]
        mdf  <-  mdf  %>%  group_by(gene_name,col)  %>%  summarise( total.val = fun(value) )
        mdf  <-  as.data.frame(mdf)


        if (verbose) { verb("\t\tcast.\n") ; show(head(mdf)) }

        amat  <-  acast( data = mdf , formula = gene_name ~ col , value.var = "total.val" )

        return(amat)

} # convert_ensembl_gene_id_to_gene_name_and_aggregate_matrix




make_design_matrix  <-  function( info , condition = "condition" , nuisance = NULL , verbose=FALSE ) {

        vars  <-  c(condition , nuisance)
        for (varx  in  vars) {
                stopifnot(varx %in% colnames(info))

                if (!is.factor(info[[varx]])) {
                        stopifnot(is.character(info[[varx]]))
                        info[[varx]]  <-  factor( info[[varx]] , levels = unique(info[[varx]]) )
                } # info
        } # varx

        dform  <-  as.formula(paste( "~0 " , paste(vars , collapse = "+" ) , sep = "+" ))
        if (verbose) { verb("\t\tdormf = %s\n", dform) }
        design  <-  model.matrix( object = dform , data = info )

        colnames(design)  <-  str_replace( colnames(design) , paste0("^",condition) , "" )

#       rownames(design)  <-  rownames(info)

        return(design)

} # make_design_matrix




make_contrast_matrix  <-  function( info , condition = "condition" , design ) {

        uniqConds  <-  unique(info[[condition]])

        contrastList <- list()
        for (i  in  1:length(uniqConds) ) {
                for (j  in  1:length(uniqConds) ) {
                        if (i < j) {
                                contrastList <- append( contrastList , sprintf("%s-%s",  uniqConds[j] , uniqConds[i] ))
                        } # i < j
                } # j
        } # i

        contrast.matrix <- makeContrasts( contrasts = contrastList , levels = design )


        return(contrast.matrix)

} # make_contrast_matrix






write_limma_de_tables  <-  function( cfit.fname , fit ,  info , condition = "condition" , group = NULL , triples = TRUE , MHadjust , FDR.thresh = 0.05 , th_log2fc = 0, map = NULL , outfbase  )  {

        if (!is.null(map)) {
                stopifnot(is.data.table(map))
                stopifnot(all(c("key","value") %in% colnames(map)))
                colnames(map)[colnames(map) == "value"]  <-  "mapped_gene"
        } #

        # cfit has to be re-loaded so that it has FDR values
        inCfit <- read.table( file = cfit.fname ,  row.names = NULL , header = TRUE , sep="\t" , stringsAsFactors = FALSE , quote = "" )
        ### Modification by Nitish; one genes name is missong <NA>. Remove rows with missing names
	inCfit <- inCfit[!is.na(inCfit$genes),]
	rownames(inCfit) <- inCfit$genes

        myExpr <- fit$coefficients
        myDev <- fit$stdev.unscaled / fit$sigma
        theGenes <- rownames(inCfit)

        uniqConds  <-  unique(info[["condition"]])

        alldedf  <-  data.frame()
        allrevdf  <-  data.frame()

        for (i  in  1:length(uniqConds) ) {
                for (j  in  1:length(uniqConds) ) {
                        if (i < j) {
                                cname1 <- as.character( uniqConds[i] )
                                cname2 <- as.character( uniqConds[j] )

                                if (!is.null(group)) {
                                        g1  <-  unique(info[[group]][info[[condition]] == cname1])
                                        g2  <-  unique(info[[group]][info[[condition]] == cname2])
                                        stopifnot(length(g1) == 1  &&  length(g2)==1)
                                        if (g1 != g2) { next }
                                } # group

                                verb("\t\t\t%s   %s\n", cname1,cname2)

                                pvalfield <- sprintf("p.value.%s.%s", cname2,cname1)
                                if (length(uniqConds) == 2) {  pvalfield  <-  "p.value"  }

                                if (MHadjust == "none") {
                                        FDRfield <- sprintf("p.value.%s.%s", cname2,cname1)
                                        if (length(uniqConds) == 2) {  FDRfield  <-  "p.value"  }
                                } else {
                                        FDRfield <- sprintf("p.value.adj.%s.%s", cname2,cname1)
                                        if (length(uniqConds) == 2) {  FDRfield  <-  "p.value.adj"  }
                                } # FDR
                                # expression
                                outData  <-  as.data.frame(fit$coefficients[ theGenes , c(cname1,cname2) ])

                                # FC and FDR
                                outData$log2FC  <-  outData[[cname2]] - outData[[cname1]]
                                outData$FDR  <-  inCfit[theGenes , FDRfield]
                                outData$p.value  <-  inCfit[theGenes , pvalfield]
                                outData  <-  outData[ order(outData$FDR) , ,drop=FALSE]

                                # long format
                                subdf  <-  as.data.frame(outData)
                                subdf$condition1  <-  cname1
                                subdf$condition2  <-  cname2
                                subdf$expression1  <-  subdf[[cname1]]
                                subdf$expression2  <-  subdf[[cname2]]
                                subdf$gene  <-  rownames(subdf)
                                subdf  <-  subdf[, c("gene","condition1","condition2","expression1","expression2","log2FC","FDR","p.value")]
                                alldedf  <-  rbind( alldedf , subdf )



                                # all
                                fname <- sprintf("%s.%s--vs--%s.all.txt",outfbase,cname1,cname2)
                                write.table( outData , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )

                                if (!is.null(map)) {
                                        subod  <-  outData
                                        mdt  <-  as.data.table(subod)
                                        mdt[, gene := rownames(subod)]
                                        mdt  <-  merge( mdt , map , by.x = "gene" , by.y = "key" , all.x = TRUE )
                                        fname <- sprintf("%s.%s--vs--%s.all.mapped.txt",outfbase,cname1,cname2)
                                        fwrite(mdt , file=fname, sep="\t", col.names = T)
                                } # map

                                # diff
                                outData  <-  outData[ order(outData$log2FC) , ,drop=FALSE]
                                logi.de <-  outData$FDR  <=  FDR.thresh  &  !is.na(outData$FDR)
                                logi.de <-  logi.de & (abs(outData$log2FC) > th_log2fc)
                                logi.up <- outData$log2FC  >  th_log2fc
                                logi.down <- outData$log2FC  <  -th_log2fc

                                fname <- sprintf("%s.%s--vs--%s.diff-all.txt",outfbase,cname1,cname2)
                                write.table( outData[logi.de, , drop=FALSE] , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )

                                if (!is.null(map)) {
                                        subod  <-  outData[logi.de, , drop=FALSE]
                                        mdt  <-  as.data.table(subod)
                                        mdt[, gene := rownames(subod)]
                                        mdt  <-  merge( mdt , map , by.x = "gene" , by.y = "key" , all.x = TRUE )
                                        fname <- sprintf("%s.%s--vs--%s.diff-all.mapped.txt",outfbase,cname1,cname2)
                                        fwrite(mdt , file=fname, sep="\t", col.names = T)
                                } # map

                                # up
                                fname <- sprintf("%s.%s--vs--%s.diff-up.txt",outfbase,cname1,cname2)
                                write.table( outData[logi.up &  logi.de, ,drop=FALSE] , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )

                                if (!is.null(map)) {
                                        subod  <-  outData[logi.up &  logi.de, , drop=FALSE]
                                        mdt  <-  as.data.table(subod)
                                        mdt[, gene := rownames(subod)]
                                        mdt  <-  merge( mdt , map , by.x = "gene" , by.y = "key" , all.x = TRUE )
                                        fname <- sprintf("%s.%s--vs--%s.diff-up.mapped.txt",outfbase,cname1,cname2)
                                        fwrite(mdt , file=fname, sep="\t", col.names = T)
                                } # map

                                # down
                                fname <- sprintf("%s.%s--vs--%s.diff-down.txt",outfbase,cname1,cname2)
                                write.table( outData[logi.down  &  logi.de, ,drop=FALSE] , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )

                                if (!is.null(map)) {
                                        subod  <-  outData[logi.down &  logi.de, , drop=FALSE]
                                        mdt  <-  as.data.table(subod)
                                        mdt[, gene := rownames(subod)]
                                        mdt  <-  merge( mdt , map , by.x = "gene" , by.y = "key" , all.x = TRUE )
                                        fname <- sprintf("%s.%s--vs--%s.diff-down.mapped.txt",outfbase,cname1,cname2)
                                        fwrite(mdt , file=fname, sep="\t", col.names = T)
                                } # map

                                for (k  in  1:length(uniqConds) ) {
                                        if (!triples) { next }
                                        if ( j < k) {
                                                cname3 <- as.character( uniqConds[k] )

                                                if (!is.null(group)) {
                                                        g1  <-  unique(info[[group]][info[[condition]] == cname1])
                                                        g3  <-  unique(info[[group]][info[[condition]] == cname2])
                                                        stopifnot(length(g1) == 1  &&  length(g3)==1)
                                                        if (g1 != g3) { next }
                                                } # group

                                                verb("\t\t\t\t%s   %s  %s\n", cname1,cname2,cname3)


                                                pvalfield23 <- sprintf("p.value.%s.%s", cname3,cname2)

                                                if (MHadjust == "none") {
                                                        FDRfield23 <- sprintf("p.value.%s.%s", cname3,cname2)
                                                } else {
                                                        FDRfield23 <- sprintf("p.value.adj.%s.%s", cname3,cname2)
                                                } # FDR



                                                log2FC12 <- sprintf("log2FC%sVS%s",cname1,cname2)
                                                log2FC23 <- sprintf("log2FC%sVS%s",cname2,cname3)
                                                FDR12 <- sprintf("FDR%sVS%s",cname1,cname2)
                                                FDR23 <- sprintf("FDR%sVS%s",cname2,cname3)
                                                pval12 <- sprintf("p.value%sVS%s",cname1,cname2)
                                                pval23 <- sprintf("p.value%sVS%s",cname2,cname3)

                                                outData  <-  as.data.frame(myExpr[theGenes, c(cname1,cname2,cname3)])
                                                outData[[log2FC12]]  <-  outData[[cname2]]  -  outData[[cname1]]
                                                outData[[log2FC23]]  <-  outData[[cname3]]  -  outData[[cname2]]
                                                outData[[FDR12]] <-  inCfit[theGenes , FDRfield]
                                                outData[[FDR23]] <-  inCfit[theGenes , FDRfield23]
                                                outData[[pval12]] <-  inCfit[theGenes , pvalfield]
                                                outData[[pval23]] <-  inCfit[theGenes , pvalfield23]
                                                outData  <-  outData[ order(pmin(outData[[FDR12]] , outData[[FDR23]])) , ,drop=FALSE]

                                                # long format
                                                subdf  <-  as.data.frame(outData)
                                                subdf$condition1  <-  cname1
                                                subdf$condition2  <-  cname2
                                                subdf$condition3  <-  cname3
                                                subdf$expression1  <-  subdf[[cname1]]
                                                subdf$expression2  <-  subdf[[cname2]]
                                                subdf$expression3  <-  subdf[[cname3]]
                                                subdf$log2FC12  <-  outData[[log2FC12]]
                                                subdf$log2FC23  <-  outData[[log2FC23]]
                                                subdf$FDR12  <-  outData[[FDR12]]
                                                subdf$FDR23  <-  outData[[FDR23]]
                                                subdf$p.value12  <-  outData[[pval12]]
                                                subdf$p.value23  <-  outData[[pval23]]
                                                subdf$gene  <-  rownames(subdf)
                                                subdf  <-  subdf[, c("gene","condition1","condition2","condition3","expression1","expression2","expression3","log2FC12","log2FC23","FDR12","FDR23","p.value12","p.value23")]
                                                allrevdf  <-  rbind( allrevdf , subdf )


                                                # all
                                                fname <- sprintf("%s.%s--vs--%s--vs--%s.all.txt",outfbase,cname1,cname2,cname3)
                                                write.table( outData , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )

                                                # rev
                                                outData  <-  outData[ order(outData[[log2FC12]]) , ,drop=FALSE]
                                                logi.de <-  outData[,FDR12]  <= FDR.thresh   &  outData[,FDR23]  <= FDR.thresh  &  !is.na(outData[[FDR12]])  &  !is.na(outData[[FDR23]])
                                                logi.de <- logi.de & (abs(outData[,log2FC12]) > th_log2fc) & (abs(outData[,log2FC23]) > th_log2fc)
                                                logi.ud <- outData[,log2FC12] > th_log2fc   &   outData[,log2FC23] < -th_log2fc
                                                logi.du <- outData[,log2FC12] < -th_log2fc   &   outData[,log2FC23] > th_log2fc
                                                logi.rev  <-  logi.de  &  (logi.ud  |  logi.du)

                                                # rev all
                                                fname <- sprintf("%s.%s--vs--%s--vs--%s.rev-all.txt",outfbase,cname1,cname2,cname3)
                                                write.table( outData[logi.rev, ,drop=FALSE] , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )

                                                # upDown
                                                fname <- sprintf("%s.%s--vs--%s--vs--%s.rev-upDown.txt",outfbase,cname1,cname2,cname3)
                                                write.table( outData[logi.ud &  logi.rev, ,drop=FALSE] , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )

                                                # downUp
                                                fname <- sprintf("%s.%s--vs--%s--vs--%s.rev-downUp.txt",outfbase,cname1,cname2,cname3)
                                                write.table( outData[logi.du  &  logi.rev, ,drop=FALSE] , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )
                                        } # j < k
                                } # k
                        } # i < j
                } # j
        } # i

        fname <- sprintf("%s.de-all.txt",outfbase,cname1,cname2)
        write.table( alldedf , file = fname , row.names = FALSE , col.names = TRUE ,  sep="\t"  , quote = FALSE )

        if (!is.null(map)) {
                mdt  <-  merge( as.data.table(alldedf) , map , by.x = "gene" , by.y = "key" , all.x = TRUE )
                fname <- sprintf("%s.de-all.mapped.txt",outfbase,cname1,cname2)
                fwrite(mdt , file=fname, sep="\t", col.names = T)
        } # map

        if (triples) {
                fname <- sprintf("%s.rev-all.txt",outfbase,cname1,cname2)
                write.table( allrevdf , file = fname , row.names = FALSE , col.names = TRUE ,  sep="\t"  , quote = FALSE )
        } # triples


} # write_limma_de_tables












####################################
#################################### mrna
####################################

# limma_voom_for_mrna
# input:
#   normMethod: c("TMM","TMMwzp","RLE","upperquartile","none")
# comment:
limma_voom_for_mrna  <-
		function(	project ,
				individual = NULL ,
				conditions = NULL ,
				nuisance = NULL ,
				strategy = NULL ,
				mtx_count = NULL ,
				normMethod = "TMM" ,
				doVoom = TRUE ,
				min.expr.counts = 5 ,
				min.expr.log2cpm = -Inf ,
				min.expr.num.samples = 2 ,
				ebayes.trend = TRUE ,
				ebayes.robust = TRUE ,
				weights = TRUE ,
				triples = FALSE ,
				MHadjust = "BH" ,
				FDR.thresh = 0.05 ,
                                th_log2fc = 0 ,
				refColumn = NULL ,
				rundate_appendix = "" ,
				outfbase = NULL , 
				gtf = NULL ,
				verbose = FALSE ) {
	
	##################### parameters
	
	
	
	######################## load
	verb("load.\n")
	
	
	
	############ project
	verb("\tproject.\n")
	
	myData.project  <-  process_project_input( project = project , individual = individual )
	rownames(myData.project)  <-  myData.project$runid
	
	myData.project  <-  myData.project[ order(myData.project$order , myData.project$runid) , ,drop=FALSE]
	
	if (is.null(individual)) {
		individual  <-  myData.project$individual[1] 
	}
	
	if (!is.null(strategy)) {
		myData.project  <-  myData.project[ myData.project$strategy == strategy , ,drop=FALSE]
	} # strategy
	
	if (!is.null(conditions)) {
		myData.project  <-  subset( myData.project , condition %in% conditions )
	} # conditions
	
	
	# basic info	
	myData.species  <-  myData.project$species[1]
	myData.strategy  <-  myData.project$strategy[1]
	myData.DNAlevel  <-  get_DNA_level_for_strategy( strategy = myData.strategy )
	
	# check
	stopifnot( length(unique(myData.project$strategy)) == 1 )
	stopifnot( length(unique(myData.project$species)) == 1 )
	stopifnot( length(unique(myData.project$individual)) == 1 )
	
	
	
	
	
	################# default output
	if (is.null(outfbase) ) {
		projname  <-  paste( unique(myData.project$project) , collapse = "_" )
                list_out <- de_analysis_type_to_directory_and_tag(type="mrna",  project)
		bdir <- list_out$dir
		odir  <-  sprintf("%s/%s%s/%s" , bdir , projname , rundate_appendix , individual )
		otag  <-  sprintf("%s.%s.%s" , projname , individual, list_out$tag)
		outfbase  <-  sprintf("%s/%s", odir, otag)
	} # outfbase
	
	outdir  <-  dirname(outfbase)
	dir.create(outdir , recursive = TRUE , showWarnings = FALSE )
	
	
	##### gtf
	verb("\t\ttranscriptome gtf.\n")
	
	if (is.null(gtf)) {
		myData.trangtf  <-  load_transcriptome_gtf( species = myData.species )
	} else {
		myData.trangtf  <-  gtf
	}
	myData.trangtf  <-  prepare_transcriptome_gtf( gtf = myData.trangtf )
	
	
	
	
	
	############################# counts
	verb("\tcounts.\n")
	
        if (is.null(mtx_count)) {
	  mtx_count  <-  load_raw_counts_and_make_matrix_project( project = myData.project , type = "mrna", legacy = FALSE )
        }
	
	fname  <-  sprintf("%s.counts.raw.gene_id.txt", outfbase )
	write.table( mtx_count , file = fname , sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	### convert to gene names
        f_meta_tags <- grepl("^__|rDNA_promoter", rownames(mtx_count))
	if (any(f_meta_tags)) {
            mtx_meta  <-  mtx_count[f_meta_tags,,drop=F]
            mtx_count  <-  mtx_count[!f_meta_tags,]
	}
        mtx_count  <-  convert_ensembl_gene_id_to_gene_name_and_aggregate_matrix( mat = mtx_count , species = myData.species , fun = sum , gtf = myData.trangtf )
        if (any(f_meta_tags)) {
            mtx_count <- rbind(mtx_meta, mtx_count)
        }
	
	fname  <-  sprintf("%s.counts.raw.txt", outfbase )
	write.table( mtx_count , file = fname , sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	
	
	
	
	
	
	
	########################### prepare
	verb("\tprepare.\n")
	
	verb("count matrix dimension:\n")
	show(dim(mtx_count))
	
	
	
	
	
	
	################### design
	verb("\n\ndesign.\n")

	design  <-  make_design_matrix( info = myData.project , condition = "condition" , nuisance = nuisance , verbose=verbose )
	
	show(design)
	write.table( design , file = sprintf("%s.design.txt",outfbase) , sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	
	
	
	
	###################### contrast
	verb("\n\ncontrast.\n")
	
	contrast.matrix  <-  make_contrast_matrix( info = myData.project , condition = "condition" , design = design )
	
	write.table( contrast.matrix , file = sprintf("%s.contrast.txt",outfbase) , sep="\t" ,  row.names=TRUE , col.names=NA , quote = FALSE  )
	
	
	
	
	
	
	
	
	
	########################## count preparation
	verb("count preparation.\n")
	
        # do not remove __no_feature, __ambiguous, __too_low_aQual, __not_aligned, __alignment_not_unique, since # of total reads need to be transfered to DGEList, but the meta tag rows from the DGEList object should be removed before downstream analyses. https://www.biostars.org/p/379690/,https://support.bioconductor.org/p/104972/
	# https://rdrr.io/bioc/edgeR/man/DGEList.html
	# default: lib.size = colSums(counts): numeric vector giving the total count (sequence depth) for each library.

	myData.counts.dge <- DGEList( counts = mtx_count  ,  group = myData.project$condition  ,  genes = row.names(mtx_count) )
	
	
	### filter
	verb("\tfiltering.\n")
	
	logi.good.count.expr  <-  rowSums(getCounts(myData.counts.dge) >= min.expr.counts)  >=  min.expr.num.samples
	
	logi.good.log2cpm.expr  <-  rowSums(cpm(myData.counts.dge , log = TRUE) >= min.expr.log2cpm) >= min.expr.num.samples
	
	logi.is.expr  <-  logi.good.log2cpm.expr  &  logi.good.count.expr

        # begin of addition by H. Kim
        f_meta_tags <- grepl("^__|rDNA_promoter", rownames(myData.counts.dge))
	logi.is.expr <- logi.is.expr & !f_meta_tags
	logi.is.nonexpr <- !logi.is.expr & !f_meta_tags
	ngene <- length(logi.is.expr) - length(which(f_meta_tags))
        # end of addition
	
        # begin of modification by H. Kim
	#verb("\t\t\tnon-expressed genes:  %d  /  %d   =  %2.2f%%\n" , sum(!logi.is.expr) , length(logi.is.expr) , sum(!logi.is.expr) / length(logi.is.expr) * 100 )
	verb("\t\t\tnon-expressed genes:  %d  /  %d   =  %2.2f%%\n" , sum(logi.is.nonexpr) , ngene , sum(logi.is.nonexpr) / ngene * 100 )
	
	### save non-expressed
	#non.expr.names  <-  myData.counts.dge$genes$gene[!logi.is.expr]
	non.expr.names  <-  myData.counts.dge$genes$gene[logi.is.nonexpr]
	# end of modification

	myData.nonexpr.mat  <-  mtx_count[ non.expr.names ,]
	
	fname  <-  sprintf("%s.nonexpr.txt" , outfbase)
	write.table( myData.nonexpr.mat , file = fname , sep="\t" , row.names = TRUE , col.names = NA , quote = FALSE )
	
	
	
	#### expressed genes only
        # the meta tag rows from the DGEList object should be removed before downstream analyses. https://www.biostars.org/p/379690/,https://support.bioconductor.org/p/104972/
	myData.counts.filt <- myData.counts.dge[logi.is.expr,, keep.lib.sizes=FALSE]
	
	
	
	################################### normalization
	
	verb("\n\n\nnormalization: [%s]\n" , normMethod )
	myData.counts.filt.normalized <- calcNormFactors(myData.counts.filt , method = normMethod , refColumn = refColumn)
	
     	# save
  	normFact.mat  <-  matrix( myData.counts.filt.normalized$samples$norm.factors , nrow = length(myData.counts.filt.normalized$samples$norm.factors) , ncol = 1 , dimnames = list( colnames(myData.counts.filt.normalized$counts) , "normFactor") )
	      write.table( normFact.mat  , file = sprintf("%s.normFactors.txt",outfbase) , sep="\t" , row.names = TRUE , col.names = NA  , quote = FALSE  )
	
	
	lmat  <-  cpm( myData.counts.filt.normalized , normalized.lib.sizes=TRUE, log=TRUE )
	fname  <-  sprintf("%s.log2cpm.samples.txt", outfbase )
	write.table( lmat , file = fname , sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	
	
	
	################################# voom stabilization
	
        pdf( sprintf("%s.limma.plots.pdf",outfbase))
	if (doVoom) {
	        verb("\n\nvoom stabilization.\n")
		verb("\tchecking voom.\n")
		the.Elist.voom  <-  voom( counts = myData.counts.filt.normalized , design = design , plot = TRUE )
		
		verb("\tvoomWithQualityWeights.\n")
		the.Elist.voom.qual  <-  voomWithQualityWeights( counts = myData.counts.filt.normalized , design = design , plot = TRUE , normalize.method = "none" , maxiter = 100 )
	
		if (weights) {	
			the.Elist  <-  the.Elist.voom.qual
		} else {
			verb("\t\tproceeding WITHOUT weights.\n")
			the.Elist  <-  the.Elist.voom
		} # weights
			
	} else {
		the.Elist  <-  list( E = as.matrix(myData.counts.filt.normalized) , genes = data.frame(genes = rownames(myData.counts.filt.normalized), stringsAsFactors = FALSE) )
		the.Elist  <-  new("EList" , the.Elist)
	}
	
	
	
	
	
	########################################### linear model analysis
	verb("linear model analysis.\n")
	# Below, use fit for log2cpm expression per condition, use cfit for FDR values.
	
	### fit
	verb("\tfit\n")
	
	fit <- lmFit( object = the.Elist , design = design)
	
	
	
	### contrast fit
	verb("\tcontrast fit.\n")
	
	cfit <- contrasts.fit( fit = fit , contrasts = contrast.matrix )
	
	
	
	
	### ebayes
	verb("\teBayes\n")
	
	cfit <- eBayes(fit = cfit , trend = ebayes.trend , robust = ebayes.robust )
	
	
	
	
	
	
	
	
	
	
	
	
	########################### save to file
	verb("\tsave\n")
	
	write.fit( fit , file = sprintf("%s.fit" , outfbase) , digits = 8 , adjust = MHadjust , method = "separate" , sep = "\t"  )
	write.fit( cfit , file = sprintf("%s.cfit" , outfbase) , digits = 8 , adjust = MHadjust , method = "separate" , sep = "\t"  )

	write.table( fit$coefficients ,  file = sprintf("%s.log2cpm.txt",outfbase)   , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE  )
	write.table( fit$stdev.unscaled * fit$sigma ,  file = sprintf("%s.dev.txt",outfbase)   , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE  )
	
	
	
	
	
	
	
	
	################## write DE tables
	verb("\n\nwrite DE tables.\n")

	tu.gene.map  <-  NULL
	if (myData.species == "mg1655") {
		verb("\t\tloading TU gene map for E.coli.\n")

		regulonDB.dirname  <-  "/home/hkim5/riboprof/data/regulonDB/regulon-9.0"
		regulonDB.tuObject.fname  <-  "tu_objects_tmp.txt"
		# tu objects
		fname  <-  sprintf("%s/%s", regulonDB.dirname , regulonDB.tuObject.fname )
		tuObjects  <-  read.table( file = fname ,  header = FALSE , row.names = NULL , sep = "\t" , quote = "" , na.strings = c("NA","NAN","","-") , stringsAsFactors = FALSE )
		colnames(tuObjects)  <-  c("tu.id" , "numtu" , "tu.left" , "tu.right" , "tu.type" , "obj.class" , "obj.id" , "obj.name" , "obj.left" , "obj.right" , "obj.strand" , "obj.colorclass" , "obj.desc" , "obj.sigma" , "obj.evidence" , "obj.ri" , "obj.type")

		tdf  <-  unique(tuObjects[ tuObjects$obj.class == "GN" , ])
		tdf  <-  tdf[,c("tu.id","obj.name")]
		colnames(tdf)  <-  c("key","value")
		tdf  <-  unique(tdf)
		tu.gene.map  <-  copy(as.data.table(tdf))
		show(tu.gene.map)
	} # species

	
	write_limma_de_tables(	cfit.fname = sprintf("%s.cfit" , outfbase) , fit = fit , info = myData.project , condition = "condition" , triples = triples , MHadjust = MHadjust , FDR.thresh = FDR.thresh , th_log2fc = th_log2fc , outfbase = outfbase , map = tu.gene.map ) # th_log2fc argument was added by H. Kim
	
	system(sprintf("gzip  -f  %s.*.txt", outfbase))
	system(sprintf("gzip  -f  %s*.cfit", outfbase))
	system(sprintf("gzip  -f  %s*.fit", outfbase))
	
	
	
	
	
	
	
	
	################# limma plots
	verb("limma plots.\n")
	
	
	
	test.res  <-  limma::decideTests( cfit ,  adjust = MHadjust , method = "separate" , p.value = FDR.thresh )
	
	# expression values from all samples
	limma::plotDensities( the.Elist , log = TRUE  , main = "density of all log2cpm gene expression from all samples" )
	
	
	
	# MA,MD
	limma::plotMA(fit , main = "MA plot of fit")
	limma::plotMD(fit , main = "MD plot of fit")
	
	# MDS
	limma::plotMDS( the.Elist , main = "MDS plot")
	
	# heat
	limma::heatDiagram( test.res , coef = cfit$coefficients  )
	
	# venn
	uniqConds  <-  unique(myData.project$condition)
	if (length(uniqConds)  <= 3) {
		verb("\t\tvenn.\n")
		for  (direc  in  c("both","up","down"))  {
			limma::vennDiagram( test.res , include = direc )
		} # direc
	} # venn
	
	dummyPrint  <-  dev.off()
	
	
	
	
	
	
	
	
	
	verb("\n\n\nDone with differential expression analysis.\n\n\n")
	

} # limma_voom_for_mrna




verb <- function(...) cat(sprintf(...), sep='', file=stdout())










####################################
#################################### rrna variants
####################################


# lmfit_method and th_log2fc arguments were added by H. Kim
limma_voom_for_rrna_variants  <-  
		function(	project ,
				individual = NULL ,
				strategy = NULL ,
				auxil.projects = NULL ,
				lmfit_method = "robust" ,
				ebayes.trend = FALSE ,
				ebayes.robust = FALSE ,
				MHadjust = "BH" ,
				FDR.thresh = 0.05 ,
				th_log2fc = 0 ,
				min.depth = 10 ,
				min.AF = 1e-3 ,
				triples = TRUE ,
				outfbase = NULL ,
				rundate_appendix = ""  ) {

	# examples:
	#		limma_voom_for_rrna_variants( project = "20170224" , individual = "NMuMG" , write = TRUE )
	#		limma_voom_for_rrna_variants( project = "170126" , individual = "NMuMG" , write = TRUE )
	#		limma_voom_for_rrna_variants( project = "170125" , individual = "BalbC" , write = FALSE )
	#		limma_voom_for_rrna_variants( project = "170124" , individual = "MCF7" , write = TRUE )
	
	
	################## DEBUG
	if (FALSE) {
	
		project = "170125" 
		outfbase = "/home/hkim5/riboprof/work/debug-new-rrna-var.again"
		individual = "BalbC"
		strategy = "RNAseq"
		auxil.projects = c("MouseGenomesProject")#,"170126")
		
		ebayes.trend = FALSE 
		ebayes.robust = FALSE
		MHadjust = "BY" 
		FDR.thresh = 0.05 
		th_log2fc = log2(1.2)
		triples = FALSE
		
		#limma_voom_for_rrna_variants( rlist = rlist , outfbase = outfbase)
	} 
	################## DEBUG
	
	
	library(Rpdb)
	
	
	
	######################## PDBs
	
	
	pdb.40s.fname  <-  "/home/hkim5/riboprof/data/pdb/Hsa_40S_RNA_full.pdb.gz"
	pdb.60s.fname  <-  "/home/hkim5/riboprof/data/pdb/Hsa_60S_RNA_full.pdb.gz"
	
	
	
	
	
	
	########################## parameters
	myParams.print.ggtheme  <-  theme_bw(base_size = 18)
	
	
	
	
	
	
	
	
	######################## load
	verb("load.\n")
	
	
	
	
	############ project
	verb("\tproject.\n")
	
	myData.project  <-  process_project_input( project = project , individual = individual )
	rownames(myData.project)  <-  myData.project$runid
	
	myData.project  <-  myData.project[ order(myData.project$order , myData.project$runid) , ,drop=FALSE]
	
	if (is.null(individual)) {
		individual  <-  myData.project$individual[1] 
	}
	
	if (!is.null(strategy)) {
		myData.project  <-  myData.project[ myData.project$strategy == strategy , ,drop=FALSE]
	} # strategy
	
	
	# basic info	
	myData.species  <-  myData.project$species[1]
	myData.strategy  <-  myData.project$strategy[1]
	myData.DNAlevel  <-  get_DNA_level_for_strategy( strategy = myData.strategy )
	
	# check
	stopifnot( length(unique(myData.project$strategy)) == 1 )
	stopifnot( length(unique(myData.project$species)) == 1 )
	stopifnot( length(unique(myData.project$individual)) == 1 )
	
	
	
	
	################# default output
	
	if (is.null(outfbase) ) {
		bdir  <-  "/home/hkim5/riboprof/out/limma-voom.rrna-var"
		# begin of modification by H. Kim
		#odir  <-  sprintf("%s/%s/%s" , bdir , project , individual )
		odir  <-  sprintf("%s/%s%s/%s" , bdir , project , rundate_appendix , individual )
		# end of modification
		dir.create(odir , recursive = TRUE , showWarnings = FALSE )
		otag  <-  sprintf("%s.%s.limma-voom.rrna-var" , project , individual)
		outfbase  <-  sprintf("%s/%s", odir, otag)
	} # outfbase
	
	
	
	
	
	
	
	########## project vcfs
	verb("\tproject vcfs.\n")
	
	myData.vcf  <-  load_combined_project_variant_call_vcfs( project = project , level = "runid" )
	
	### annotate
	subinfo  <-  myData.project[ , c("runid","condition") ,drop=FALSE]
	myData.vcf  <-  merge( myData.vcf , subinfo , by.x = "library" , by.y = "runid" , all.x = TRUE )
	
	
	
	
	### brc
	verb("\tproject brc.\n")
	
	myData.brc  <-  load_combined_project_variant_brc_vcfs( project = project , level = "runid" )
	
	
	
	
	### prototype fasta
	verb("\tprototype fasta.\n")
	
	myData.proto.fa  <-  load_prototype_fasta_for_strategy( species = myData.species , strategy = myData.strategy )
	myData.proto.arch  <-  load_prototype_arch_for_strategy( species = myData.species , strategy = myData.strategy )
	
	
	
	
	
	
	### pdb
	verb("\tpdb.\n")
	
	myData.pdb  <-  list()
	
	myData.pdb$small  <-  read.pdb( file = pdb.40s.fname )
	myData.pdb$large  <-  read.pdb( file = pdb.60s.fname  )
	
	
	
	
	
	
	
	
	
	
	
	
	##################################### prepare
	verb("\n\n\nprepare.\n")
	
	
	################ PDB
	verb("\tpdb.\n")
	
	### full pdb
	myData.pdb$both  <-  merge( myData.pdb$small , myData.pdb$large , reindex = FALSE )
	
	### color subunits
	myData.pdb$both$atoms  <-  pdb_default_color_subunit( pdb = myData.pdb$both$atoms , style = "subunit" )
	
	### fix atom IDs
	myData.pdb$both$atoms  <-  pdb_fix_atom_ids( myData.pdb$both$atoms )
	myData.pdb$both$atoms  <-  pdb_fix_resids( myData.pdb$both$atoms )
	
	
	### pdb sequences
	myData.pdb.fa  <-  pdb_get_sequences( myData.pdb$both$atoms  )
	
	
	
	
	#### map prototype rrna to PDB
	verb("\tmap prototype rrna to PDB.\n")
	
	if (myData.DNAlevel == "DNA") {
		archrrna  <-  myData.proto.arch[ myData.proto.arch$score %in% c("lsu","ssu") , ,drop=FALSE]
		subfa  <-  extract_sequences_from_prototype_by_arch( proto = myData.proto.fa , arch = archrrna )
	
		myData.pairmatch.rrna.arch  <-  archrrna
		myData.pairmatch.rrna.fa  <-  subfa
	
	} else if (myData.DNAlevel == "RNA") {
	
		myData.pairmatch.rrna.arch  <-  myData.proto.arch
		myData.pairmatch.rrna.fa  <-  myData.proto.fa
		
	} # level
	
	myData.pairmatch.rrna.pairs  <-  pairwiseMatchSequences( query.fa = myData.pdb.fa , ref.fa = myData.pairmatch.rrna.fa )
	
	
	
	
	
	
	
	
	
	
	
	
	########### prepare VCFs
	verb("prepare VCFs.\n")
	
	#### allele freq and indel
	myData.vcf  <-  prepare_vcf( vcf = myData.vcf )
	myData.brc  <-  prepare_vcf( vcf = myData.brc )
	
	
	### pass
	myData.vcf  <-  myData.vcf[ myData.vcf$filter == "PASS" , ,drop=FALSE]
	
	
	#### correct AF
	myData.vcf  <-  parse_read_depths_from_vcf_DP4( vcf = myData.vcf , groups = c("library") , correct.AF = TRUE )
	myData.brc  <-  parse_read_depths_from_vcf_DP4( myData.brc , groups = "library" , correct.AF = FALSE )
	
	
	#### vcf names
	myData.vcf  <-  make_vcf_names( vcf = myData.vcf)
	myData.brc  <-  make_vcf_names( vcf = myData.brc )
	
	
	
	
	
	############ threshold
	verb("threshold.\n")

	myData.vcf  <-  myData.vcf[ myData.vcf$num.var.forward >= min.depth  &  myData.vcf$num.var.reverse >= min.depth , ,drop=FALSE]
	
	
	####### basic depth
	verb("basic depth.\n")
	
	rdmat  <-  acast( data = myData.brc , formula = vname ~ library , value.var = "ref.depth" , fill = 0 )
	vdmat  <-  acast( data = myData.brc , formula = vname ~ library , value.var = "var.depth" , fill = 0 )
	plusdmat  <-  acast( data = myData.brc , formula = vname ~ library , value.var = "num.var.forward" , fill = 0 )
	minusdmat  <-  acast( data = myData.brc , formula = vname ~ library , value.var = "num.var.reverse" , fill = 0 )
	
	fname  <-  sprintf("%s.count.ref-depth.txt", outfbase)
	write.table( rdmat , fname ,  sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	fname  <-  sprintf("%s.count.var-depth.txt", outfbase)
	write.table( vdmat , fname ,  sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	fname  <-  sprintf("%s.count.var-forward.txt", outfbase)
	write.table( plusdmat , fname ,  sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	fname  <-  sprintf("%s.count.var-reverse.txt", outfbase)
	write.table( minusdmat , fname ,  sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	
	
	
	
	
	
	
	
	good.vnames  <-  c(myData.vcf$vname)
	
	
	##### auxiliary projects
	if (!is.null(auxil.projects)) {
		verb("auxiliary projects.\n")
	
		subarch  <-  myData.proto.arch
		subarch  <-  subarch[ subarch$score %in% c("lsu","ssu") , ,drop=FALSE]	
	
		auxdf  <-  data.frame()
		for (projx  in  auxil.projects) {
			verb("\t%s\n", projx )
	
			prodf  <-  get_project_info( project = projx )
			indf  <-  load_combined_project_variant_call_vcfs( project = projx )
			indf  <-  make_vcf_names(indf)	
	
			# we only need a ujnique set of vnames
			indf  <-  indf[ !duplicated(indf$vname) , ,drop=FALSE]
			show(dim(indf))
	
			mapvcf  <-  map_mixed_library_vcf_to_prototype_coordinates( to.proto = myData.proto.fa , to.strategy = myData.strategy , to.arch = subarch ,
										from.vcf = indf , from.project = prodf )
	
			mapvcf  <-  make_vcf_names(mapvcf)
			good.vnames  <-  c( good.vnames , mapvcf$vname )
		} # projx
	
		good.vnames  <-  unique(good.vnames)
	
	
		######## depth
		verb("\texclude low read depth.\n")
		
		bad.vnames  <-  rownames(rdmat)[ rowSums((rdmat+vdmat) >= min.depth) < ncol(rdmat) ]
		verb("\t\tbad.vnames:  %d\n", length(bad.vnames))
		
		good.vnames  <-  setdiff( good.vnames , bad.vnames )
	
	} # auxil
	
	verb("good.vnames:  %d\n", length(good.vnames))
	
	
	
	
	
	########### impute
	verb("impute.\n")
	
	myData.impute  <-  impute_variants_across_libraries( vcf = myData.vcf , brc = myData.brc , libraries =  myData.project$runid , vnames = good.vnames )

	imat  <-  acast( myData.impute , vname ~ library , value.var = "AF" , fill = 0 )
	# same order as project
	imat  <-  imat[, myData.project$runid  ,drop=FALSE]
	
	fname  <-  sprintf("%s.af.impute.txt", outfbase)
	write.table( imat , fname ,  sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	
	
	
	
	
	
	
	
	
	##### make AF matrix
	verb("make AF matrix.\n")
	
	myData.af  <-  acast( myData.impute , vname ~ library , value.var = "AF" , fill = 0 )
	
	# same order as project
	myData.af  <-  myData.af[, myData.project$runid  ,drop=FALSE]
	
	
	
	######## exclude low AF
	if (!is.null(auxil.projects)) {
		verb("exclude low AF.\n")
		# i.e. must have non-trivial AF for all replicates in at least one condition
		
		splv  <-  split( myData.project , myData.project$condition)
		
		good.vnames  <-  c()
		for (splx  in  splv) {
			num.reps  <-  rowSums(myData.af[, splx$runid , drop=FALSE] >= min.AF)
			new.good  <-  rownames(myData.af)[ num.reps == nrow(splx) ]
			good.vnames  <-  c( good.vnames , new.good )
		} # splx
		good.vnames  <-  unique(good.vnames)
	
	
		verb("\t%d low AF variants excluded.\n" , nrow(myData.af) - length(good.vnames) )
		
		myData.af  <-  myData.af[ good.vnames , ,drop=FALSE]
		
	} # auxil
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	################ reproducible calls
	verb("reproducible calls.\n")
	
	splv  <-  split( myData.project , myData.project$condition)
	
	good.vnames  <-  c()
	for (splx  in  splv) {
		num.reps  <-  rowSums(myData.af[, splx$runid , drop=FALSE] > 0)
		new.good  <-  rownames(myData.af)[ num.reps == nrow(splx) ]
		good.vnames  <-  c( good.vnames , new.good )
	} # splx
	good.vnames  <-  unique(good.vnames)
	
	myData.reprod.af  <-  myData.af[ good.vnames , ,drop=FALSE]
	myData.irreprod.af  <-  myData.af[ setdiff(rownames(myData.af) , good.vnames) , ,drop=FALSE]
	
	
	
	
	#### save
	verb("\tsave.\n")
	
	fname  <-  sprintf("%s.af.reproducible.txt", outfbase)
	write.table( myData.reprod.af , fname ,  sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	fname  <-  sprintf("%s.af.irreproducible.txt", outfbase)
	write.table( myData.irreprod.af , fname ,  sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	
	
	
	
	################# set AF mat
	verb("set AF mat.\n")
	
	myData.af  <-  myData.reprod.af
	
	fname  <-  sprintf("%s.af.txt", outfbase)
	write.table( myData.af , fname ,  sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	
	sv  <-  vdmat[ rownames(myData.af) , ,drop=FALSE]
	sr  <-  rdmat[ rownames(myData.af) , ,drop=FALSE]
	
	fname  <-  sprintf("%s.af.varcounts.txt", outfbase)
	write.table( sv , fname ,  sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	fname  <-  sprintf("%s.af.refcounts.txt", outfbase)
	write.table( sr , fname ,  sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	
	
	
	
	
	
	
	
	
	
	
	############ all variants pdb
	verb("\tall variants pdb.\n")
	
	dummyde  <-  data.frame( gene = rownames(myData.af) , a = -2 , b = -10 , log2FC = 5 , FDR = 0 )
	rownames(dummyde)  <-  rownames(myData.af)
	
	newpdb  <-  make_de_pdb( pdb = myData.pdb$both$atoms , pdb.fa = myData.pdb.fa ,
				rrna.fa = myData.proto.fa , pairs = myData.pairmatch.rrna.pairs ,
				flag = "vcf"  , df = dummyde , conditions = c("a","b") , annot.type = "binary" )
	
	fname  <-  sprintf("%s.all.pdb" , outfbase )
	write.pdb( pdb( atoms = newpdb ) , file = fname )
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	################### design
	verb("\n\n\tdesign.\n")
	
	design  <-  make_design_matrix( info = myData.project , condition = "condition" )
	
	show(design)
	write.table( design , file = sprintf("%s.design.txt",outfbase) , sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
	
	
	
	
	
	
	
	
	######################## contrast
	verb("\tcontrast.\n")
	
	contrast.matrix  <-  make_contrast_matrix( info = myData.project , condition = "condition" , design = design )
	
	write.table( contrast.matrix , file = sprintf("%s.contrast.txt",outfbase) , sep="\t" ,  row.names=TRUE , col.names=NA , quote = FALSE  )
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	########################################### linear model analysis
	verb("\tlinear model analysis.\n")
	
	# Below, use fit for log2cpm expression per condition, use cfit for FDR values.
	
	
	### log2
	myData.laf  <-  log2(myData.af)
	myData.laf[myData.laf == -Inf]  <-  -50
	
	show(dim(myData.laf))
	
	the.Elist  <-  list( E = myData.laf , genes = data.frame(genes = rownames(myData.laf) , stringsAsFactors = FALSE) )
	the.Elist  <-  new("EList" , the.Elist)
	
	
	
	
	### fit
	verb("\t\tfit\n")
	
        # begin of modification by H. Kim
	#fit <- lmFit( object = the.Elist , design = design)
	fit <- lmFit( object = the.Elist , design = design, method = lmfit_method, maxit=200)
	# end of modification
	
	
	
	### contrast fit
	verb("\t\tcontrast fit.\n")
	
	cfit <- contrasts.fit( fit = fit , contrasts = contrast.matrix )
	
	
	
	### ebayes
	verb("\t\teBayes\n")
	
	cfit <- eBayes(fit = cfit , trend = ebayes.trend , robust = ebayes.robust )
	
	
	
	
	
	
	
	################# limma plots
	verb("limma plots.\n")
	
	pdf( sprintf("%s.limma.plots.pdf",outfbase))
	
	
	### reproducibility
	verb("\treproducibility.\n")
	
	for (rx  in  1:nrow(myData.project)) {
		runid  <-  myData.project$runid[rx]
		strategy  <-  myData.project$strategy[rx]
		condition  <-  myData.project$condition[rx]
		verb("\t\t%s\n", runid )
	
		subdf  <-  myData.vcf[ myData.vcf$runid == runid , , drop=FALSE]
		subdf$is.reproducible  <-  subdf$vname  %in%  rownames(myData.af)
	
		gp  <-  ggplot( data = subdf , aes( x = vname , fill = is.reproducible ) )   +
			geom_bar()   +
			theme_bw()   +
			ggtitle(sprintf("Reproducibility of variant allele calls\n%s  %s  %s" , runid , condition , strategy ))
		print(gp)
	} # rx
	
	
	# heatmap
	heatmap.2(	x = myData.af ,
			Colv = FALSE ,
			dendrogram = "row" ,
			trace = "none" ,
			srtCol = 15 ,
			col = colorpanel( n = 256 , low = "white" , high=  "red" ) ,
			main = sprintf("rRNA variant AF") )
	
	
	
	test.res  <-  limma::decideTests( cfit ,  adjust = MHadjust , method = "separate" , p.value = FDR.thresh )
	
	# expression values from all samples
	limma::plotDensities( the.Elist , log = TRUE  , main = "density of all log2cpm gene expression from all samples" )
	
	
	
	# MA,MD
	limma::plotMA(fit , main = "MA plot of fit")
	limma::plotMD(fit , main = "MD plot of fit")
	
	# MDS
	limma::plotMDS( the.Elist , main = "MDS plot")
	
	# heat
	limma::heatDiagram( test.res , coef = cfit$coefficients  )
	
	# venn
	uniqConds  <-  unique(myData.project$condition)
	if (length(uniqConds)  <= 3) {
		verb("\t\tvenn.\n")
		for  (direc  in  c("both","up","down"))  {
			limma::vennDiagram( test.res , include = direc )
		} # direc
	} # venn
	
	
	dummyPrint  <-  dev.off()
	
	
	
	
	
	
	
	########################### save to file
	verb("\tsave\n")
	
	write.fit( fit , file = sprintf("%s.fit" , outfbase) , digits = 8 , adjust = MHadjust , method = "separate" , sep = "\t"  )
	write.fit( cfit , file = sprintf("%s.cfit" , outfbase) , digits = 8 , adjust = MHadjust , method = "separate" , sep = "\t"  )
	
	write.table( fit$coefficients ,  file = sprintf("%s.log2cpm.txt",outfbase)   , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE  )
	write.table( fit$stdev.unscaled * fit$sigma ,  file = sprintf("%s.dev",outfbase)   , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE  )
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	################## write DE tables
#	verb("\n\nwrite DE tables.\n")
#	
#	write_limma_de_tables(  cfit.fname = sprintf("%s.cfit" , outfbase) , fit = fit ,
#				info = myData.project , condition = "condition" ,
#				triples = triples , MHadjust = MHadjust , FDR.thresh = FDR.thresh , th_log2fc = th_log2fc , 
#				outfbase = outfbase ) # th_log2fc argument was added by H. Kim
	
	
	
	
	########## plots
#	verb("plots.\n")
	
#	compatibility.proj.df  <-  myData.project
#	compatibility.proj.df$sname  <-  compatibility.proj.df$runid
	
#	for (i  in  1:length(uniqConds) ) {
#		for (j  in  1:length(uniqConds) ) {
#			if (i < j) {
#				cname1 <- as.character( uniqConds[i] )
#				cname2 <- as.character( uniqConds[j] )
#	
#				condvec  <-  c(cname1,cname2)
#				condstr  <-  paste( condvec , collapse = "--vs--")
#				de.fname  <-  sprintf("%s.%s.all.txt",outfbase , condstr)
#	
#				verb("\t\t\t\t%s\n", condstr)
#	
#				### plot
#				fname  <-  sprintf("%s.%s.de-plots.pdf" , outfbase , condstr )
#				pdf(file = fname)
	
#				### volcano
#				verb("\t\t\t\t\tvolcano.\n")
	
#				volcano_plot( df = de.fname , conditions = condvec , FDR.thresh = FDR.thresh , th_log2fc = th_log2fc) # argement th_log2fc was added by H. Kim
	
	
				#### scatter
#				verb("\t\t\t\t\tscatter.\n")
#	
#				fold_change_scatter_plot( df = de.fname , conditions = condvec , FDR.thresh = FDR.thresh , th_log2fc = th_log2fc) # argement th_log2fc was added by H. Kim
	
	
				#### single entity
#				verb("\t\t\t\t\tsingle entity.\n")
	
#				plot_datapoints_for_single_gene( counts = myData.af , info = compatibility.proj.df , conditions = condvec , df = de.fname , topN = 50 )
	
#				#### pdb
#				verb("\t\t\t\t\tpdb.\n")
#	
#				newpdb  <-  make_de_pdb( pdb = myData.pdb$both$atoms , pdb.fa = myData.pdb.fa ,
#							rrna.fa = myData.proto.fa , pairs = myData.pairmatch.rrna.pairs , 
#							flag = "vcf"  , df = de.fname , conditions = condvec , annot.type = "binary" )
#	
#				fname  <-  sprintf("%s.%s.pdb" , outfbase , condstr )
#				write.pdb( pdb( atoms = newpdb ) , file = fname )
#	
#	
#				dummy.dev  <-  dev.off()
#	
#	
				##################### triple
#				if (triples) {
#					for (k  in  1:length(uniqConds) ) {
#						if ( j < k) {
#							cname3  <- as.character( uniqConds[k] )
#							condvec  <-  c(cname1,cname2,cname3)
#							condstr  <-  paste( condvec , collapse = "--vs--")
#							de.fname  <-  sprintf("%s.%s.all.txt",outfbase , condstr)
#							verb("\t\t\t\t%s\n", condstr)
#	
	
							### plot
#							fname  <-  sprintf("%s.%s.de-plots.pdf" , outfbase , condstr )
#							pdf(file = fname)
	
							
							### volcano
#							verb("\t\t\t\t\tvolcano.\n")
#			
#							volcano_plot( df = de.fname , conditions = condvec , FDR.thresh = FDR.thresh , th_log2fc = th_log2fc) # argement th_log2fc was added by H. Kim
	
	
							#### scatter
#							verb("\t\t\t\t\tscatter.\n")
	
#							fold_change_scatter_plot( df = de.fname , conditions = condvec , FDR.thresh = FDR.thresh , th_log2fc = th_log2fc) # argement th_log2fc was added by H. Kim
				
	
							#### single entity
#							verb("\t\t\t\t\tsingle entity.\n")
	
#							plot_datapoints_for_single_gene( counts = myData.af , info = compatibility.proj.df , conditions = condvec , df = de.fname , topN = 50 )
	
	
							#### pdb
#							verb("\t\t\t\t\tpdb.\n")
#	
#							newpdb  <-  make_de_pdb( pdb = myData.pdb$both$atoms , pdb.fa = myData.pdb.fa ,
#										rrna.fa = myData.proto.fa , pairs = myData.pairmatch.rrna.pairs , 
#										 flag = "vcf"  , df = de.fname , conditions = condvec , annot.type = "binary" )
#				
#							fname  <-  sprintf("%s.%s.pdb" , outfbase , condstr )
#							write.pdb( pdb( atoms = newpdb ) , file = fname )
#	
#							dummy.dev  <-  dev.off()
#	
#						} # j < k
#					} # k
#				} # triples
#			} # i < j
#		} # j
#	} # i
	
#	system(sprintf("gzip  -f  %s.*.txt", outfbase))
#	system(sprintf("gzip  -f  %s.*.pdb", outfbase))
#	
#	system(sprintf("gzip  -f  %s*.cfit", outfbase))
#	system(sprintf("gzip  -f  %s*.fit", outfbase))
#		
	
#	verb("\n\n\nDone with differential expression analysis.\n\n\n")

}  # limma_voom_for_rrna_variants


############Modified by Nitish

write_limma_de_tables  <-  function( cfit.fname , fit ,  info , condition = "condition" , group = NULL , triples = TRUE , MHadjust , FDR.thresh = 0.05 , th_log2fc = 0, map = NULL , outfbase  )  {
  
  if (!is.null(map)) {
    stopifnot(is.data.table(map))
    stopifnot(all(c("key","value") %in% colnames(map)))
    colnames(map)[colnames(map) == "value"]  <-  "mapped_gene"
  } #
  
  # cfit has to be re-loaded so that it has FDR values
  inCfit <- read.table( file = cfit.fname ,  row.names = NULL , header = TRUE , sep="\t" , stringsAsFactors = FALSE , quote = "" )
  ### Modification by Nitish; one genes name is missong <NA>. Remove rows with missing names
  inCfit <- inCfit[!is.na(inCfit$genes),]
  rownames(inCfit) <- inCfit$genes
  
  myExpr <- fit$coefficients
  myDev <- fit$stdev.unscaled / fit$sigma
  theGenes <- rownames(inCfit)
  
  uniqConds  <-  unique(info[["condition"]])
  
  alldedf  <-  data.frame()
  allrevdf  <-  data.frame()
  
  for (i  in  1:length(uniqConds) ) {
    for (j  in  1:length(uniqConds) ) {
      if (i < j) {
        cname1 <- as.character( uniqConds[i] )
        cname2 <- as.character( uniqConds[j] )
        
        if (!is.null(group)) {
          g1  <-  unique(info[["group"]][info[["condition"]] == cname1])
          g2  <-  unique(info[["group"]][info[["condition"]] == cname2])
          stopifnot(length(g1) == 1  &&  length(g2)==1)
          if (g1 != g2) { next }
        } # group
        
        verb("\t\t\t%s   %s\n", cname1,cname2)
        
        pvalfield <- sprintf("P.value.%s.%s", cname2,cname1)
        if (length(uniqConds) == 2) {  pvalfield  <-  "P.value"  }
        
        if (MHadjust == "none") {
          FDRfield <- sprintf("P.value.%s.%s", cname2,cname1)
          if (length(uniqConds) == 2) {  FDRfield  <-  "p.value"  }
        } else {
          FDRfield <- sprintf("P.value.adj.%s.%s", cname2,cname1)
          if (length(uniqConds) == 2) {  FDRfield  <-  "P.value.adj"  }
        } # FDR
        # expression
        outData  <-  as.data.frame(fit$coefficients[ theGenes , c(cname1,cname2) ])
        
        # FC and FDR
        outData$log2FC  <-  outData[[cname2]] - outData[[cname1]]
        outData$FDR  <-  inCfit[theGenes , FDRfield]
        outData$p.value  <-  inCfit[theGenes , pvalfield]
        outData  <-  outData[ order(outData$FDR) , ,drop=FALSE]
        
        # long format
        subdf  <-  as.data.frame(outData)
        subdf$condition1  <-  cname1
        subdf$condition2  <-  cname2
        subdf$expression1  <-  subdf[[cname1]]
        subdf$expression2  <-  subdf[[cname2]]
        subdf$gene  <-  rownames(subdf)
        subdf  <-  subdf[, c("gene","condition1","condition2","expression1","expression2","log2FC","FDR","p.value")]
        alldedf  <-  rbind( alldedf , subdf )
        
        
        
        # all
        fname <- sprintf("%s.%s--vs--%s.all.txt",outfbase,cname1,cname2)
        write.table( outData , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )
        
        if (!is.null(map)) {
          subod  <-  outData
          mdt  <-  as.data.table(subod)
          mdt[, gene := rownames(subod)]
          mdt  <-  merge( mdt , map , by.x = "gene" , by.y = "key" , all.x = TRUE )
          fname <- sprintf("%s.%s--vs--%s.all.mapped.txt",outfbase,cname1,cname2)
          fwrite(mdt , file=fname, sep="\t", col.names = T)
        } # map
        
        # diff
        outData  <-  outData[ order(outData$log2FC) , ,drop=FALSE]
        logi.de <-  outData$FDR  <=  FDR.thresh  &  !is.na(outData$FDR)
        logi.de <-  logi.de & (abs(outData$log2FC) > th_log2fc)
        logi.up <- outData$log2FC  >  th_log2fc
        logi.down <- outData$log2FC  <  -th_log2fc
        
        fname <- sprintf("%s.%s--vs--%s.diff-all.txt",outfbase,cname1,cname2)
        write.table( outData[logi.de, , drop=FALSE] , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )
        
        if (!is.null(map)) {
          subod  <-  outData[logi.de, , drop=FALSE]
          mdt  <-  as.data.table(subod)
          mdt[, gene := rownames(subod)]
          mdt  <-  merge( mdt , map , by.x = "gene" , by.y = "key" , all.x = TRUE )
          fname <- sprintf("%s.%s--vs--%s.diff-all.mapped.txt",outfbase,cname1,cname2)
          fwrite(mdt , file=fname, sep="\t", col.names = T)
        } # map
        
        # up
        fname <- sprintf("%s.%s--vs--%s.diff-up.txt",outfbase,cname1,cname2)
        write.table( outData[logi.up &  logi.de, ,drop=FALSE] , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )
        
        if (!is.null(map)) {
          subod  <-  outData[logi.up &  logi.de, , drop=FALSE]
          mdt  <-  as.data.table(subod)
          mdt[, gene := rownames(subod)]
          mdt  <-  merge( mdt , map , by.x = "gene" , by.y = "key" , all.x = TRUE )
          fname <- sprintf("%s.%s--vs--%s.diff-up.mapped.txt",outfbase,cname1,cname2)
          fwrite(mdt , file=fname, sep="\t", col.names = T)
        } # map
        
        # down
        fname <- sprintf("%s.%s--vs--%s.diff-down.txt",outfbase,cname1,cname2)
        write.table( outData[logi.down  &  logi.de, ,drop=FALSE] , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )
        
        if (!is.null(map)) {
          subod  <-  outData[logi.down &  logi.de, , drop=FALSE]
          mdt  <-  as.data.table(subod)
          mdt[, gene := rownames(subod)]
          mdt  <-  merge( mdt , map , by.x = "gene" , by.y = "key" , all.x = TRUE )
          fname <- sprintf("%s.%s--vs--%s.diff-down.mapped.txt",outfbase,cname1,cname2)
          fwrite(mdt , file=fname, sep="\t", col.names = T)
        } # map
        
        for (k  in  1:length(uniqConds) ) {
          if (!triples) { next }
          if ( j < k) {
            cname3 <- as.character( uniqConds[k] )
            
            if (!is.null(group)) {
              g1  <-  unique(info[["group"]][info[["condition"]] == cname1])
              g3  <-  unique(info[["group"]][info[["condition"]] == cname2])
              stopifnot(length(g1) == 1  &&  length(g3)==1)
              if (g1 != g3) { next }
            } # group
            
            verb("\t\t\t\t%s   %s  %s\n", cname1,cname2,cname3)
            
            
            pvalfield23 <- sprintf("p.value.%s.%s", cname3,cname2)
            if (MHadjust == "none") {
              FDRfield23 <- sprintf("p.value.%s.%s", cname3,cname2)
            } else {
              FDRfield23 <- sprintf("p.value.adj.%s.%s", cname3,cname2)
            } # FDR
            
            
            
            log2FC12 <- sprintf("log2FC%sVS%s",cname1,cname2)
            log2FC23 <- sprintf("log2FC%sVS%s",cname2,cname3)
            FDR12 <- sprintf("FDR%sVS%s",cname1,cname2)
            FDR23 <- sprintf("FDR%sVS%s",cname2,cname3)
            pval12 <- sprintf("p.value%sVS%s",cname1,cname2)
            pval23 <- sprintf("p.value%sVS%s",cname2,cname3)
            
            outData  <-  as.data.frame(myExpr[theGenes, c(cname1,cname2,cname3)])
            outData[["log2FC12"]]  <-  outData[[cname2]]  -  outData[[cname1]]
            outData[["log2FC23"]]  <-  outData[[cname3]]  -  outData[[cname2]]
            outData[["FDR12"]] <-  inCfit[theGenes , FDRfield]
            outData[["FDR23"]] <-  inCfit[theGenes , FDRfield23]
            outData[["pval12"]] <-  inCfit[theGenes , pvalfield]
            outData[["pval23"]] <-  inCfit[theGenes , pvalfield23]
            outData  <-  outData[ order(pmin(outData[["FDR12"]] , outData[["FDR23"]])) , ,drop=FALSE]
            
            # long format
            subdf  <-  as.data.frame(outData)
            subdf$condition1  <-  cname1
            subdf$condition2  <-  cname2
            subdf$condition3  <-  cname3
            subdf$expression1  <-  subdf[[cname1]]
            subdf$expression2  <-  subdf[[cname2]]
            subdf$expression3  <-  subdf[[cname3]]
            subdf$log2FC12  <-  outData[["log2FC12"]]
            subdf$log2FC23  <-  outData[["log2FC23"]]
            subdf$FDR12  <-  outData[["FDR12"]]
            subdf$FDR23  <-  outData[["FDR23"]]
            subdf$p.value12  <-  outData[["pval12"]]
            subdf$p.value23  <-  outData[["pval23"]]
            subdf$gene  <-  rownames(subdf)
            subdf  <-  subdf[, c("gene","condition1","condition2","condition3","expression1","expression2","expression3","log2FC12","log2FC23","FDR12","FDR23","p.value12","p.value23")]
            allrevdf  <-  rbind( allrevdf , subdf )
            
            
            # all
            fname <- sprintf("%s.%s--vs--%s--vs--%s.all.txt",outfbase,cname1,cname2,cname3)
            write.table( outData , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )
            
            # rev
            outData  <-  outData[ order(outData[["log2FC12"]]) , ,drop=FALSE]
            logi.de <-  outData[,FDR12]  <= FDR.thresh   &  outData[,FDR23]  <= FDR.thresh  &  !is.na(outData[["FDR12"]])  &  !is.na(outData[["FDR23"]])
            logi.de <- logi.de & (abs(outData[,log2FC12]) > th_log2fc) & (abs(outData[,log2FC23]) > th_log2fc)
            logi.ud <- outData[,log2FC12] > th_log2fc   &   outData[,log2FC23] < -th_log2fc
            logi.du <- outData[,log2FC12] < -th_log2fc   &   outData[,log2FC23] > th_log2fc
            logi.rev  <-  logi.de  &  (logi.ud  |  logi.du)
            
            # rev all
            fname <- sprintf("%s.%s--vs--%s--vs--%s.rev-all.txt",outfbase,cname1,cname2,cname3)
            write.table( outData[logi.rev, ,drop=FALSE] , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )
            
            # upDown
            fname <- sprintf("%s.%s--vs--%s--vs--%s.rev-upDown.txt",outfbase,cname1,cname2,cname3)
            write.table( outData[logi.ud &  logi.rev, ,drop=FALSE] , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )
            
            # downUp
            fname <- sprintf("%s.%s--vs--%s--vs--%s.rev-downUp.txt",outfbase,cname1,cname2,cname3)
            write.table( outData[logi.du  &  logi.rev, ,drop=FALSE] , file = fname , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE )
          } # j < k
        } # k
      } # i < j
    } # j
  } # i
  
  fname <- sprintf("%s.de-all.txt",outfbase,cname1,cname2)
  write.table( alldedf , file = fname , row.names = FALSE , col.names = TRUE ,  sep="\t"  , quote = FALSE )
  
  if (!is.null(map)) {
    mdt  <-  merge( as.data.table(alldedf) , map , by.x = "gene" , by.y = "key" , all.x = TRUE )
    fname <- sprintf("%s.de-all.mapped.txt",outfbase,cname1,cname2)
    fwrite(mdt , file=fname, sep="\t", col.names = T)
  } # map
  
  if (triples) {
    fname <- sprintf("%s.rev-all.txt",outfbase,cname1,cname2)
    write.table( allrevdf , file = fname , row.names = FALSE , col.names = TRUE ,  sep="\t"  , quote = FALSE )
  } # triples
  
  
} # write_limma_de_tables











######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
###################   helper  ########################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################















































######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
###################   plot  ##########################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################



# begin of modification by H. Kim
#annotate_limma_df  <-  function( df , FDR.thresh = 0.05 , conditions )  {
annotate_limma_df  <-  function( df , FDR.thresh = 0.05 , th_log2fc = 0, conditions )  {
# end of modification

	### load, if necessary
	if (is.vector(df)  &&  is.character(df)  &&  length(df) == 1) {
		df  <-  read.table( file = df , header = TRUE , row.names = 1 ,  sep="\t" , stringsAsFactors = FALSE , quote = "" )	
	} # load

	df  <-  as.data.frame(df)

	if (!("direc"  %in%  colnames(df))) {
		if (length(conditions) == 2) {
	
			logi.sig  <-  df$FDR <= FDR.thresh
			# begin of modification by H. Kim
			#logi.up  <-  df$log2FC > 0
			logi.up  <-  df$log2FC > th_log2fc
			logi.dn  <-  df$log2FC < -th_log2fc
	
			df$direc  <-  "notSig"
			df$direc[logi.sig & logi.up]  <-  "up"
			#df$direc[logi.sig & !logi.up]  <-  "down"
			df$direc[logi.sig & logi.dn]  <-  "down"
			# enf of modification
	
		} else if (length(conditions) == 3) {

			fdr.12  <-  paste( "FDR" ,  conditions[1] , "VS" , conditions[2] , sep="" )
			fdr.23  <-  paste( "FDR" ,  conditions[2] , "VS" , conditions[3] , sep="" )
			l2fc.12  <-  paste( "log2FC" ,  conditions[1] , "VS" , conditions[2] , sep="" )
			l2fc.23  <-  paste( "log2FC" ,  conditions[2] , "VS" , conditions[3] , sep="" )

			stopifnot( all(c(fdr.12,fdr.23,l2fc.12,l2fc.23) %in% colnames(df)) )
	
			logi.sig2  <-  df[[fdr.12]] <= FDR.thresh  &  df[[fdr.23]] <= FDR.thresh
			# begin of modification by H. Kim
			#logi.upDown  <-  df[[l2fc.12]] > 0  &  df[[l2fc.23]] < 0
			#logi.downUp  <-  df[[l2fc.12]] < 0  &  df[[l2fc.23]] > 0
			logi.upDown  <-  df[[l2fc.12]] > th_log2fc  &  df[[l2fc.23]] < -th_log2fc
			logi.downUp  <-  df[[l2fc.12]] < -th_log2fc  &  df[[l2fc.23]] > th_log2fc
			# end of modification
	
			df$direc  <-  "notSig"
			df$direc[logi.sig2 & logi.upDown]  <-  "upDown"
			df$direc[logi.sig2 & logi.downUp]  <-  "downUp"
		
			# max log2FC
			df$max.abs.log2FC  <-  pmax( abs(df[[l2fc.12]]) , abs(df[[l2fc.23]]) )
			df$max.abs.log2FC[df$direc == "downUp"]  <-  -df$max.abs.log2FC[df$direc == "downUp"]
	
			df$min.FDR  <-  pmin(df[[fdr.12]] , df[[fdr.23]])
	
			# place holder
			df$log2FC  <-  df$max.abs.log2FC
			df$FDR  <-  df$min.FDR
		} # num conds
	} # direc
	
	return(df)

} # annotate_limma_df








basic_count_info_limma_df  <-  function( df ) {

	counts  <-  table(df$direc)
	num.tot  <-  sum(counts)
	perc  <-  counts / num.tot * 100

	maxval  <-  NULL
	if ("log2FC"  %in%  colnames(df)) {
		maxval  <-  max(abs(df$log2FC))
	} # log2FC

	retval  <-  list(count = counts , perc = perc , maxval = maxval )

	return(retval)

} # basic_count_info_limma_df







make_info_text  <-  function( info ) {


	direcs  <-  names(info$count)

	str  <-  paste( direcs , " = " , info$count , " (" , formatC(info$perc , digits = 1 , format = "f") , "%)" , sep = "" , collapse = " , ")
	str  <-  sprintf("%s\ntotal = %d\n", str , sum(info$count) )


	return(str)
} # make_info_text






# begin of modification by H. Kim
#volcano_plot  <-  function( df , conditions , FDR.thresh = 0.05 ) {
volcano_plot  <-  function( df , conditions , FDR.thresh = 0.05 , th_log2fc = 0) {
# end of modification

	library(ggrepel)

	df  <-  annotate_limma_df( df , FDR.thresh = FDR.thresh , th_log2fc = th_log2fc , conditions = conditions ) # argement th_log2fc was added by H. Kim
	info  <-  basic_count_info_limma_df( df )
	df$nlog10FDR  <-  -log10(df$FDR)

	df.sig  <-  df[ df$direc != "notSig" , , drop=FALSE]
	df.notsig  <-  df[ df$direc == "notSig" , , drop=FALSE]


	gp  <-  ggplot()
	if (nrow(df.notsig) > 0) {
		gp  <-  gp  +  geom_point( data = df.notsig , aes( x = log2FC , y = nlog10FDR , color = direc ) )
	}
	if (nrow(df.sig) > 0) {
		rdf  <-  df.sig
		rdf  <-  rdf[order(-abs(rdf$log2FC)) ,,drop=FALSE]
		rdf$label  <-  rdf$gene
		rdf$label[1:nrow(rdf) > 10]  <-  NA
		gp  <-  gp  +  geom_point( data = rdf , aes( x = log2FC , y = nlog10FDR , color = direc ) )   +
				geom_text_repel( data = rdf[!is.na(rdf$label),,drop=FALSE] , aes( x = log2FC , y = nlog10FDR , label = label ) )
	}

	gp  <-  gp   +
		geom_hline(yintercept=0)   +
		geom_vline(xintercept=0)   +
		xlim(c(-info$maxval,info$maxval))   +
		theme_bw()   +
		myParams.limmaFuncs.gg.scale.color   +
		ggtitle(sprintf("Volcano plot\n%s\n%s" , paste(conditions,collapse="  vs  ") , make_info_text(info)  ))
	print(gp)

} # volcano_plot









# begin of modification by H. Kim
#fold_change_scatter_plot  <-  function( df , conditions , FDR.thresh = 0.05 ) {
fold_change_scatter_plot  <-  function( df , conditions , FDR.thresh = 0.05 , th_log2fc = 0) {
# end of modification

	df  <-  annotate_limma_df( df , FDR.thresh = FDR.thresh , th_log2fc = th_log2fc , conditions = conditions ) # argement th_log2fc was added by H. Kim
	info  <-  basic_count_info_limma_df( df )

	df.sig  <-  df[ df$direc != "notSig" , , drop=FALSE]
	df.notsig  <-  df[ df$direc == "notSig" , , drop=FALSE]

	if (length(conditions) == 3) {
		xname  <-  paste( "log2FC" ,  conditions[1] , "VS" , conditions[2] , sep="" )
		yname  <-  paste( "log2FC" ,  conditions[2] , "VS" , conditions[3] , sep="" )
		df$FDR  <-  df$min.FDR
	} else {
		xname  <-  conditions[1]
		yname  <-  conditions[2]
	} # 


	gp  <-  ggplot()
	if (nrow(df.notsig) > 0) {
		gp  <-  gp  +  geom_point( data = df.notsig , aes_string( x = xname , y = yname , color = "direc" ) )
	}
	if (nrow(df.sig) > 0) {
		gp  <-  gp  +  geom_point( data = df.sig , aes_string( x = xname , y = yname , color = "direc" ) )
	}

	if (length(conditions) == 3) {
		gp  <-  gp   +
			xlim(c(-info$maxval,info$maxval))   +
			ylim(c(-info$maxval,info$maxval))
	} # conds

	gp  <-  gp   +
		geom_hline(yintercept=0)   +
		geom_vline(xintercept=0)   +
		geom_abline(intercept=0,slope=1)   +
		coord_fixed()   +
		theme_bw()   +
		myParams.limmaFuncs.gg.scale.color   +
		ggtitle(sprintf("log2 fold-change scatter plot\n%s\n%s" ,  paste(conditions,collapse="  vs  ") , make_info_text(info) ))
	print(gp)

} # fold_change_scatter_plot


































# begin of modification
#limma_go_analysis  <-  function( df , species , conditions , outfbase = NULL , FDR.thresh = 0.05 , num.top.GO = 4e2 , abs.log2FC.thresh = 0 )  {
limma_go_analysis  <-  function( df , species , conditions , outfbase = NULL , rundate_appendix = "" , FDR.thresh = 0.05 , num.top.GO = 4e2 , th_log2fc = 0 )  {
# end of modification

	### species 2 letter abbreviation
	org.code  <-  species_to_org_code( species = species )

	verb("\t\t\tabbreviated to %s\n", org.code)


	### conditions string
	all.conds.str  <-  paste( conditions , collapse = "--vs--" )

	### process df
	df  <-  annotate_limma_df( df , FDR.thresh = FDR.thresh , th_log2fc = th_log2fc , conditions = conditions ) # argement th_log2fc was added by H. Kim

	### log2FC thresh
        # begin of modification
	#stopifnot(abs.log2FC.thresh >= 0)
	#df  <- subset( df , abs(log2FC) >= abs.log2FC.thresh )
	stopifnot(th_log2fc >= 0)
	df  <- subset( df , abs(log2FC) >= th_log2fc )
        # end of modification


	uniq.direcs  <-  setdiff( unique(df$direc) , "notSig" )

	df$entrez.id  <-  convert_genes_to_entrez( genes = rownames(df) , species = species )

	uniq.direcs  <-  unique(c( uniq.direcs , "sig" ))
	outdf  <-  list()
	for (direcx  in  uniq.direcs) {
		verb("\t\t%s\n", direcx)

		logi.na  <-  is.na(df$direc)  |  is.na(df$entrez.id)

		if (direcx == "sig") {
			logi.direc  <-  df$direc != "notSig"  &  !logi.na
			logi.other  <-  df$direc == "notSig"  &  !logi.na
		} else {
			logi.direc  <-  df$direc == direcx  &  !logi.na
			logi.other  <-  df$direc != direcx  &  !logi.na
		} # direcx

		de.names  <-  df$entrez.id[logi.direc]		
		other.names  <-  df$entrez.id[logi.other]
		de.names  <-  unique(de.names)
		other.names  <-  unique(other.names)

		### GO
		go.res  <-  goana( de = de.names , universe = c(de.names,other.names) , species = org.code , plot = TRUE )
		go.res$GOID  <-  rownames(go.res)
		go.res  <-  go.res[ order(go.res$P.DE) , ,drop=FALSE]
		go.res  <-  go.res[ 1:min(nrow(go.res),num.top.GO) , ,drop=FALSE]
		go.res  <-  subset( go.res  , P.DE <= 0.01 )

		### save
		if (!is.null(outfbase)) {
			fname  <-  sprintf("%s.%s.%s.GO.txt" , outfbase , all.conds.str , direcx )
			write.table( go.res , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
		} # outfbase

		outdf[[direcx]]  <-  go.res
	} # direcx

	return(outdf)

} # limma_go_analysis






































# begin of modifcation
#plot_differential_trna_gene_heatmap  <-  function( df , bed , level , conditions , FDR.thresh = 0.05 ) {
plot_differential_trna_gene_heatmap  <-  function( df , bed , level , conditions , FDR.thresh = 0.05 , th_log2fc = 0 ) {
# end of modification
if (level != "aminoAcid") {

	### load bed, if necessary
	if (is.vector(bed)  &&  is.character(bed)  &&  length(bed) == 1) {
		bed  <-  read.bed(file = bed , header = TRUE )	
	} # load
	bed$gene  <-  bed$name

	# basic info
	df  <-  annotate_limma_df( df , FDR.thresh = FDR.thresh , th_log2fc = th_log2fc , conditions = conditions ) # argement th_log2fc was added by H. Kim
	info  <-  basic_count_info_limma_df( df )

	df[[level]]  <-  rownames(df)

	if (level == "gene") {
		df  <-  merge( df , bed[,c("gene","antiCodon","aminoAcid")] , by = level , all = TRUE )
	} else if (level == "antiCodon") {
		df  <-  merge( df , bed[,c("antiCodon","aminoAcid")] , by = level , all = TRUE )
	} 

	# assign log2FC
	if (length(conditions) == 3) { 
		df$log2FC  <-  df$max.abs.log2FC
	}

	df$log2FC[df$direc == "notSig"]  <-  0
	logi.na  <-  is.na(df$FDR)
	df$log2FC[logi.na]  <-  0
	df$FDR[logi.na]  <-  1
	df$direc[logi.na]  <-  "notSig"
	


	## exclude stop codons
	# df  <-  df[ !(df$aminoAcid %in% c("Amb","Och","Sup","Sec")) , ,drop=FALSE]

	df$sig.order  <-  ifelse( df$direc != "notSig" , 1 , 2 )
	# get amino acid property, number of genes
	df  <-  merge( df , get_amino_acid_properties() , by = "aminoAcid" , all.x = TRUE )
	df  <-  df[ order(df$sig.order , df$log2FC) , ,drop=FALSE]
	df  <-  df  %>%  group_by(aminoAcid)  %>%  mutate(num.genes = n() , num.antiCodons = n_distinct(antiCodon))
	df  <-  df  %>%  group_by(aminoAcid,antiCodon)  %>%  mutate(num.genes.ac = n())
	df  <-  as.data.frame(df)

	# order by AA property, num genes
	if (level == "gene") {
		#df  <-  df[ order(df$AAprop,-df$num.genes,-df$num.genes.ac,df$antiCodon,df$gene) , ,drop=FALSE]
		df  <-  df[ order(df$AAprop,df$aminoAcid,df$antiCodon,df$gene) , ,drop=FALSE]
		#df  <-  df[ order(df$AAprop,-df$num.genes,-df$num.genes.ac,df$sig.order,df$log2FC) , ,drop=FALSE]
	} else if (level == "antiCodon") {
		df  <-  df[ order(df$AAprop,-df$num.antiCodons,df$antiCodon) , ,drop=FALSE]
	} # level

	df$aminoAcid  <-  factor( df$aminoAcid , levels = unique(df$aminoAcid) )
	df  <-  df  %>%  group_by(aminoAcid)  %>%  mutate( x = row_number() )
	df  <-  as.data.frame(df)

	df$label  <-  df[[level]]
	if (level == "gene") {
		df$label  <-  str_replace( df$label , "^..." , "" )
	}

	textSize  <-  3
	angle  <-  0
	if (level == "gene") {
		textSize  <-  2
		angle  <-  270
	 } # level

	# with NA
	gp  <-  ggplot( data = df , aes( x = x , y = aminoAcid , fill = log2FC , label = label ) )   +
		geom_raster()   +
		scale_fill_gradient2( low = muted("red") , mid = "white", high = muted("green") , midpoint = 0 , na.value = "grey54" , limits = c(-3,3) )   +
		geom_text( hjust = "center" , vjust = "center" , size = textSize , angle = angle )   +
		ggtitle(sprintf("differential tRNA %s expression\n%s", level , paste(conditions,collapse="  vs  ")))   +
		theme_cowplot()   +
		theme(  axis.title.x=element_blank() , axis.text.x=element_blank() , axis.ticks.x=element_blank()  )
	print(gp)

	# with NA, white
	gp  <-  ggplot( data = df , aes( x = x , y = aminoAcid , fill = log2FC , label = label ) )   +
		geom_raster()   +
		scale_fill_gradient2( low = muted("red") , mid = "white", high = muted("green") , midpoint = 0 , na.value = "white" , limits = c(-3,3) )   +
		geom_text( hjust = "center" , vjust = "center" , size = textSize , angle = angle )   +
		ggtitle(sprintf("differential tRNA %s expression\n%s", level , paste(conditions,collapse="  vs  ")))   +
		theme_cowplot()   +
		theme(  axis.title.x=element_blank() , axis.text.x=element_blank() , axis.ticks.x=element_blank()  )
	print(gp)

	# no NA
	subm  <-  df[ !is.na(df$FDR) & !is.na(df$log2FC) ,,drop=FALSE]
	if (nrow(subm) > 0) {
		gp  <-  ggplot( data = subm , aes( x = x , y = aminoAcid , fill = log2FC , label = label ) )   +
			geom_raster()   +
			scale_fill_gradient2( low = muted("red") , mid = "white", high = muted("green") , midpoint = 0 , na.value = "grey54" , limits = c(-3,3) )   +
			geom_text( hjust = "center" , vjust = "center" , size = textSize  , angle = angle )   +
			ggtitle(sprintf("differential tRNA %s expression\n%s", level , paste(conditions,collapse="  vs  ")))   +
			theme_cowplot()   +
			theme(  axis.title.x=element_blank() , axis.text.x=element_blank() , axis.ticks.x=element_blank()  )
		print(gp)
	} # > 0

} # level != aminoAcid
} # plot_differential_trna_gene_heatmap



























#' Plot count datapoints for single entity
#'
#' Prints a dot plot of expression sample datapoints for a single gene across conditions.
#'
#' @param  counts   a matrix of gene expression values (e.g. normalized counts, log2cpm).  rownames = genes.  colnames = sample names.
#' @param  info   info df/matrix about sample, e.g. rlist.
#' @param  conditions   [character vector] a subset of conditions to be plotted.  if given, counts and info will be restricted to these conditions.
#' @param  genes  [char vec] a subset of genes to plot.  if absent, will topN significant genes from df.
#' @param  df   [filename or df] limma differential expression results df or filename.  If give.
#'
#' @return  nothin.g  a plot is printed to current graphics device.
#'
#' @export
# begin of modification
#plot_datapoints_for_single_gene  <-  function( counts , info , conditions = NULL , genes = NULL , df = NULL , topN = 50 )  {
plot_datapoints_for_single_gene  <-  function( counts , info , conditions = NULL , genes = NULL , df = NULL , topN = 50 , FDR.thresh = 0.05 , th_log2fc = 0)  {
# end of modification

	if (!is.null(df)) {
		# basic info
		df  <-  annotate_limma_df( df , FDR.thresh = FDR.thresh , th_log2fc = th_log2fc , conditions = conditions ) # argement th_log2fc was added by H. Kim

		if (is.null(genes)) {
			df  <-  df[ order(df$FDR) , ,drop=FALSE]
			df  <-  df[ df$direc != "notSig"  &  !is.na(df$direc) , ,drop=FALSE]
			sigs  <-  rownames(df)
			sigs  <-  sigs[1:min(length(sigs),topN)]
			genes  <-  sigs
		} # genes
	} # df

	if (!is.null(conditions)) {
		info  <-  info[ info$condition %in% conditions , ,drop=FALSE]
	} else {
		conditions  <-  unique(info$condition)
	} # conditions

	subdf  <-  as.data.frame(counts[ rownames(counts) %in% genes , ,drop=FALSE])
	subdf$gene  <-  rownames(subdf)
	mdf  <-  reshape2::melt(subdf , id.vars = "gene" , measure.vars = info$sname , variable.name = "sname" , value.name = "normcount" )
	mdf  <-  merge( mdf , info[,c("sname","condition")] , by = "sname" , all.x = TRUE )
	mdf  <-  mdf[ mdf$condition  %in%  conditions , ,drop=FALSE]
	mdf$condition  <-  factor( mdf$condition , levels = conditions )

	for (genex  in  genes) {
		verb("\t\t\t%s\n", genex)

		gdf  <-  mdf[ mdf$gene == genex , ,drop=FALSE]

		gp  <-  ggplot( data = gdf , aes_string( x = "condition" , y = "normcount" )  )   +
			geom_jitter(height = 0 , width = 0.15)   +
			theme_bw()   +
			ggtitle(sprintf("%s" , genex))
		print(gp)
	
	} # genex

} # plot_datapoints_for_single_gene










#'  Convert "gene" names to sequence coordinates
#'
#'  Parses the name of a tested entity to obtain sequence or genomic coordinates.
#'
#' @param  df		a dataframe, e.g. from load_DE_tables().  rownames are entity names
#' @param  flag		[vcf, posLength, cov]  string indicating format of the entity name
#' @param  names	if given, will override the rownames as entity names
#'
#' @return  data frame with additional columns, according to the "flag" provided.
#'		if flag == "vcf", then adds columns "chr","position","ref","alt" as per VCF format.
#'		if flag == "posLength" or "cov", then adds columns "chr","start","end" as per BED format.
#'
#' @export
convert_entity_names_to_sequence_coordinates  <-  function( df , flag , names = NULL ) {

	if (is.null(names)) {
		names  <-  rownames(df)
	}

	if (flag == "vcf") {
		# vname to vcf format
		splv  <-  str_split( names , "_" )
		df$chr  <-  sapply( splv , FUN = function(x) x[1] )
		df$position  <-  as.numeric(sapply( splv , FUN = function(x) x[2] ))
		df$ref  <-  sapply( splv , FUN = function(x) x[3] )
		df$alt  <-  sapply( splv , FUN = function(x) x[4] )
	} else if (flag == "posLength") {
		# to bed format
		splv  <-  str_split( names , "\\." )
		df$chr  <-  sapply( splv , FUN = function(x) x[1] )
		df$start  <-  as.numeric(sapply( splv , FUN = function(x) x[2] )) - 1
		df$end  <-  as.numeric(sapply( splv , FUN = function(x) x[3] ))
	} else if (flag == "cov") {
		# to bed format
		splv  <-  str_split( names , "\\." )
		df$chr  <-  sapply( splv , FUN = function(x) x[1] )
		df$start  <-  as.numeric(sapply( splv , FUN = function(x) x[2] )) - 1
		df$end  <-  df$start + 1
	} # flag

	return( df )

} # convert_entity_names_to_sequence_coordinates









#' Make a PDB for differential variants
#'
#' Creates a pdb with differential positions colored by fold change.
#'
#' @param  df		[filename or df] limma differential expression results df or filename.  If give.
#' @param  conditions	[character vector] a subset of conditions to be plotted.  if given, counts and info will be restricted to these conditions.
#' @param  annot.type	[binary, ternary, full] "binary" = red if DE, blue o/w.  ternary = red if up, blue if down, green o/w.  full = gradient color by fold change.  
#' @param  flag		["rrnaVariant","rrnaPosLengthCov"]  will execute special processing via the "df" interpreting "gene names" as rrna variant names, rrna poslength cov, or rrna pos cov.
#' @param  pairs	the result of "pairwiseMatchSequences" with pdb.fa as QUERY and rrna.fa as SUBJECT.
#'
#' @return dataframe   a pdb data frame
#'
#' @export
# begin of modification
#make_de_pdb  <-  function( pdb , pdb.fa , rrna.fa , pairs , flag = NULL , df , conditions , FDR.thresh = 0.05 , annot.type ) {
make_de_pdb  <-  function( pdb , pdb.fa , rrna.fa , pairs , flag = NULL , df , conditions , FDR.thresh = 0.05 , th_log2fc = 0 , annot.type ) {
# end of modificaiton

	# basic info
	df  <-  annotate_limma_df( df , FDR.thresh = FDR.thresh , th_log2fc = th_log2fc , conditions = conditions ) # argement th_log2fc was added by H. Kim
	info  <-  basic_count_info_limma_df( df )

	if (!any(df$direc != "notSig")) {
		return(pdb)
	} # df


	pcols  <-  pdb_get_colors()

	vdf  <-  df[ df$direc != "notSig"  &  !is.na(df$direc) , ,drop=FALSE ]

	if (!is.null(flag)) {
		prepdf  <-  convert_entity_names_to_sequence_coordinates( df = vdf , flag = flag )
	} # flag 

	if (annot.type == "full") {
		logi.down  <-  (vdf$direc == "down"  |  vdf$direc == "downUp")  &  !is.na(vdf$direc)
		logi.up  <-  (vdf$direc == "up"  |  vdf$direc == "upDown")  &  !is.na(vdf$direc)
	
		maxl2  <-  max(abs(vdf$log2FC))
		absdev  <-  abs(vdf$log2FC) / maxl2 * 50

			
		direcdev  <-  rep(50, times = nrow(vdf))
		direcdev[logi.up]  <-  50 + absdev
		direcdev[logi.down]  <-  50 - absdev
		color  <-  direcdev

	} else if (annot.type == "binary") {

		color  <-  rep(pcols["red"] , times = nrow(vdf))

	} else if (annot.type == "ternary") {

		color  <-  rep(pcols["red"] , times = nrow(vdf))
		logi.down  <-  (vdf$direc == "down"  |  vdf$direc == "downUp")  &  !is.na(vdf$direc)
		color[logi.down]  <-  pcols["green"]

	} # annot.type

	# must be within acceptable range
	logi.bad.col  <-  color < 0  |  color > 100  |  is.na(color)
	if (any(logi.bad.col)) {
		verb("\n\n\nERROR!  invalid color values!\n\n")
		show(color[logi.bad.col])
		stop()
	} # check

	# restrict to acceptable range
	color  <-  pmax(color,0)
	color  <-  pmin(color,99)

	# begin of modification by H. Kim
	#show(annot.type)
	#show(hist(color))
	# end of modification

	newpdb  <-  pdb_color_by_annotation( df = prepdf , pdb = pdb , pdb.fa = pdb.fa , seq.fa = rrna.fa , color = color , pairs = pairs , verbose = FALSE )

	return(newpdb)
} # make_de_pdb















# # # 			source("/home/hkim5/riboprof/src/differential_expression_analysis.R")

test_property_enrichment  <-  function( groups ,  prop , outfbase = NULL , rundate_appendix = "" , plot = FALSE , topN = 10 ) {
### DESCRIPTION
#	This function will test various combinations of the groupings for differential property enrichment via the finite median test
#
### INPUT
#	groups	= a named vector.  names are ids.  values are groupings.  can be a factor.
#	prop	= a data frame or matrix where rownames are the ids, and each column is a property.  NAs will be properly ignored per column.
#			Internalyl, "prop" is split by column, and NAs are removed within each split column, separately.  Then, for each
#			NA-removed split column of "prop" separately, "groups" will be forced to be a subset of whatever is left in that column.
#	outfbase	= if provided, will write the results to <outfbase>.txt
#	plot	= if TRUE, will make a box plot of results.
#	topN	= the topN most significant results will be plotted if plot is true.

	## NA not allowed
	if (any(is.na(names(groups)))) {
		verb("\n\n\n\nERROR!  NA values found in groups!\n")
		show(groups[is.na(names(groups))])
		stop()
	} #no NA

	groups  <-  groups[!is.na(groups)]

	# factor
	if (!is.factor(groups)) { groups  <-  factor(groups) }


	### pop to list of 1-column data frames
	prop  <-  as.data.frame(prop)
	prn  <-  rownames(prop)
	prop  <-  as.list(prop)
	for (propx  in  names(prop)) {
		df  <-  data.frame( value = prop[[propx]] , row.names = prn )
		colnames(df)  <-  propx
		logi.na  <-  is.na(df[[propx]])
		df  <-  df[ !logi.na , ,drop=FALSE]
		prop[[propx]]  <-  df 
	} # propx
	spldf  <-  split( names(groups) , groups , drop = TRUE )


	### test whole poulation
	verb("test population.\n")

	fmdb  <-  list()
	for (propx in names(prop)) {
		### observations must be within this NA-removed prop column,
		subobs  <-  sapply(spldf , FUN = function(z) z[z %in% rownames(prop[[propx]])],simplify=FALSE)
		len.so  <-  sapply( subobs , length )
		subobs  <-  subobs[len.so > 1]
		if (length(subobs) > 0) {
			fmdb[[propx]]  <-  finiteMedianTest( pop = prop[[propx]] , obs = subobs )
		}
	} # propx

	#fmdb  <-  sapply( prop , FUN = function(propx)  finiteMedianTest( pop = propx , obs = sapply(spldf , FUN = function(z) z[z %in% rownames(propx)],simplify=FALSE) ) , simplify = FALSE )
	if (length(fmdb) > 0) {
		fmdb  <-  do.call( rbind , fmdb )
		fmdb$condition1  <-  "population"
		fmdb$condition2  <-  fmdb$sample
	}

	verb("test pairs.\n")
	### pairs
	fm.list  <-  data.frame()

	udirecs  <-  intersect( levels(groups) ,  unique(groups) )
	dpairs  <-  combn(udirecs, 2)
	for (dpx  in  1:ncol(dpairs)) {
		g1  <-  dpairs[1,dpx]
		g2  <-  dpairs[2,dpx]
		verb("\t\t%s  vs  %s\n", g1,g2)

		subde  <-  groups[ groups %in% c(g1,g2) ]
		subpop  <-  sapply( prop , FUN = function(x) x[rownames(x) %in% names(subde) , ,drop=FALSE] , simplify = FALSE )
		len.pop <-  sapply( subpop , nrow )
		subpop  <-  subpop[len.pop > 1]
		if (length(subpop) > 0) {
			obs  <-  names(subde)[subde == g2]
			for (propx  in  names(subpop)) {
				subobs  <-  names(subde)[ subde == g2 ]
				subobs  <-  subobs[ subobs %in% rownames(subpop[[propx]]) ]
				if (length(subobs) > 0) {
					pfmdb  <-  finiteMedianTest( pop = subpop[[propx]] , obs = subobs )
					pfmdb$condition1  <-  g1
					pfmdb$condition2  <-  g2
					fm.list  <-  rbind(fm.list , pfmdb )
				} # obs > 0
			} # propx
		} # subpop > 0
	} # dpx
	all.fmdb  <-  fm.list
	all.fmdb  <-  rbind( all.fmdb , fmdb )
	all.fmdb$p.value  <-  exp(all.fmdb$logp)
	all.fmdb  <-  all.fmdb[ order(all.fmdb$p.value) , ,drop=FALSE]
	all.fmdb$sample  <-  NULL

	### save
	if (!is.null(outfbase)) {
		odb  <-  all.fmdb
		odb$FDR  <-  p.adjust( odb$p.value , method = "BH" )

		fname  <-  sprintf("%s.txt" , outfbase )
		write.table( odb , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
	} # save


	### plot
	topdf  <-  all.fmdb[ all.fmdb$p.value <= 0.05 , ,drop=FALSE]
	if (plot &&  nrow(topdf) > 0) {
		verb("plot.\n")

		topdf  <-  topdf[ order(topdf$p.value) , ,drop=FALSE ]
		topdf  <-  topdf[ 1:min(topN,nrow(topdf)) , , drop=FALSE]

		for (rx  in  1:nrow(topdf)) {
			g1  <-  topdf$condition1[rx]
			g2  <-  topdf$condition2[rx]
			featx  <-  topdf$feature[rx]
			pval  <-  topdf$p.value[rx]	
			verb("\t\t\t\t%s  %s  %s  %g\n", g1 , g2 , featx , pval )

			if (g1 == "population") { 
				subde  <-  groups
				subpop  <-  prop[[featx]]
			} else {
				subde  <-  groups[ groups %in% c(g1,g2) ]
				subpop  <-  prop[[featx]]
				subpop  <-  subpop[ rownames(subpop) %in% names(subde) , ,drop=FALSE]
			} # g1

			# pop df
			mdf  <- data.frame( id = rownames(subpop) , variable = featx , value = subpop[[featx]] )
			colnames(mdf)[3]  <-  featx

			dedf  <-  data.frame( id = names(subde) , group = subde )
			mdf  <-  merge( mdf , dedf , by = "id" , all = TRUE )

			if (g1 == "population") {
				mdf$group[is.na(mdf$group)]  <-  g1
			}

			gp  <-  ggplot( data = mdf , aes_string( x = "group" , y = featx ) )   +
				geom_jitter( width = 0.65 , height = 0 , alpha = 0.05 )   +
				geom_violin( scale = "count" , alpha = 0.3 )   +
				geom_boxplot( width = 0.15 , fill = "white" , alpha = 0.3 , outlier.shape = NA )   +
				theme_bw()   +
				ggtitle(sprintf("%s, p = %g\n%s  vs  %s", featx , pval , g1,g2 ))
			print(gp)
		} # rx
	} # plot


	return(all.fmdb)

} # test_property_enrichment










plot_codon_usage_bias_heatmap  <-  function( fmdb , codon , groups ) {
if (nrow(fmdb) > 0) {

	aaprops  <-  get_amino_acid_properties()
	aamap  <-  get_codon_aminoAcid_map()

	codon  <-  str_replace_all(codon,"U","T")
	fmdb$codon  <-  codon

	### split
	spldb  <-  split( fmdb , groups )

	for (splx  in  names(spldb))  {
		verb("\t%s\n", splx)

		subdb  <-  spldb[[splx]]

		# get amino acid, property, number of genes
		subdb  <-  merge( subdb , aamap , by = "codon" , all = TRUE )
		subdb  <-  merge( subdb , aaprops , by="aminoAcid" , all.x = TRUE )
		subdb  <-  subdb  %>%  group_by(aminoAcid)  %>%  mutate(num.codons = n_distinct(codon))
		subdb  <-  as.data.frame(subdb)
		
		# order by AA property
		subdb  <-  subdb[ order(subdb$AAprop,-subdb$num.codons,subdb$codon) , ,drop=FALSE]
	
		# give row number by amino acid
		subdb$aminoAcid  <-  factor( subdb$aminoAcid , levels = unique(subdb$aminoAcid) )
		subdb  <-  subdb  %>%  group_by(aminoAcid)  %>%  mutate( x = row_number() )
		subdb  <-  as.data.frame(subdb)
	
		# FDR value within group
		subdb$p.value  <-  exp(subdb$logp)
		subdb$FDR  <-  p.adjust(subdb$p.value , method = "BH" )
		subdb$log2FC[subdb$FDR > 0.05]  <-  0
	
		### plot
		textSize  <-  3
		angle  <-  0

	
		gp  <-  ggplot( data = subdb , aes( x = x , y = aminoAcid , fill = log2FC , label = codon ) )   +
			geom_raster()   +
			scale_fill_gradient2( low = muted("red") , mid = "white", high = muted("green") , midpoint = 0 , na.value = "grey54" , limits = c(-1,1) )   +
			geom_text( hjust = "center" , vjust = "center" , size = textSize  , angle = angle )   +
			ggtitle(sprintf("codon usage bias\n%s" , splx))   +
			theme_cowplot()   +
			theme(  axis.title.x=element_blank() , axis.text.x=element_blank() , axis.ticks.x=element_blank()  )
		print(gp)
	}
} # subdb > 0
} # plot_codon_usage_bias_heatmap











#' Annotate dataframe with GO terms.
#'
#' Will append all relevant GO term IDs and defintions (but not full description) to each gene in a dataframe
#'
#' @param df   input data frame to which GO terms will be appended.  must have a column "gene"
#' @param godb  data frame> the GO association database base.  If null, will load from file using "species".
#' @param godefs  definition of each GO term.  loaded from file if null.
#' @param species used to load godb from file if godb is null.
#'
#' @return dataframe with GO ids and defintions.  WIll have multiple lines per input row if multiple GO ids pertain 
#'	to a single gene.
#'
#' export
add_GO_terms_to_gene_list  <-  function( df , godb = NULL , godefs = NULL , species = NULL ) {

	if (is.null(godb)) {
		godb  <-  load_species_GO_db( species = species )
	} # load godb

	if (is.null(godefs)) {
		godefs  <-  load.GO.definitions()
	} # load godefs

	df  <-  merge( godb[, c("gene","GO.id")] , df , by = "gene" , all.y = TRUE )
	df  <-  merge( df , godefs[,c("go_id","Term","Ontology")] , by.x = "GO.id" , by.y = "go_id" , all.x = TRUE )
	df  <-  unique(df)

	return(df)

} # add_GO_terms_to_gene_list






























################################ load
################################ load
################################ load
################################ load
################################ load
################################ load
################################ load
################################ load
################################ load
################################ load


load_raw_counts_and_make_matrix_project  <-  function( project , type = "mrna" , legacy = FALSE ) {

	project  <-  process_project_input( project = project )
        str_project <- project$project[1]

	counts  <-  data.frame()
	for (rx  in  1:nrow(project)) {
		pname  <-  project$project[rx]
		indiv  <-  project$individual[rx]
		runid  <-  project$runid[rx]
		strat  <-  project$strategy[rx]
		spec  <-  project$species[rx]
		verb("\t\treading %s  %s  %s  %s  %s\n", pname , indiv , runid , strat , spec )

		fname  <-  NULL
                id.name  <-  "gene_id"
                idtag  <-  "genes"
                if (legacy) {
                        cdir  <-  sprintf("./out/rsem-mrna/%s/%s" , pname , runid)
                        ctag  <-  sprintf("%s.%s.rsem.mrna" , pname , runid , idtag )
                        fname  <-  sprintf("%s/%s.%s.results.gz" , cdir , ctag , idtag )
		} else {
			cdir  <-  sprintf("./%s_rdna_rn18s/%s/%s" , pname , indiv , runid)
			ctag  <-  sprintf("htseq-annot.%s.%s.%s.%s" , pname , indiv ,  runid , idtag )
			fname  <-  sprintf("%s/%s.results.gz" , cdir , ctag  )
		} # if

		# remove fraction name
		fname <- gsub("wcl\\.|nucleoplasm\\.|chromatin\\.", '.', fname)
		verb('\t\t\t%s\n', fname)
		indf  <-  read.table( file = fname , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = ""  )
		colnames(indf) <- c(id.name, "expected_count")
                # do not remove __no_feature, __ambiguous, __too_low_aQual, __not_aligned, __alignment_not_unique, since # of total reads need to be transfered to DGEList, but the meta tag rows from the DGEList object should be removed before downstream analyses. https://www.biostars.org/p/379690/,https://support.bioconductor.org/p/104972/
		#indf <- indf[!grepl('^__', indf[,1]),]
		indf$runid  <-  runid
		counts  <-  rbind( counts , indf )
	} # rx


	### merge
	verb("merge.\n")

	counts$expected_count  <-  as.numeric(counts$expected_count)
	form  <-  as.formula(paste( id.name , "runid" , sep = "~" ))
	cmat  <-  acast( counts , formula = form , value.var = "expected_count",  fill = 0 )
	stopifnot(ncol(cmat) == nrow(project))

	colnames(cmat)  <-  as.character(colnames(cmat)) # to be safe, in case runids are integers

	stopifnot(all(colnames(cmat) %in% project$runid))
	stopifnot(all(project$runid %in% colnames(cmat)))

	# order properly
	cmat  <-  cmat[, as.character(project$runid) , drop=FALSE]

	stopifnot(identical(colnames(cmat) , project$runid))

	return(cmat)

} # load_rsem_counts_and_make_matrix_project




# begin of addition by H. Kim
load_star_final_log  <-  function( project , type = "mrna" , level = "gene" , legacy = FALSE ) {

        project  <-  process_project_input( project = project )
        # begin of modification by H. Kim
        str_project <- project$project[1]
        #if (level == "gene") {
        f_htseq <- FALSE
        if (level == "protein_spectral_count" || level == "protein_total_peptides") {
          cdir  <-  sprintf("/home/hkim5/riboprof/out/protein_spectral_count/%s" , str_project)
          fname <- sprintf("%s/%s_protein_spectral_count.txt", cdir, str_project)
          cmat <- read.table(fname, sep='\t', header=T, quote='')

          # Group: 100127.1
          # Reference:  sp|P62908|RS3_MOUSE
          # Description
          # 40S ribosomal protein S3 OS=Mus musculus GN=Rps3 PE=1 SV=1
          df1 <- stringr::str_match(cmat[,2],'\\|(.+)\\|(.+)')
          df2 <- stringr::str_match(cmat[,3],'(.+) OS=(.+) GN=(.+) PE=(.+) SV=(.+)')
          rownames(cmat) <- sprintf("%s_%s", df1[,2],df2[,4])
          cols <- colnames(cmat)
          switch(level,
                "protein_spectral_count"={ cmat <- cmat[, grepl("SC",cols)] },
                "protein_total_peptides"={ cmat <- cmat[, grepl("TP",cols)] },
                {}
          )
          return(cmat)
        }

        if ( grepl('^htseq', level) ) f_htseq <- TRUE
        if (level == "gene" || level == "htseq_gene") {
        # end of modification
                id.name  <-  "gene_id"
                idtag  <-  "genes"
        } else {
                id.name  <-  "transcript_id"
                idtag  <-  "isoforms"
        }

        counts  <-  data.frame()
        for (rx  in  1:nrow(project)) {
                pname  <-  project$project[rx]
                indiv  <-  project$individual[rx]
                runid  <-  project$runid[rx]
                strat  <-  project$strategy[rx]
                spec  <-  project$species[rx]
                verb("\t\treading %s  %s  %s  %s  %s\n", pname , indiv , runid , strat , spec )

                fname  <-  NULL
                #if (strat == "RNAseq"  &&  spec == "mg1655") {
                if (spec == "mg1655") {
                        cdir  <-  sprintf("/home/hkim5/riboprof/out/rsem-bacteria-rnaseq-transcriptome/%s/%s/%s" , pname , indiv , runid )
                        ctag  <-  sprintf("rsem.transcriptome.%s.%s.%s" , pname , indiv , runid)
                        fname  <-  sprintf("%s/%s.%s.results.gz" , cdir , ctag , idtag )
                        if (!file.exists(fname)) {
                                verb("\n\nERROR!  file = [%s] does not exist!\n", fname)
                                stop()
                                #cdir  <-  sprintf("/home/hkim5/riboprof/out/rsem-mrna/%s/%s" , pname ,  runid )
                                #ctag  <-  sprintf("rsem.transcriptome.%s.%s" , pname  , runid)
                                #fname  <-  sprintf("%s/%s.%s.results.gz" , cdir , ctag , id.name )
                        } # old
                # begin of addition by H. Kim
                } else if ( f_htseq ) {
                        switch (type,
                          "edger.mrna.fusion"={
                            cdir  <-  sprintf("/home/hkim5/riboprof/out/rnaseq-fusion-gene-detection/%s/%s/%s" , pname , indiv , runid)
                          },
                          {
                            cdir  <-  sprintf("/home/hkim5/riboprof/out/rsem-rnaseq-annotated-genome/%s/%s/%s" , pname , indiv , runid)
                          }
                        ) # switch
                        ctag  <-  sprintf("star-annot.%s.%s.%s" , pname , indiv ,  runid )
                        fname  <-  sprintf("%s/%s.Log.final.out" , cdir , ctag  )
                # end of addition
                } else {
                        if (legacy) {
                                cdir  <-  sprintf("/home/hkim5/riboprof/out/rsem-mrna/%s/%s" , pname , runid)
                                ctag  <-  sprintf("%s.%s.rsem.mrna" , pname , runid , idtag )
                                fname  <-  sprintf("%s/%s.%s.results.gz" , cdir , ctag , idtag )
                        } else {
                                switch (type,
                                  "edger.mrna.fusion"={
                                    cdir  <-  sprintf("/home/hkim5/riboprof/out/rnaseq-fusion-gene-detection/%s/%s/%s" , pname , indiv , runid)
                                  },
                                  {
                                    cdir  <-  sprintf("/home/hkim5/riboprof/out/rsem-rnaseq-annotated-genome/%s/%s/%s" , pname , indiv , runid)
                                  }
                                ) # switch
                                ctag  <-  sprintf("rsemRnaseq-annot.%s.%s.%s" , pname , indiv ,  runid )
				# star_for_rsem.Log.final.out
                                fname  <-  sprintf("%s/%s.star_for_rsem.Log.final.out" , cdir , ctag  )
                        } # legacy
                } # species

		fname <- gsub("wcl\\.|nucleoplasm\\.|chromatin\\.", '.', fname)
                indf  <-  read.table( file = fname , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = "" , fill = T )
		colnames(indf) <- c(id.name, "expected_count")
		indf <- indf[!grepl("READS:", indf[,1]),]
		indf[,1] <- trimws(indf[,1])
		indf[,1] <- gsub("[ ]*\\|[ ]*", '', indf[,1])
                indf$runid  <-  runid
                counts  <-  rbind( counts , indf )
        } # rx


        ### merge
        verb("merge.\n")

        form  <-  as.formula(paste( id.name , "runid" , sep = "~" ))
        df  <-  reshape2::dcast( counts , formula = form , value.var = "expected_count",  fill = 0 )
        rownames(df) <- df[,1]
        df <- df[, 2:ncol(df)]
        df <- df[indf[,1], ] # reorder rows
        stopifnot(ncol(df) == nrow(project))

        colnames(df)  <-  as.character(colnames(df)) # to be safe, in case runids are integers

        stopifnot(all(colnames(df) %in% project$runid))
        stopifnot(all(project$runid %in% colnames(df)))

        # order properly
        df  <-  df[, as.character(project$runid) , drop=FALSE]

        stopifnot(identical(colnames(df) , project$runid))

        return(df)

} 
# end of addition









































#' Plot correlation among runids within each condition
#'
#' This function plots the correlation of data (e.g. gene expression log2cpm)
#' between every pair of runids pertaining to the same condition, for each condition.
#'
#' @param  project	data frame of project information
#' @param  mat		numeric matrix. columns = runid names, rows = element names (e.g. gene names).
#'
#' @return	nothing. prints plot to current device
#'
#' @export
plot_correlation_per_runid  <-  function( project , mat , verbose = FALSE ) {

	if (verbose) { verb("\n\nplot_correlation_per_runid\n") }

	#mat  <-  load_limma_results( project = project , type = type , level = level , rundate_appendix = rundate_appendix , result = result , verbose = verbose )

	for (condx  in  unique(project$condition)) {
		if (verbose) { verb("\t\t%s\n" , condx) }

		rids  <-  project$runid[project$condition == condx]
		rids  <-  make.names(rids)
		subdf  <-  as.data.frame(mat[, rids ,drop=FALSE])
		
		if (ncol(subdf) < 2) { next }
	
		gp  <-  ggscatmat( data = subdf , alpha = 0.5 )
		gp  <-  gp  +  
			ggtitle(sprintf("%s", condx ))
		print(gp)
	
	} # condx
} # plot_correlation_per_runid
	














