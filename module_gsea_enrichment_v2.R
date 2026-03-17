#===============================================================================
# Module Gene Extraction and Enrichment
#
# Functions:
#   - extract_module_genes()
#   - enrich_modules_fast()
#
#' @author Charles Opara
#===============================================================================
#' Extract genes (SYMBOLs/Ensembl IDs) for each module
#'
#' @param mod_data data.frame or data.table with module assignments and gene info.
#' @param trait_class_keep character, which trait class to keep (default "RNA").
#' @param module_col character, name of column containing module IDs.
#' @param trait_col character, name of column containing trait class.
#' @param symbol_col character, name of column containing gene symbols.
#' @param ensembl_col character, name of column containing Ensembl IDs.
#' @param kME_col character, name of column containing module eigengene connectivity (kME).
#' @param min_kME numeric or NULL, minimum kME threshold; if NULL, all genes kept.
#'
#' @return data.table with one row per module containing:
#'   - symbols (list-column of gene SYMBOLs)
#'   - ensembl_ids (list-column of Ensembl IDs)
#'   - n_symbols, n_ensembl (counts)
#'   - median_kME
#'   - hub gene info (hub_symbol, hub_ensembl, hub_kME)
#' @export
#-------------------------------------------------------------------------------
#' Run GO enrichment analysis per module
#'
#' @param module_tables data.table or *named list* of tables produced by
#'   extract_module_genes(). Each table must have `module` and `symbols` columns.
#' @param ont character, GO ontology to use ("BP","MF","CC"). Default: "BP".
#' @param p_cut numeric, p-value cutoff. Default: 0.05.
#' @param q_cut numeric, q-value cutoff (FDR). Default: 0.1.
#' @param min_set integer, minimum number of genes (ENTREZ IDs) per module to test.
#' @param minGSSize integer, minimum gene set size for enrichment. Default: 5.
#' @param maxGSSize integer, maximum gene set size. Default: 5000.
#' @param use_universe logical, whether to restrict enrichment to the module’s
#'   own background universe. Default: FALSE.
#' @param readable logical, whether to map ENTREZ IDs back to gene symbols.
#'   Default: TRUE.
#' @param workers integer, number of parallel workers to use.
#'
#' @return A list with three elements:
#'   - diagnostics: summary of module sizes per dataset
#'   - enrichment_results: list of tables with raw enrichments (per module)
#'   - enriched_terms: flattened data.table with enriched terms across modules
#'
#' @export


# Required libraries
# data.table, clusterProfiler, org.Mm.eg.db, GO.db, future.apply, parallelly
suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(GO.db)
  library(future.apply)
  library(parallelly)
})

#===============================================================================
# 1) Extract module gene SYMBOLS into a table
#    Input format (example):
#    mod_data: columns = module | trait class | Ensembl id.Symbol | Ensembl id | Symbol | kME_module
#===============================================================================
extract_module_genes <- function(mod_data,
                                 trait_class_keep = "RNA",
                                 module_col = "module",
                                 trait_col = "trait class",
                                 symbol_col = "Symbol",
                                 ensembl_col = "Ensembl id",
                                 kME_col = "kME_module",
                                 min_kME = NULL) {
  
  dt <- as.data.table(mod_data)
  
  # checks
  need <- c(module_col, trait_col, symbol_col, ensembl_col, kME_col)
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("Missing required column(s): ", paste(miss, collapse = ", "))
  
  # filter to requested trait class (e.g., RNA)
  dt <- dt[get(trait_col) == trait_class_keep]
  
  # optional kME threshold
  if (!is.null(min_kME)) {
    if (!is.numeric(min_kME) || length(min_kME) != 1L) stop("'min_kME' must be a single number.")
    dt <- dt[!is.na(get(kME_col)) & get(kME_col) >= min_kME]
  }
  
  # normalize types
  dt[, (module_col) := as.character(get(module_col))]
  dt[, (symbol_col) := as.character(get(symbol_col))]
  dt[, (ensembl_col) := as.character(get(ensembl_col))]
  
  # drop empty/NA symbols
  dt <- dt[!is.na(get(symbol_col)) & nzchar(get(symbol_col))]
  
  # compute hubs (max kME per module)
  hubs <- dt[ , .SD[which.max(get(kME_col))][1L],
              by = get(module_col),
              .SDcols = c(symbol_col, ensembl_col, kME_col)]
  setnames(hubs, old = c("get", symbol_col, ensembl_col, kME_col),
           new = c("module", "hub_symbol", "hub_ensembl", "hub_kME"))
  
  # assemble "hotspot-like" table:
  # one row per module; list of symbols; list of ensembl ids; counts
  out <- dt[ , .(
    symbols = list(unique(get(symbol_col))),
    ensembl_ids = list(unique(get(ensembl_col))),
    n_symbols = uniqueN(get(symbol_col)),
    n_ensembl = uniqueN(get(ensembl_col)),
    median_kME = median(get(kME_col), na.rm = TRUE)
  ), by = get(module_col)]
  setnames(out, "get", "module")
  
  # join hub info
  out <- hubs[out, on = "module"]
  
  # for full analogy to hotspot_tables, you already have the crucial pieces:
  # - a list column 'symbols'
  # - a count column 'n_symbols'
  # You can add any other module metadata here if needed.
  setorder(out, module)
  out[]
}

#===============================================================================
# 2) GO enrichment per module (fast), modeled after enrich_hotspots_fast()
#    Accepts either:
#      - a single data.table from extract_module_genes(), or
#      - a *named list* of such tables (e.g., multiple tissues/datasets)
#    Required columns in each table: 'module', 'symbols' (list of SYMBOLs)
#===============================================================================
enrich_modules_fast <- function(module_tables,
                                ont = "BP",
                                p_cut = 0.05, q_cut = 0.10,
                                min_set = 5, # min ENTREZ per module to test
                                minGSSize = 5, maxGSSize = 5000,
                                use_universe = FALSE,
                                readable = TRUE,
                                workers = max(1, parallelly::availableCores() - 1)) {
  
  # helper: ensure list input with names
  as_named_list <- function(x) {
    if (is.data.table(x) || is.data.frame(x)) {
      return(list(modules = as.data.table(x)))
    } else if (is.list(x)) {
      if (is.null(names(x))) stop("If passing a list of tables, it must be *named*.")
      return(lapply(x, as.data.table))
    } else {
      stop("module_tables must be a data.table/data.frame or a *named* list of them.")
    }
  }
  mod_list <- as_named_list(module_tables)
  
  # Ensure required cols in each table
  for (nm in names(mod_list)) {
    dt <- mod_list[[nm]]
    req <- c("module", "symbols")
    miss <- setdiff(req, names(dt))
    if (length(miss)) stop(sprintf("Table '%s' is missing columns: %s", nm, paste(miss, collapse=", ")))
    # normalize
    if (!is.list(dt$symbols)) stop(sprintf("Table '%s': 'symbols' must be a list column.", nm))
    dt[, module := as.character(module)]
    mod_list[[nm]] <- dt
  }
  
  ##-----------------------------
  ## Build GO maps once (for ont)
  ##-----------------------------
  keys_entrez <- keys(org.Mm.eg.db, keytype = "ENTREZID")
  ann <- as.data.table(AnnotationDbi::select(
    org.Mm.eg.db,
    keys = keys_entrez,
    columns = c("GO","ONTOLOGY"),
    keytype = "ENTREZID"
  ))
  ann <- ann[!is.na(GO) & ONTOLOGY == ont, .(GO, ENTREZID)]
  TERM2GENE <- unique(data.frame(
    term = ann$GO,
    gene = ann$ENTREZID,
    stringsAsFactors = FALSE
  ))
  
  t2n <- AnnotationDbi::select(
    GO.db,
    keys = unique(ann$GO),
    columns = "TERM",
    keytype = "GOID"
  )
  TERM2NAME <- unique(data.frame(
    term = t2n$GOID,
    name = t2n$TERM,
    stringsAsFactors = FALSE
  ))
  
  # Map all symbols -> ENTREZ once across ALL modules
  all_syms <- unique(unlist(
    lapply(mod_list, \(dt) unlist(dt$symbols, use.names = FALSE)),
    use.names = FALSE
  ))
  all_syms <- unique(all_syms[!is.na(all_syms) & nzchar(all_syms)])
  sym_map <- suppressMessages(
    bitr(all_syms,
         fromType = "SYMBOL",
         toType = "ENTREZID",
         OrgDb = org.Mm.eg.db,
         drop = TRUE)
  )
  sym2entrez <- setNames(sym_map$ENTREZID, sym_map$SYMBOL)
  
  map_syms <- \(syms) {
    x <- unique(unname(sym2entrez[as.character(syms)]))
    x[!is.na(x)]
  }
  universe_for <- \(dt) {
    if (!isTRUE(use_universe)) return(NULL)
    u <- unique(unlist(dt$symbols, use.names = FALSE))
    u <- unique(unname(sym2entrez[as.character(u)]))
    u[!is.na(u)]
  }
  
  # Generic filter + flatten helpers
  filter_trim <- function(dt, obj_col = "go") {
    if (!nrow(dt)) return(dt)
    objs <- dt[[obj_col]]
    objs <- lapply(objs, function(ego) {
      if (is.null(ego) || is.null(ego@result) || !nrow(ego@result)) return(NULL)
      res <- ego@result
      keep <- (!is.na(res$Count)) & res$Count >= 5 &
        (!is.na(res$qvalue)) & res$qvalue <= 0.05
      res2 <- res[keep, , drop = FALSE]
      if (!nrow(res2)) return(NULL)
      ego@result <- res2
      ego
    })
    dt[[obj_col]] <- objs
    keep_mod <- vapply(dt[[obj_col]], function(x)
      !is.null(x) && !is.null(x@result) && nrow(x@result) > 0L, logical(1))
    dt[keep_mod]
  }
  
  flatten <- function(dt, obj_col = "go") {
    # one row per enriched term per module
    if (!nrow(dt)) return(data.table())
    objs <- dt[[obj_col]]
    rows <- lapply(seq_len(nrow(dt)), function(i) {
      ego <- objs[[i]]
      if (is.null(ego)) return(NULL)
      out <- as.data.table(ego@result)
      out[, `:=`(
        module = dt$module[i],
        n_syms = lengths(dt$symbols)[i],
        n_entrez = dt$n_entrez[i]
      )]
      out
    })
    rbindlist(rows, fill = TRUE)
  }
  
  diag_tbl <- function(name, dt) {
    if (is.null(dt) || !nrow(dt)) {
      return(data.table(dataset=name, n_modules=0,
                        med_syms=NA_real_, med_entrez=NA_real_))
    }
    ns <- lengths(dt$symbols)
    ne <- vapply(dt$symbols, \(s) length(map_syms(s)), integer(1))
    data.table(dataset=name,
               n_modules = nrow(dt),
               med_syms = if (any(ns>0)) median(ns[ns>0]) else NA_real_,
               med_entrez = if (any(ne>0)) median(ne[ne>0]) else NA_real_)
  }
  
  # Parallel processing
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)
  
  # Per-dataset work
  results <- lapply(names(mod_list), function(nm) {
    base_dt <- as.data.table(mod_list[[nm]])
    
    if (!nrow(base_dt)) {
      return(list(
        name = nm,
        enr = data.table(),
        enr_kegg = data.table(),
        diag = diag_tbl(nm, base_dt),
        terms = data.table(),
        terms_kegg = data.table()
      ))
    }
    
    uni <- universe_for(base_dt)
    entrez_list <- lapply(base_dt$symbols, map_syms)
    idx <- which(vapply(entrez_list, length, 0L) >= min_set)
    
    # Lists to hold GO and KEGG objects
    go_list <- vector("list", nrow(base_dt))
    kegg_list <- vector("list", nrow(base_dt))  ## KEGG
    
    if (length(idx)) {
      # GO enrichment
      go_list[idx] <- future_lapply(idx, function(i) {
        gene <- entrez_list[[i]]
        ego <- tryCatch(
          enricher(
            gene = gene,
            universe = uni,
            TERM2GENE = TERM2GENE,
            TERM2NAME = TERM2NAME,
            pAdjustMethod = "BH",
            pvalueCutoff = p_cut,
            qvalueCutoff = q_cut,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize
          ),
          error = function(e) NULL
        )
        if (!is.null(ego) && readable) {
          ego <- tryCatch(
            setReadable(ego, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"),
            error = function(e) ego
          )
        }
        ego
      }, future.seed = TRUE)
      
      # KEGG
      kegg_list[idx] <- future_lapply(idx, function(i) {
        gene <- entrez_list[[i]]
        kk <- tryCatch(
          enrichKEGG(
            gene = gene,
            organism = "mmu",        # mouse
            universe = uni,
            pAdjustMethod = "BH",
            pvalueCutoff = p_cut,
            qvalueCutoff = q_cut,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize
          ),
          error = function(e) NULL
        )
        if (!is.null(kk) && readable) {
          kk <- tryCatch(
            setReadable(kk, OrgDb = org.Mm.eg.db, keyType = "ENTREZID"),
            error = function(e) kk
          )
        }
        kk
      }, future.seed = TRUE)
    }
    
    # Build GO and KEGG module tables separately
    n_entrez <- vapply(entrez_list, length, 0L)
    
    dt_go <- base_dt[, .(module, symbols)]
    dt_go[, `:=`(go = go_list, n_entrez = n_entrez)]
    dt_go_f <- filter_trim(dt_go, obj_col = "go")
    
    dt_kegg <- base_dt[, .(module, symbols)]
    dt_kegg[, `:=`(kegg = kegg_list, n_entrez = n_entrez)]
    dt_kegg_f <- filter_trim(dt_kegg, obj_col = "kegg")
    
    list(
      name = nm,
      enr = dt_go,
      enr_kegg = dt_kegg,
      diag = diag_tbl(nm, base_dt),
      terms = flatten(dt_go_f, obj_col = "go"),
      terms_kegg = flatten(dt_kegg_f, obj_col = "kegg")
    )
  })
  
  # Assemble outputs
  list(
    diagnostics = rbindlist(lapply(results, `[[`, "diag"), fill = TRUE),
    enrichment_results = setNames(lapply(results, `[[`, "enr"),
                                      lapply(results, `[[`, "name")),
    enriched_terms = setNames(lapply(results, `[[`, "terms"),
                                      lapply(results, `[[`, "name")),
    # KEGG outputs
    enrichment_results_kegg = setNames(lapply(results, `[[`, "enr_kegg"),
                                       lapply(results, `[[`, "name")),
    enriched_terms_kegg = setNames(lapply(results, `[[`, "terms_kegg"),
                                       lapply(results, `[[`, "name"))
  )
}



#===============================================================================
# Example usage
#===============================================================================
# Read in data
mod_data <- fread("/Users/copara/Downloads/Module_membership_md20_beta12_4enrichment.csv")

# Build module tables (this extracts genes lists from each module into a table.. 
# .... it requires a row label 'RNA' indicating that the trait is a gene)
module_tbl <- extract_module_genes(mod_data, trait_class_keep = "RNA", min_kME = 0.5)

# Enrichment per module (this performs the enrichment on each mod...*set you desired thresholds)
mod_enrichment <- enrich_modules_fast(module_tbl, ont = "BP", p_cut = 0.05, q_cut = 0.1)

# If you have multiple datasets (e.g., different tissues), pass a *named list*:
# module_tables <- list(liver = module_tbl, adipose = module_tbl2, muscle = module_tbl3)
# mod_enrichment <- enrich_modules_fast(module_tables)

# GO output
fwrite(mod_enrichment$enriched_terms$modules, file = "./modules_enrichment_GO.csv", sep = ",")

# KEGG output
fwrite(mod_enrichment$enriched_terms_kegg$modules, file = "./modules_enrichment_KEGG.csv", sep = ",")

