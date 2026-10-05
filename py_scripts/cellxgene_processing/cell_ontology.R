library(ontologyIndex)

# Scarica e carica la Cell Ontology (.obo format) dal repository ufficiale
# (Potrebbe impiegare qualche secondo per il download)
cl_onto <- readRDS("/projects/shared/cellxgene_split_by_dataset/cell_ontology.rds")

example_dataset <- "/projects/shared/cellxgene_split_by_dataset/intestine/2adb1f8a-a6b1-4909-8ee8-484814e2d4bf.h5ad"

sce <- anndataR::read_h5ad(example_dataset)

cl_onto
unique(sce$cell_type_ontology_term_id)
unique(sce$cell_type)

# ID generici/radici da ignorare per evitare che le cellule diventino semplicemente "cell"
generic_roots <- c("CL:0000000", "CL:0000003", "CL:0000001") # cell, native cell, etc.

auto_reduce_cell_type <- function(term_id, onto) {
  # 1. FORZA A TESTO: risolve il problema dei numeri
  term_id <- as.character(term_id) 
  
  if (is.na(term_id) || !term_id %in% names(onto$name)) {
    return("Unknown")
  }
  
  # Trova i genitori diretti
  parents <- onto$parents[[term_id]]
  
  # Filtra: mantieni solo genitori CL validi, non generici e non obsoleti
  valid_parents <- parents[
    grepl("^CL:", parents) & 
      !parents %in% generic_roots & 
      !grepl("obsolete", onto$name[parents], ignore.case = TRUE)
  ]
  
  # Se esiste un genitore più generale ma valido, usa quello; altrimenti tieni l'originale
  if (length(valid_parents) > 0) {
    target_id <- valid_parents[1]
  } else {
    target_id <- term_id
  }
  
  return(onto$name[target_id])
}

# 2. Applica la funzione convertendo esplicitamente in character
sce$auto_broad_cell_type <- sapply(
  as.character(sce$cell_type_ontology_term_id), 
  auto_reduce_cell_type, 
  onto = cl_onto
)

# 3. Controlla il risultato pulito
unique(sce$auto_broad_cell_type, useNA = "always")

# no non funziona, usiamo ontoProc

library(ontoProc)

cl <- getOnto("cellOnto", "2023")
opi <- graph_from_graphnel(onto_plot2(cl_onto))
