suppressMessages(library(motifStack))
source_lines <- readLines("threshold_test.R")
# reuse helper definitions (everything before the main loop)
eval(parse(text = source_lines[1:grep("^res <- list", source_lines)[1] - 1]))
inv <- read.csv(file.path(motif_dir, "motif_inventory_human_only.csv"), check.names = FALSE, fileEncoding = "UTF-8-BOM")
g <- ifelse(inv$category == "PPAR:RXR heterodimer", "BOTH", toupper(inv$name))
inv <- inv[order(match(g, c("BOTH","PPARA","PPARD","PPARG","RXRA","RXRB","RXRG")), match(sub(" .*","",inv$database), c("JASPAR","HOCOMOCO","CIS-BP")), inv$motif_id), ]
set.seed(1)
orders <- list(sheet = names(pfms), display = inv$motif_id, reversed = rev(names(pfms)), random1 = sample(names(pfms)), random2 = sample(names(pfms)), random3 = sample(names(pfms)))
for (o in names(orders)) {
  lst <- c(pfms[orders[[o]]], list(ppre)); names(lst) <- sapply(lst, function(m) m@name)
  al <- DNAmotifAlignment(lst, threshold = 0.5, minimalConsensus = 6)
  r <- evaluate(al); cat(sprintf("order=%-9s %s\n", o, paste(names(r), r, collapse = " | ")), flush = TRUE)
}
