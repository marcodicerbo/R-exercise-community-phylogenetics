### a bit of community phylogenetics 

# https://cran.r-project.org/web/packages/fishtree/vignettes/community-analysis.html

# Exercise

# Vedremo se le specie di pesci con pinne raggiate associate alla barriera corallina sono raggruppate o eccessivamente disperse negli 
# oceani Atlantico, Pacifico e Indiano

install.packages("fishtree")
install.packages("rfishbase")
library(fishtree)
library(rfishbase)
library(geiger)

reef_species <- rfishbase::species(fields = c("SpecCode", "DemersPelag"))  # tagliamo la tabella prendendo solo le colonne con il codice 
# della sp. e demerspelag
reef_species <- reef_species[reef_species$DemersPelag=="reef-associated",] # prendiamo solo le specie di reef

# prendiamo le specie endemiche o native degli oceani Atlantico, Pacifico e Indiano
eco <- rfishbase::ecosystem(species_list = reef_species$SpecCode)
valid_idx <- eco$Status %in% c("native","endemic")& eco$EcosystemName %in% c("Atlantic Ocean", "Pacific Ocean", "Indian Ocean")
eco <- eco[valid_idx,c("Species","EcosystemName")]

phy <- fishtree_phylogeny(species = eco$Species) # Recupera la filogenesi delle sole specie autoctone della barriera corallina in tutti e 
# tre gli oceani

# E' necessario ripulire i dati prima di effettuare l'analisi con picante. Innanzitutto, convertire il dataframe delle specie per sito in  
# una matrice di presenza-assenza. Utilizzeremo base::table per questo e utilizzeremo unclass per convertire la tabella in una matrice.

sample_matrix <- unclass(table(eco))
dimnames(sample_matrix)$Species <- gsub(" ","_",dimnames(sample_matrix)$Species,fixed = T)

# Successivamente, utilizzeremo geiger::name.check per garantire che le etichette dei tips della filogenesi e le righe della matrice di 
# dati corrispondano tra loro

nc <- geiger::name.check(phy,sample_matrix)
sample_matrix <- sample_matrix[!rownames(sample_matrix) %in% nc$data_not_tree,]

# Infine, genereremo la matrice cofenetica basata sulla filogenesi e trasporremo la matrice presenza-assenza poiché a picante piace che le 
# sue colonne siano specie e le sue righe siano siti

cophen <- cophenetic(phy) # Una matrice cofenetica è una rappresentazione numerica delle distanze filogenetiche tra le coppie
# di unità tassonomiche (ad esempio, specie, individui, sequenze di DNA) in un albero filogenetico. 
# Definizione: Per un albero filogenetico, la matrice cofenetica è una matrice quadrata in cui ciascun elemento 
# (i, j) rappresenta la distanza cofenetica tra il taxon i e il taxon j. La distanza cofenetica è definita come 
# la distanza dall'antenato comune più recente (MRCA - Most Recent Common Ancestor) di i e j alla radice dell'albero.

sample_matrix <- t(sample_matrix)

# Eseguiremo picante::ses.mpd e picante::ses.mntd con solo 100 iterazioni, per velocizzare l'analisi. Per un'analisi reale aumenteremmo questo 
# valore a 1000 ed eventualmente testeremo altri null models se i set di dati contengono, ad esempio, informazioni sull'abbondanza.

picante::ses.mpd(sample_matrix,cophen,null.model = "taxa.labels", runs = 99)
picante::ses.mntd(sample_matrix,cophen,null.model = "taxa.labels",runs = 99)

# Gli oceani Atlantico e Indiano sono sovradispersi utilizzando la metrica MPD e tutti e tre gli oceani sono raggruppati sotto la metrica MNTD. 
# Si ritiene che l'MPD sia più sensibile ai modelli più vicini alla radice dell'albero, mentre si ritiene che l'MNTD rifletta più da vicino i 
# modelli verso le punte della filogenesi.
# Possiamo confermare visivamente questi modelli eseguendo il codice seguente, che traccerà la filogenesi e aggiungerà punti colorati 
# (rosso, verde e blu) per indicare se una punta è associata a uno specifico bacino oceanico.

plot(phy,show.tip.label = F,no.margin = T)
obj <- get("last_plot.phylo", .PlotPhyloEnv)
matr <- t(sample_matrix)[phy$tip.label,]
xx <- obj$xx[1:obj$Ntip]
yy <- obj$yy[1:obj$Ntip]
cols <- c("#1b9e77", "#d95f02", "#7570b3")
for (ii in 1:ncol(matr)) {
  present_idx <- matr[,ii]==1
  points(xx[present_idx]+ii,yy[present_idx],col=cols[ii],cex=0.1)
}