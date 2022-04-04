install.packages("phangorn")
install.packages("readr")
install.packages("ape")
install.packages("TreeDist")
library(phangorn)
library(readr)
library(ape)
library(TreeDist)

# load viridian simulations trees
simulated_art_phylogeny <- readr::read_file("~/Documents/GitHub/viridian_simulations_workflow/simulated_phylogenies/ART_assemblies/uncondensed-final-tree.nh")
simulated_art_phylogeny <- ape::read.tree(text=simulated_art_phylogeny)
artic_art_phylogeny <- readr::read_file("~/Documents/GitHub/viridian_simulations_workflow/artic_phylogenies/ART_assemblies/uncondensed-final-tree.nh")
artic_art_phylogeny <- ape::read.tree(text=artic_art_phylogeny)
viridian_art_phylogeny <- readr::read_file("~/Documents/GitHub/viridian_simulations_workflow/viridian_phylogenies/ART_assemblies/uncondensed-final-tree.nh")
viridian_art_phylogeny <- ape::read.tree(text=viridian_art_phylogeny)
# calculate tree distances
artic_simulated <- treedist(simulated_art_phylogeny, artic_art_phylogeny, check.labels = TRUE)
viridian_simulated <- treedist(simulated_art_phylogeny, viridian_art_phylogeny, check.labels = TRUE)
artic_viridian <- treedist(artic_art_phylogeny, viridian_art_phylogeny, check.labels = TRUE)
# calculate Robinson-Foulds distance
artic_simulated_RF <- phangorn::RF.dist(simulated_art_phylogeny, artic_art_phylogeny, check.labels=TRUE, rooted=FALSE)
viridian_simulated_RF <- phangorn::RF.dist(simulated_art_phylogeny, viridian_art_phylogeny, check.labels=TRUE, rooted=FALSE)
artic_viridian_RF <- phangorn::RF.dist(artic_art_phylogeny, viridian_art_phylogeny, check.labels=TRUE, rooted=FALSE)
# calculate Kendall–Colijn tree distance
artic_simulated_KC <- TreeDist::KendallColijn(simulated_art_phylogeny, artic_art_phylogeny)
viridian_simulated_KC <- TreeDist::KendallColijn(simulated_art_phylogeny, viridian_art_phylogeny)
artic_viridian_KC <- TreeDist::KendallColijn(artic_art_phylogeny, viridian_art_phylogeny)


# load SA analysis trees
gisaid_illumina_phylogeny <- readr::read_file("~/Documents/GitHub/vridian_SA_analysis/gisaid_phylogenies/illumina/uncondensed-final-tree.nh")
gisaid_illumina_phylogeny <- ape::read.tree(text=gisaid_illumina_phylogeny)
viridian_illumina_phylogeny <- readr::read_file("~/Documents/GitHub/vridian_SA_analysis/viridian_phylogenies/illumina/uncondensed-final-tree.nh")
viridian_illumina_phylogeny <- ape::read.tree(text=viridian_illumina_phylogeny)
gisaid_nanopore_phylogeny <- readr::read_file("~/Documents/GitHub/vridian_SA_analysis/gisaid_phylogenies/nanopore/uncondensed-final-tree.nh")
gisaid_nanopore_phylogeny <- ape::read.tree(text=gisaid_nanopore_phylogeny)
viridian_nanopore_phylogeny <- readr::read_file("~/Documents/GitHub/vridian_SA_analysis/viridian_phylogenies/nanopore/uncondensed-final-tree.nh")
viridian_nanopore_phylogeny <- ape::read.tree(text=viridian_nanopore_phylogeny)
# calculate tree distances
viridian_gisaid_illumina <- treedist(gisaid_illumina_phylogeny, viridian_illumina_phylogeny, check.labels = TRUE)
viridian_gisaid_nanopore <- treedist(gisaid_nanopore_phylogeny, viridian_nanopore_phylogeny, check.labels = TRUE)
# calculate Robinson-Foulds distance
viridian_gisaid_illumina_RF <- phangorn::RF.dist(gisaid_illumina_phylogeny, viridian_illumina_phylogeny, check.labels=TRUE, rooted=FALSE)
viridian_gisaid_nanopore_RF <- phangorn::RF.dist(gisaid_nanopore_phylogeny, viridian_nanopore_phylogeny, check.labels=TRUE, rooted=FALSE)
# calculate Kendall–Colijn tree distance
gisaid_viridian_illumina_KC <- TreeDist::KendallColijn(gisaid_illumina_phylogeny, viridian_illumina_phylogeny)
gisaid_viridian_nanopore_KC <- TreeDist::KendallColijn(gisaid_nanopore_phylogeny, viridian_nanopore_phylogeny)

