import mandos
import sys

if len(sys.argv) != 4:
    print "python "+sys.argv[0]+" <tree> <alignment> <partitions file>"
    sys.exit(0)

infile = sys.argv[1]
tree = mandos.tree_utils2.read_tree(infile)

seqs = mandos.tree_utils2.read_phylip_file(sys.argv[2])
parts = sys.argv[3]
sitels = mandos.tree_utils2.read_partition_file(parts)

#mandos.tree_utils2.make_ancestor(tree,"B,D")

print "optimizing branch lengths on: "+ tree.get_newick_repr()
morphopt = mandos.tree_likelihood_calculator.calc_mk_like(sitels,tree,seqs,True,False)
morpholike = morphopt[0]
nparams = morphopt[1]
print "Morphological log-likelihood: "+str(morpholike)
print "AIC: "+str(mandos.tree_utils2.aic(morpholike,nparams))

print tree.get_newick_repr(True)
