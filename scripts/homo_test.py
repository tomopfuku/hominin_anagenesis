import mandos
import sys

if len(sys.argv) != 5:
    print "python "+sys.argv[0]+" <tree> <alignment> <partitions file> <range data file>"
    sys.exit(0)

infile = sys.argv[1]
tree = mandos.tree_utils2.read_tree(infile)

seqs = mandos.tree_utils2.read_phylip_file(sys.argv[2])
parts = sys.argv[3]
sitels = mandos.tree_utils2.read_partition_file(parts)

ranges = mandos.tree_utils2.read_strat(sys.argv[4])
mandos.tree_utils2.match_strat(tree,ranges)
mandos.tree_utils2.init_heights_strat(tree)

#mandos.tree_utils2.make_ancestor(tree,"Homo_heidelbergensis")
taxa = ['Homo_sapiens', 'Homo_neanderthalensis','Homo_heidelbergensis','Homo_antecessor']
subtree = mandos.tree_utils2.prune_SA(tree,taxa)



print "calculating AIC of fully bifurcating arrangement"
print subtree.get_newick_repr(True)
morphopt = mandos.tree_likelihood_calculator.calc_mk_like(sitels,subtree,seqs,True)
morpholike = morphopt[0]
nparams = morphopt[1]
opt = mandos.stratoML.optim_lambda_heights(subtree,ranges)
nparams += opt[2]
print subtree.get_newick_repr(True)
strat_like = -opt[1][1]
print strat_like,morpholike
combined = strat_like + morpholike
print nparams
print (2*(nparams))-(2*combined)
print mandos.tree_utils2.aic(combined,nparams)

"""
print "making H. heidelbergensis a SA and calculating AIC"
subtree = mandos.tree_utils2.make_ancestor(subtree,"Homo_heidelbergensis")
print subtree.get_newick_repr(True)
morpholike = mandos.tree_likelihood_calculator.calc_mk_like(sitels,subtree,seqs,True,True)
opt = mandos.stratoML.optim_lambda_heights(subtree,ranges)
print subtree.get_newick_repr(True)
strat_like = -opt[1][1]
print strat_like,morpholike
combined = strat_like + morpholike
print mandos.tree_utils2.tree_AIC(subtree,combined,6)
"""


print "making H. antecessor a SA and calculating AIC"
subtree = mandos.tree_utils2.make_ancestor(subtree,"Homo_heidelbergensis")
print subtree.get_newick_repr(True)

morphopt = mandos.tree_likelihood_calculator.calc_mk_like(sitels,subtree,seqs,True,True)
morpholike = morphopt[0]
nparams = morphopt[1]
opt = mandos.stratoML.optim_lambda_heights(subtree,ranges,True)
nparams += opt[2]
print subtree.get_newick_repr(True)
strat_like = -opt[1][1]
print strat_like,morpholike
combined = strat_like + morpholike
print (2*(nparams))-(2*combined)
print mandos.tree_utils2.aic(combined,nparams)


