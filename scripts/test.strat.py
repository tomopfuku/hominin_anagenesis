import sys
import mandos

if len(sys.argv) != 3:
    print "usage: "+sys.argv[0]+ " <newick file> <range data file>"
    sys.exit(0)

tree = mandos.tree_utils2.read_tree(sys.argv[1])
ranges = mandos.tree_utils2.read_strat(sys.argv[2])
mandos.tree_utils2.match_strat(tree,ranges)
ranges = None

mandos.tree_utils2.init_heights_strat(tree)
#mandos.tree_utils2.make_ancestor(tree,"Homo_heidelbergensis")
#print tree.get_newick_repr(True)


opt = mandos.stratoML.optim_lambda_heights(tree,ranges)

#print opt
t1= tree.get_newick_repr(True)
t2= opt[0].get_newick_repr(True)
print opt[1][1]
print t1
