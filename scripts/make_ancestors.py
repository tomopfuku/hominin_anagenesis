import mandos
import sys

if len(sys.argv) != 3:
    print "python "+sys.argv[0]+" <tree> <comma-separated tips>"
    sys.exit(0)

infile = sys.argv[1]
tree = mandos.tree_utils2.read_tree(infile)

tree = mandos.tree_utils2.make_ancestor(tree,sys.argv[2])

print tree.get_newick_repr()

