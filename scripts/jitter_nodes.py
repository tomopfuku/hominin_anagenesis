import mandos
import sys

if len(sys.argv) != 2:
    print "python "+sys.argv[0]+" <tree> "
    sys.exit(0)

infile = sys.argv[1]
tree = mandos.tree_utils2.read_tree(infile)

for i in tree.iternodes():
    if i == tree:
        continue
    if i.length < 0.02:
        i.length  = 0.02

print tree.get_newick_repr(True)

