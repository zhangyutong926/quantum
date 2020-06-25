import math
import dimod as dm
import dwavebinarycsp as dmcsp

N = 6
M = 3
d = [6, 12]
e = [6, 6, 6]
D = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
E = range(1, 5) + range(1, 6) + ([0] * 22)

def avarname(i, j, k):
    return "a^" + str(i) + "_{" + str(j) + "," + str(k) + "}"
def bvarname(i, j, k):
    return "b^" + str(i) + "_{" + str(j) + "," + str(k) + "}"
def cvarname(i, n, m):
    return "c^" + str(i) + "_{" + str(n) + "," + str(m) + "}"
def gtvarname(i):
    return "gt^" + str(i)
def otarname(i):
    return "ot^" + str(i)
def ptarname(i):
    return "pt^" + str(i)
def yvarname(i, l):
    return "y^" + str(i) + "_{" + str(l) + "}"
varlist = [avarname(i, j, k) for k in range(N) for j in len(d) for i in range(ceiling(math.log2(d[j])))] +
          [bvarname(i, j, k) for k in range(N) for j in len(e) for i in range(ceiling(math.log2(e[j])))] +
          [bvarname(i, n, m) for n in range(N) for m in range(M) for j in len(d) + len(e)] + # FIXME Q(n) -> Q(log_2(n))
          [yvarname(i, l) for l in len(D) for i in range(ceiling(math.log2(e[j])))] +
          [gtvarname(i) for i in ceiling(math.log2(range(N)))] + [otvarname(i) for i in ceiling(math.log2(range(N)))] + [ptvarname(i) for i in ceiling(math.log2(range(N)))]
def todec(*argv):
    return int("".join([str(x) for x in argv]), 2)

bqm = dm.BinaryQuadraticModel({}, {}, 0.0, dimod.SPIN)
csp = dmcsp.ConstraintSatisfactionProblem([], varlist, dwavebinarycsp.BINARY)

for n in range(N):
    for m in range(M):
        var = [bvarname(i, n, m) for j in len(d) + len(e)]
        csp.add_constraint(
            dm.generators.constraints.combinations(var, M), var
        )

for l in range(len(E)):
    csp.add_constraint(
        lambda *argv: E[l] <= todec(argv),
        [yvarname(i, l) for i in range(ceiling(math.log2(e[j])))]
    )

# FIXME BEGIN boilerplate, hardcoding -> decouplation
for l in range(1, 11):
    csp.add_constraint(
        lambda *argv: todec(argv[ceiling(math.log2(range(N))):]) <= l 
                      if todec(argv[:ceiling(math.log2(range(N)))]) is 0
                      else todec(argv[ceiling(math.log2(range(N))):]) > l
        [gtvarname(i) for i in range(N)] + \
        [yvarname(i, l) for i in range(ceiling(math.log2(e[0])))]
    )
    csp.add_constraint(
        lambda *argv: todec(argv[:ceiling(math.log2(e[0]))]) > 0 or todec(argv[ceiling(math.log2(e[0])):]) == 0,
        [yvarname(i, l - 1) for i in range(ceiling(math.log2(e[0])))] + \
        [yvarname(i, l) for i in range(ceiling(math.log2(e[0])))]
    )
for l in range(12, 22):
    csp.add_constraint(
        lambda *argv: todec(argv[ceiling(math.log2(range(N))):]) <= l - 11 
                      if todec(argv[:ceiling(math.log2(range(N)))]) is 0
                      else todec(argv[ceiling(math.log2(range(N))):]) > l - 11
        [gtvarname(i) for i in range(N)] + \
        [yvarname(i, l) for i in range(ceiling(math.log2(e[1])))]
    )
    csp.add_constraint(
        lambda *argv: todec(argv[:ceiling(math.log2(e[1]))]) > 0 or todec(argv[ceiling(math.log2(e[1])):]) == 0,
        [yvarname(i, l - 1) for i in range(ceiling(math.log2(e[1])))] + \
        [yvarname(i, l) for i in range(ceiling(math.log2(e[1])))]
    )
for l in range(23, 33):
    csp.add_constraint(
        lambda *argv: todec(argv[ceiling(math.log2(range(N))):]) <= l - 22
                      if todec(argv[:ceiling(math.log2(range(N)))]) is 0
                      else todec(argv[ceiling(math.log2(range(N))):]) > l - 22
        [gtvarname(i) for i in range(N)] + \
        [yvarname(i, l) for i in range(ceiling(math.log2(e[2])))]
    )
    csp.add_constraint(
        lambda *argv: todec(argv[:ceiling(math.log2(e[2]))]) > 0 or todec(argv[ceiling(math.log2(e[2])):]) == 0,
        [yvarname(i, l - 1) for i in range(ceiling(math.log2(e[2])))] + \
        [yvarname(i, l) for i in range(ceiling(math.log2(e[2])))]
    )
    csp.add_constraint(
        lambda *argv: todec(argv[:ceiling(math.log2(e[2]))]) < todec(argv[ceiling(math.log2(e[2])):]),
        [yvarname(i, l - 1) for i in range(ceiling(math.log2(e[2])))] + \
        [yvarname(i, l) for i in range(ceiling(math.log2(e[2])))]
    )
# FIXME END


# TODO CSP to BQM

for l in range(len(D)):
    for i in range(ceiling(math.log2(e[l]))):
        bqm.add_variable(yvarname(i, l), D[l])

# TODO dwave solve
