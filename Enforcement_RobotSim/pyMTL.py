import mtl
import sys
import json
#from sys import argv

s = set()
for i in sys.argv[1:]:
    s.add(i)
formula = sys.argv[1]


a, b = mtl.parse('a'), mtl.parse('b')


# Assumes piece wise constant interpolation.


dataType = {
    'a': [(0, True), (1, False), (3, False)],
    'b': [(0, False), (0.2, True), (3, False)]
}

js = json.dumps(dataType)
file = open('/test.txt', 'w')
file.write(js)
file.close()



file = open('/test.txt', 'r')
js = file.read()
data = json.loads(js)
#print(data)
file.close()


phi = mtl.parse(formula)
print(phi(data, quantitative=False, time= 4))
# output: True

#phi = mtl.parse('F(a | b)')
#print(phi(data, quantitative=False))
# output: True

# Note, quantitative parameter defaults to False

# Evaluate at t=3.
#print(phi(data, quantitative=False, time=3))
# output: False

# Compute sliding satisifaction.
#print(phi(data, time=None))
# output: [(0, True), (0.2, True), (4, False)]

# Evaluate with discrete time
#phi = mtl.parse('X b')
#print(phi(data, dt=0.2))
# output: True
print ("the string variable is:  G (Speed_x(robot_i) < 3 & Robot_dis(robot_i) > 20)")
