
class Employee:

    empCount = 0

    def __init__(self, name, salary):
        self.name = name
        self.salary = salary
        Employee.empCount += 1

    def displayCount(self):
        print "Total employee %d" % Employee.empCount

    def displayEmployee(self):
        print "Name: ", self.name, ", Salary: ", self.salary


class Point:

    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def __del__(self):
        class_name = self.__class__.__name__
        print class_name, "BALETED!"


pt1 = Point()
pt2 = pt1
pt3 = pt1

print id(pt1), id(pt2), id(pt3)
del pt1
del pt2
del pt3

str1 = "this is string example....wow!!!";
str2 = "exam";

print str1.find(str2);
print str1.find(str2, 10);
print str1.find(str2, 40);

# "This would create first object of Employee class"
# emp1 = Employee("Zara", 2000)
# "This would create second object of Employee class"
# emp2 = Employee("Manni", 5000)
# emp1.displayEmployee()
# emp2.displayEmployee()
# print "Total Employee %d" % Employee.empCount
#
# print "Employee.__doc__:", Employee.__doc__
# print "Employee.__name__:", Employee.__name__
# print "Employee.__module__:", Employee.__module__
# print "Employee.__bases__:", Employee.__bases__
# print "Employee.__dict__:", Employee.__dict__



def de_bruijn(k, n):
    import operator
    """
    De Bruijn sequence for alphabet k
    and subsequences of length n.
    """
    try:
        # let's see if k can be cast to an integer;
        # if so, make our alphabet a list
        _ = int(k)
        alphabet = list(range(k))
        alphabet = map(str, alphabet)

    except (ValueError, TypeError):
        alphabet = k
        k = len(k)

    a = [0] * k * n
    sequence = []
    def db(t, p):
        if t > n:
            if n % p == 0:
                for j in range(1, p + 1):
                    sequence.append(a[j])
        else:
            a[t] = a[t - p]
            db(t + 1, p)
            for j in range(a[t - p] + 1, k):
                a[t] = j
                db(t + 1, t)
    db(1, 1)
    return "".join(map(alphabet.__getitem__, sequence))

print(de_bruijn(2, 3))
print(de_bruijn("abcd", 2))



# Get the edge elements.
DBG_edge_elmts = set()
for kmer in k_mers:
	DBG_edge_elmts.add(kmer)
	DBG_edge_elmts.add(RevComp(kmer))

# Create the edges.
k = len(k_mers[0])
edge = lambda elmmt: '('+elmt[0:k-1]+', '+elmt[1:k]+')'
DBG_edges = [edge(elmt) for elmt in DBG_edge_elmts]

# Write and save the adjacency list.
print '\n'.join(DBG_edges)
with open('output/057_DBRU.txt', 'w') as output_file:
	output_file.write('\n'.join(DBG_edges))