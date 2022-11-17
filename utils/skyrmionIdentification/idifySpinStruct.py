#!/usr/bin/env python3
# full name: identifySpinStructure.py


import csv
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import shelve
import sys


# definition of a point:
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        return

    def __eq__(self, otherPoint):
        if self.x == otherPoint.x and self.y == otherPoint.y:
            return True
        return False

    def __str__(self):
        return "(" + str(self.x) + ", " + str(self.y) + ")"

class Segment:
    def __init__(self, point1=0, point2=0, spin1=0, spin2=0, orientation=0):
        self.p1 = point1
        self.p2 = point2
        self.s1 = spin1
        self.s2 = spin2
        self.orientation = orientation # 0 for horizontal segment, 1 for vertical one
        return

    def __eq__(self, otherSegment):
        if self.p1 == otherSegment.p1 and self.p2 == otherSegment.p2:
            return True
        if self.p1 == otherSegment.p2 and self.p2 == otherSegment.p1:
            return True
        return False

    def __hash__(self):
        global superCellSize
        return int(self.p1.x + self.p2.x)*superCellSize + int(self.p1.y + self.p2.y) + self.orientation

    def __str__(self):
        return str(self.p1) + " " + str(self.p2)

    def isAdjacent(self, otherSegment):
        if self.p1 == otherSegment.p1 or self.p1 == otherSegment.p2:
            return True
        if self.p2 == otherSegment.p1 or self.p2 == otherSegment.p2:
            return True
        return False

    def getInPlaneMedSpin(self):
        return [(self.s1[0]+self.s2[0])/2, (self.s1[1]+self.s2[1])/2]

    def swapPoints(self):
        tmp = self.p1
        self.p1 = self.p2
        self.p2 = tmp

        tmp = self.s1
        self.s1 = self.s2
        self.s2 = tmp
        return

    def plot(self, ax, color='k'):
        ax.plot([self.p1.x, self.p2.x], [self.p1.y, self.p2.y], color=color)

    def povrayPlot(self, fileName):
        f_handle = open(fileName, "a")
        f_handle.write("cylinder{ < " + str(-self.p1.x) + ", "
                                      + str(self.p1.y) + ", -0.25 >, < "
                                      + str(-self.p2.x) + ", "
                                      + str(self.p2.y) + ", -0.25 >, 0.1 }\n")
        f_handle.write("Spin(" + str((self.p1.x+self.p2.x)/2) + ", "
                               + str((self.p1.y+self.p2.y)/2) + ", "
                               + str(0.0) + ", "
                               + str((self.s1[0]+self.s2[0])/2)+ ", "
                               + str((self.s1[1]+self.s2[1])/2) + ", "
                               + str((self.s1[2]+self.s2[2])/2) + ")\n")
        f_handle.close()

class Curve:
    def __init__(self, initialSegment=None):
        if initialSegment is None:
            self.charge = None
            return
        self.segments = []
        self.segments.append(initialSegment)
        self.startSeg = initialSegment
        self.endSeg = initialSegment
        self.charge = None
        self.innerRegion = None
        return

    def __str__(self):
        retString = "Curve{\n"
        for seg in self.segments:
            retString += "\t" + str(seg) + "\n"
        retString += "}"
        return retString

    def lightCopy(self):
        copy = Curve()
        copy.charge = self.charge
        return copy

    def isEquivalent(self, other):
        if self.getCharge() == other.getCharge():
            return True
        return False

    def size(self):
        return len(self.segments)

    def addSeg(self, seg):
        if self.size() == 1:
            if self.startSeg.p1 == seg.p1 or self.startSeg.p1 == seg.p2:
                self.startSeg.swapPoints()
        if seg.p1 != self.endSeg.p2:
            seg.swapPoints()
        self.segments.append(seg)
        self.endSeg = seg
        return

    def isClosed(self):
        if self.startSeg.isAdjacent(self.endSeg) and len(self.segments) > 2:
            return True
        return False

    def isConnected(self, seg):
        if self.endSeg.isAdjacent(seg):
            return True
        return False

    def getRightMostSegment(self):
        retSeg = None
        for seg in self.segments:
            if seg.orientation != 1: # if not a vertical segment
                continue
            if retSeg is None:
                retSeg = seg
            if seg.p1.x > retSeg.p1.x: # p1.x and p2.x are equivalent since vertical segments
                retSeg = seg
        return retSeg

    # return  1 if counterclockwise
    # return -1 if clock-wise
    def getCirculationSign(self):
        seg = self.getRightMostSegment()
        if seg.p2.y > seg.p1.y:
            return 1
        return -1

    def getChirality(self):
        global spinLat
        inner = self.getInnerRegion()
        elem = iter(inner).__next__()
        rightSeg = self.getRightMostSegment()
        x = (rightSeg.p1.x + rightSeg.p2.x)/2 - 0.5
        y = (rightSeg.p1.y + rightSeg.p2.y)/2
        inner_magnetization = np.sign(spinLat[int(x), int(y)][2])
        outer_magnetization = - inner_magnetization
        chirality = (inner_magnetization - outer_magnetization) / 2
        return chirality

    def computeCharge(self):
        theta = []
        s1 = self.segments[-1].getInPlaneMedSpin()
        phi1 = np.arctan2(s1[1], s1[0])
        for i in range(len(self.segments)):
            s2 = self.segments[i].getInPlaneMedSpin()
            phi2 = np.arctan2(s2[1], s2[0])
            theta_i = phi2 - phi1
            if theta_i > np.pi:
                theta_i = theta_i - 2*np.pi
            if theta_i < -np.pi:
                theta_i = theta_i + 2*np.pi
            theta.append(theta_i)
            s1 = s2
            phi1 = phi2

        phaseTot = 0
        for theta_i in theta:
            phaseTot += theta_i
        sign = self.getCirculationSign()
        chirality = self.getChirality()
        # print("sign: " + str(sign))
        # print("chirality: " + str(chirality))
        # print("phaseTot/(2*pi): " + str(phaseTot/(2*np.pi)))
        self.charge = round(sign * chirality * phaseTot / (2*np.pi))
        return

    def getCharge(self):
        if self.charge == None:
            self.computeCharge()
        return self.charge

    def setCharge(self, charge):
        self.charge = charge
        return

    def printInPlaneMedSpin(self):
        for seg in self.segments:
            print(seg.getInPlaneMedSpin())
        return

    def plot(self, ax):
        if round(self.getCharge()) == -1:
            color = 'r'
        elif round(self.getCharge()) == 1:
            color = 'b'
        else:
            color = 'k'
        for seg in self.segments:
            seg.plot(ax, color)

    def povrayPlot(self, fileName):
        for seg in self.segments:
            seg.povrayPlot(fileName)
        return

    def addNeighbor(self, queue, point):
        x = point[0]
        y = point[1]
        queue.append((x-1, y))
        queue.append((x+1, y))
        queue.append((x, y-1))
        queue.append((x, y+1))
        return

    def determineInnerRegion(self):
        # building of the contour
        # contour: set of all points forming the contour
        contour = set()
        circulationSign = self.getCirculationSign()
        for seg in self.segments:
            normal_x = circulationSign * (seg.p2.y - seg.p1.y)
            normal_y = - circulationSign * (seg.p2.x - seg.p1.x)
            x = (seg.p1.x + seg.p2.x)/2 + normal_x * 0.5
            y = (seg.p1.y + seg.p2.y)/2 + normal_y * 0.5
            contour.add((x, y))

        # init inner Region
        innerRegion = set()
        rightSeg = self.getRightMostSegment()
        x = (rightSeg.p1.x + rightSeg.p2.x)/2 - 0.5
        y = (rightSeg.p1.y + rightSeg.p2.y)/2
        startPoint = (x, y)
        innerRegion.add(startPoint)
        stack = []
        self.addNeighbor(stack, startPoint)

        # loop
        while len(stack) != 0:
            p = stack.pop()
            if p not in innerRegion and p not in contour:
                innerRegion.add(p)
                self.addNeighbor(stack, p)
        self.innerRegion = innerRegion
        return

    def getInnerRegion(self):
        if self.innerRegion == None:
            self.determineInnerRegion()
        return self.innerRegion

    def plotInnerRegion(self, ax):
        # converting the set in list
        X = list()
        Y = list()
        for p in self.getInnerRegion():
            X.append(p[0])
            Y.append(p[1])
        # plot
        if round(self.getCharge()) == -1:
            color = 'r'
        elif round(self.getCharge()) == 1:
            color = 'b'
        else:
            color = 'k'
        ax.scatter(X, Y, color=color, marker='s')
        return

class TreeLeaf:
    def __init__(self, curve=None):
        self.curve = curve
        self.childs = []
        return

    def __eq__(self, other):
        if self.curve.getCharge() != other.curve.getCharge():
            return False
        if len(self.childs) != len(other.childs):
            return False
        return False

    def addChild(self, leaf):
        self.childs.append(leaf)
        return

    def isEquivalent(self, other):
        if self.curve is not None:
            if not self.curve.isEquivalent(other.curve):
                return False
        if len(self.childs) != len(other.childs):
            return False

        mask = [0]*len(self.childs)
        for c1 in self.childs:
            i = 0
            for c2 in other.childs:
                if c1.isEquivalent(c2) and mask[i] == 0:
                    mask[i] = 1
                    break
                i += 1
            if i == len(self.childs):
                return False
        return True

    def print(self, degree=0):
        if degree != 0 and self.curve is not None:
            for i in range(degree-1):
                print(" ", end='')
            print("Tree leaf of degree " + str(degree) + ", charge is " + str(self.curve.getCharge()) + " with " + str(len(self.childs)) + " child tree leaves")
        for c in self.childs:
            c.print(degree=degree+1)
        return

    def lightCopy(self):
        if self.curve is None:
            copy = TreeLeaf()
        else:
            copy = TreeLeaf(self.curve.lightCopy())
        for c in self.childs:
            copy.addChild(c.lightCopy())
        return copy

    def plot(self, ax):
        if self.curve is None:
            for child in self.childs:
                child.plot(ax)
            return

        seg = self.curve.getRightMostSegment()
        x = (seg.p1.x + seg.p2.x) / 2
        y = (seg.p1.y + seg.p2.y) / 2
        for child in self.childs:
            # print("Hello there")
            segChild = child.curve.getRightMostSegment()
            xChild = (segChild.p1.x + segChild.p2.x) / 2
            yChild = (segChild.p1.y + segChild.p2.y) / 2
            ax.annotate('', xy=(xChild, yChild), xytext=(x, y),
                        arrowprops=dict(facecolor='black', shrink=0.05))
            child.plot(ax)
        return

class TreeStruc:
    # file database located at script home directory
    simuRoot = os.path.dirname(__file__)
    dbFile = simuRoot + '/spinStruct'

    def __init__(self, curves=None):
        if curves is None or len(curves) == 0:
            self.root = None
            self.degree = None
            return
        # print(curves)
        # print("")
        containedBy = []
        i = 0
        for c1 in curves:
            j = 0
            containedBy.append([])
            for c2 in curves:
                if c1.segments == c2.segments:
                    # print("hi: " + str(i) + ", " + str(j))
                    continue
                if c1.getInnerRegion() < c2.getInnerRegion():
                    # print("hello: " + str(i) + ", " + str(j))
                    containedBy[i].append(c2)
                    # print(containedBy)
                j += 1
            i += 1

        self.degree = max([len(val) for val in containedBy]) + 1 + 1
        # double +1 to take into account: the deepest node itself and the root node

        # init root of tree
        self.root = TreeLeaf()
        level = []
        for i in range(len(curves)):
            c = curves[i]
            if len(containedBy[i]) == 0:
                leaf = TreeLeaf(c)
                self.root.addChild(leaf)
                level.append(leaf)

        # populating tree
        nextLevel = []
        for d in range(1, self.degree):
            # print("level: " + str(level))
            for i in range(len(curves)):
                c1 = curves[i]
                if len(containedBy[i]) != d:
                    continue
                for j in range(len(containedBy[i])):
                    c2 = containedBy[i][j]
                    for ancestor in level:
                        # print("Hoho")
                        if ancestor.curve.segments == c2.segments:
                            # print("Hello!")
                            leaf = TreeLeaf(c1)
                            ancestor.addChild(leaf)
                            nextLevel.append(leaf)
                            break
            level = nextLevel.copy()
            nextLevel.clear()
        return

    def __eq__(self, other):
        return False

    def print(self):
        # print("Tree of degree: " + str(self.degree))
        if self.root is None:
            print("Empty tree")
            return
        self.root.print(degree=0)
        return

    def lightCopy(self):
        copy = TreeStruc()
        copy.degree = self.degree
        if self.root is not None:
            copy.root = self.root.lightCopy()
        return copy

    def plot(self, ax):
        self.root.plot(ax)
        return

    def isEquivalent(self, other):
        if other.root is None and self.root is None:
            return True
        elif self.root is None or other.root is None:
            return False
        if self.root.isEquivalent(other.root):
            return True
        return False

    def compareToDatabase(self):
        database = shelve.open(self.dbFile, flag='c', protocol=5)
        i = 0
        for name, tree in database.items():
            if self.isEquivalent(tree):
                database.close()
                return name
            i += 1

        name = str(i)
        database[name] = self.lightCopy()
        database.close()
        return name

    @classmethod
    def printDatabase(cls):
        database = shelve.open(cls.dbFile, flag='c', protocol=5)
        for name, tree in database.items():
            print("Tree name: " + str(name))
            tree.print()
            print("\n")
        database.close()
        return


def getSpinLat(fileName):
    f_handle = open(fileName, "r")
    reader = csv.reader(f_handle, delimiter=' ', quoting=csv.QUOTE_NONE,
                        skipinitialspace=True)
    superCellSize = 0
    for _ in reader:
        superCellSize += 1
    superCellSize = int(np.sqrt(superCellSize))
    f_handle.seek(0)
    spinLat = np.empty((superCellSize, superCellSize, 3))
    for i in range(superCellSize):
        for j in range(superCellSize):
            spinLat[j, i, :] = np.asarray(reader.__next__(), dtype=np.float32)

    f_handle.close()
    return spinLat

def createSegments(spinLat):
    superCellSize = len(spinLat[:,0,0])
    signLat = np.heaviside(spinLat[:,:,2], 0)
    # print(signLat)
    # fig = plt.figure(figsize=(12,9), tight_layout=True)
    # ax = fig.subplots(1, 1)
    # im = ax.imshow(signLat)
    # plt.show()

    segmentSet = set()

    # sweep from left to right
    orientation = 1 # vertical segment
    for i in range(superCellSize-1):
        for j in range(superCellSize):
            if signLat[i, j] != signLat[i+1, j]:
                p1 = Point(i+0.5, j-0.5)
                p2 = Point(i+0.5, j+0.5)
                seg = Segment(p1, p2, spinLat[i, j, :], spinLat[i+1, j, :], orientation)
                segmentSet.add(seg)

    # sweep from down to top
    orientation = 0 # horizontal segment
    for i in range(superCellSize):
        for j in range(superCellSize-1):
            if signLat[i, j] != signLat[i, j+1]:
                p1 = Point(i-0.5, j+0.5)
                p2 = Point(i+0.5, j+0.5)
                seg = Segment(p1, p2, spinLat[i, j, :], spinLat[i, j+1, :], orientation)
                segmentSet.add(seg)

    return segmentSet

def constructCurves(segmentSet, scripting):
    curves = []
    i = 0
    seg = Segment()
    deleted = set()
    while len(segmentSet) != 0:
        curves.append(Curve(segmentSet.pop()))
        while not curves[i].isClosed():
            segFound = False
            for seg in segmentSet:
                if seg in deleted:
                    continue
                if curves[i].isConnected(seg):
                    curves[i].addSeg(seg)
                    deleted.add(seg)
                    segFound = True
                    break # we break out of 'for seg in segmentSet'
            if segFound is False:
                break # we break out of 'while not curves[i].isClosed()'
        segmentSet -= deleted
        deleted.clear()
        if segFound is True:
            i += 1
        else:
            del curves[i]

    for c in curves:
        c.computeCharge()
    curves = list(filter(lambda c: c.getCharge() != 0, curves))

    return curves

if __name__ == "__main__":

    fileName = 'magnetic_end.dat'
    scripting = False
    verbose = False
    printDatabase = False
    for i in range(1, len(sys.argv)):
        if sys.argv[i] == "--scripting" or sys.argv[i] == "-s":
            scripting = True
        if sys.argv[i] == "--verbose" or sys.argv[i] == "-v":
            verbose = True
        if sys.argv[i] == "--print-database" or sys.argv[i] == "-pdb":
            printDatabase = True

    if printDatabase == True:
        TreeStruc.printDatabase()
        exit()

    spinLat = getSpinLat(fileName)
    superCellSize = len(spinLat[:,0,0])
    segmentSet = createSegments(spinLat)
    curves = constructCurves(segmentSet, scripting)
    for c in curves:
        c.determineInnerRegion()
    # tree-shaped struct
    tree = TreeStruc(curves)
    name = tree.compareToDatabase()
    if not scripting:
        print("The structure is " + name)

    # printing results to screen
    if not scripting and verbose:
        i = 0
        for c in curves:
            print("Curve: " + str(i))
            print("Curve length: " + str(len(c.segments)))
            print("charge: " + str(c.getCharge()))
            print("")
            i += 1

    if not scripting and verbose:
        tree.print()

    # plotting
    fig = plt.figure(constrained_layout=True, figsize=[6.4, 6.4])
    ax = fig.subplots(1, 1)
    ax.tick_params(direction='in')
    ax.tick_params(top=True)
    ax.tick_params(right=True)
    # ax.set_aspect('equal', adjustable='box')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # ax.axis('equal')
    ax.axis('square')
    # ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([0, superCellSize-1])
    ax.set_ylim([0, superCellSize-1])
    for c in curves:
        if c.getCharge() != 0:
            c.plot(ax)
            # c.povrayPlot("Visualization/idify_curves.dat")
            # c.plotInnerRegion(ax)
    # curves[0].plotInnerRegion(ax)
    # tree.plot(ax)
    plt.savefig("plot_" + name + ".pdf")
    if not scripting and verbose:
        plt.show()

    if scripting:
        print(name)
