import numpy as np

"""
Verbatim translation of routine found in tauola.
"""

def calcR(p1,p2,p3,e):
    return np.sqrt(e**2 - p1**2 - p2**2 - p3**2)/2

def spherd(R):
    r1, r2 = np.random.rand(2)
    costh = 2*r1 -1
    sinth = np.sqrt(1 - costh**2)
    return [R*sinth*np.cos(2*np.pi*r2), R*sinth*np.sin(2*np.pi*r2), R*costh, R]

def boost(v, p):
    amv = np.sqrt(abs(v[3]**2 - v[0]**2 -v[1]**2 -v[2]**2))
    t = (p[0]*v[0] + p[1]*v[1] + p[2]*v[2] + p[3]*v[3])/amv
    wsp = (t + p[3])/(v[3] + amv)
    return [p[0] + wsp*v[0], p[1] + wsp*v[1], p[2] + wsp*v[2], t]

def decay(p1, p2, p3, e):
    R = calcR(p1, p2, p3, e)
    X = spherd(R)
    Y = [-X[0], -X[1], -X[2], R]

    phot1 = boost( [p1,p2,p3,e], X)
    phot2 = boost( [p1,p2,p3,e], Y)
    return phot1, phot2

if __name__ == "__main__":
    import sys

    assert(sys.argv[1] != sys.argv[2])

    g = open(sys.argv[2], "w")

    with open(sys.argv[1]) as f:
        comments=f.readline()
        g.write(comments)
        for line in f:
            l = line.split(",")
            pid = l[1]
            if pid == "111":
                p1 = float(l[3])
                p2 = float(l[4])
                p3 = float(l[5])
                e  = float(l[2])
                PH1, PH2 = decay(p1,p2,p3,e)
                l1 = [x for x in l]
                l2 = [x for x in l]
                l1[1] = "22"
                l1[2] = str(PH1[-1])
                l1[3] = str(PH1[0])
                l1[4] = str(PH1[1])
                l1[5] = str(PH1[2])
                l2[1] = "22"
                l2[2] = str(PH2[-1])
                l2[3] = str(PH2[0])
                l2[4] = str(PH2[1])
                l2[5] = str(PH2[2])
                g.write(",".join(l1))
                g.write(",".join(l2))
            else:
                g.write(line)
    g.close()
