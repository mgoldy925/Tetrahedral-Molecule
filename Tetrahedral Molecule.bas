REM Tetrahedral Molecule Simulation
REM Maxwell Goldberg 4/18/21
REM Inspired by Goodie

DIM SHARED xEye AS INTEGER
DIM SHARED yEye AS INTEGER
DIM SHARED zEye AS INTEGER
DIM SHARED yplane AS INTEGER
CONST MAXINT = 2147483647

DefineProj

SCREEN 12
WINDOW (-10, -10)-(10, 10)

n = 4 ' number of hydrogens
DIM SHARED h(n + 1, 3) ' hydrogens, carbon is last particle
DIM SHARED q(n) ' assuming all positive charges on hydrogens
DIM SHARED m(n) ' assuming all have mass of 1kg

dt = 0.001
bl = 2 ' eq bond length
C = 10000 ' k in Coulomb's Law
k = 1 ' spring constant
rc = 1 ' radius of hydrogens
b = 0.1 ' friction

InitializeAtoms 4, n

t = 0
DIM r(n) ' holds radii from carbon
DIM F(n, 3) ' array for net force of each hydrogen
DIM v(n, 3) ' velocities
'DIM oh(n, 3) ' position from previous iteration to erase
DIM order(n + 1, 2) ' order y values for plotting correctly, holds y values and index
DIM d(n, n) ' Array for distances, will be nxn to tell which particles distance is between
DIM index(2) ' Hold two particles w/ shortest distance

press$ = "q"
check = 0
GOTO drawing ' I wanna draw the atoms at the beginning without copy and pasting code
afterBeg:

DO

    press$ = INKEY$

    FOR i = 0 TO n - 1

        FOR k = 0 TO 2 ' adding in friction
            F(i, k) = -1 * SGN(v(i, k)) * b * (ABS(v(i, k))) ^ 2 ' need to add negative sign for friction to work, based on sign of velocity
        NEXT k

        FOR j = 0 TO n ' calculate force acting on ith particle

            IF i = j THEN _CONTINUE

            dx = h(i, 0) - h(j, 0) ' doing i - j gives repulsive component aka negative component
            dy = h(i, 1) - h(j, 1)
            dz = h(i, 2) - h(j, 2)
            r = SQR((dx) ^ 2 + (dy) ^ 2 + (dz) ^ 2)
            r(i) = r ' This will put the distance between ith particle and nth particle, carbon, in array, this is risky, be careful

            IF j = n THEN ' carbon is 5th particle aka index 4
                F = -1 * n * k * (r - bl)
            ELSE
                F = C * q(i) * q(j) / r ^ 2 ' repulsion force
            END IF

            F(i, 0) = F(i, 0) + getComp(dx, dy, dz, r, F, "X")
            F(i, 1) = F(i, 1) + getComp(dx, dy, dz, r, F, "Y")
            F(i, 2) = F(i, 2) + getComp(dx, dy, dz, r, F, "Z")
        NEXT j

        FOR d = 0 TO 2

            a = F(i, d) / m(i)
            dv = a * dt
            v(i, d) = v(i, d) + dv
            dd = v(i, d) * dt
            'oh(i, d) = h(i, d)
            h(i, d) = h(i, d) + dd

        NEXT d

    NEXT i

    drawing:
    ' have to draw particles in back before front
    ' do selection sort

    FOR i = 0 TO n ' go to n bc including carbon
        order(i, 0) = i ' index
        order(i, 1) = h(i, 1) ' value
    NEXT i

    FOR i = 0 TO (n - 1) ' n-1 bc including carbon

        k = i
        min = MAXINT
        index = -1

        FOR j = i TO n

            IF order(j, 1) < min THEN
                k = j
                index = order(j, 0)
                min = order(j, 1)
            END IF

        NEXT j

        order(k, 0) = order(i, 0)
        order(i, 0) = index

        order(k, 1) = order(i, 1)
        order(i, 1) = min

    NEXT i

    ' Find distances between hydrogens
    min = MAXINT

    FOR i = 0 TO n - 1
        FOR j = i + 1 TO n - 1
            d(i, j) = SQR((h(i, 0) - h(j, 0)) ^ 2 + (h(i, 1) - h(j, 1)) ^ 2 + (h(i, 2) - h(j, 2)) ^ 2)
            IF d(i, j) < min THEN
                min = d(i, j)
                index(0) = i
                index(1) = j
            END IF
        NEXT j
    NEXT i

    xcp = Projection(h(n, 0), h(n, 1), h(n, 2), "X") ' These are projected coords for carbon
    yc = h(n, 1) '                                     This is the y coord
    zcp = Projection(h(n, 0), h(n, 1), h(n, 2), "Z") ' They will be used for drawing bonds in the next loop
    size = 2 ' How many times bigger carbon should be

    IF check = 1 THEN CLS
    FOR a = n TO 0 STEP -1 ' remeber, bigger y is further away, so index from n to 0

        i = order(a, 0)

        IF i = n THEN ' give carbon bigger radius and different color
            r = size * rc
            CO = 1
        ELSE
            r = rc
            CO = 12
        END IF

        ' This part is kinda repetitive from before, might make it into a function
        dx = h(i, 0) - h(n, 0)
        dy = h(i, 1) - h(n, 1)
        dz = h(i, 2) - h(n, 2)

        'oxp = Projection(oh(i, 0), oh(i, 1), oh(i, 2), "X") ' old coordinates to paint black over
        'ozp = Projection(oh(i, 0), oh(i, 1), oh(i, 2), "Z")
        xp = Projection(h(i, 0), h(i, 1), h(i, 2), "X")
        zp = Projection(h(i, 0), h(i, 1), h(i, 2), "Z")
        rp = RadiusProjection(h(i, 1), h(i, 2), zp, r)

        IF i <> n THEN ' drawing bonds between surfaces of atoms
            csx = getComp(dx, dy, dz, r(i), r * size, "X")
            csy = getComp(dx, dy, dz, r(i), r * size, "Y")
            csz = getComp(dx, dy, dz, r(i), r * size, "Z")
            hsx = getComp(dx, dy, dz, r(i), r(i) - r, "X")
            hsy = getComp(dx, dy, dz, r(i), r(i) - r, "Y")
            hsz = getComp(dx, dy, dz, r(i), r(i) - r, "Z")

            cspx = Projection(csx, csy, csz, "X")
            cspz = Projection(csx, csy, csz, "Z")
            hspx = Projection(hsx, hsy, hsz, "X")
            hspz = Projection(hsx, hsy, hsz, "Z")

            LINE (cspx, cspz)-(hspx, hspz), 15

        END IF

        'PAINT (oxp, ozp), 0  ' commented out because clearing screen instead of painting black circle over last circle
        CIRCLE (xp, zp), rp, CO
        PAINT (xp, zp), CO

    NEXT a


    IF check = 1 THEN ' If you'd like to change theta to not track the smallest angle
        'h1 = index(0) ' Change this to 0
        'h2 = index(1) ' Change this to 1
        h1 = 0
        h2 = 1
        theta = ARCCOS((r(h1) ^ 2 + r(h2) ^ 2 - d(h1, h2) ^ 2) / (2 * r(h1) * r(h2))) * 180 / _PI

        PRINT "Bond Angle:  "; Rounder(theta, 2) ' NOTE: This angle is the smallest angle in the bond, so the angle will jump and won't hover around true bond angle but instead approach it
        PRINT "Bond Length:  "; Rounder(r(h1), 2) ', Rounder(r(1), 2), Rounder(r(2), 2), Rounder(r(3), 2)
    END IF

    _DELAY 0.003
    t = t + dt

LOOP WHILE press$ <> "q"

IF check = 0 THEN
    check = 1
    PRINT "Press q at any time to stop the simulation and see the bond angle."
    PRINT "Press any button to continue."
    SLEEP
    CLS
    GOTO afterBeg
END IF


END


FUNCTION Rounder (num, dec)
    Rounder = _ROUND(num * 10 ^ dec) / 10 ^ dec
END FUNCTION

' Credit to qb64 wiki for this function
FUNCTION ARCCOS (x) ' Inverse Cosine
    IF x < 1 THEN ARCCOS = (2 * ATN(1)) - ATN(x / SQR(1 - x * x))
END FUNCTION


SUB DefineProj
    xEye = 0
    yEye = -40
    zEye = 0
    yplane = -20
END SUB



FUNCTION getComp (dx, dy, dz, r, F, mode$) ' mode$ denotes which comp of force to get, made r a parameter bc I already calculated it in program, now using this for not forces, oops

    IF mode$ = "X" THEN
        comp = F * dx / r
    ELSEIF mode$ = "Y" THEN
        comp = F * dy / r
    ELSEIF mode$ = "Z" THEN
        comp = F * dz / r
    ELSE
        comp = 0
    END IF

    getComp = comp

END FUNCTION



FUNCTION RadiusProjection (y, z, zproj, rc)

    RadiusProjection = rc * SQR(((yplane - yEye) ^ 2 + (zproj - zEye) ^ 2) / ((y - yEye) ^ 2 + (z - zEye) ^ 2))

END FUNCTION



FUNCTION Projection (x, y, z, mode$)

    IF mode$ = "X" THEN
        proj = xEye + (yplane - yEye) * (x - xEye) / (y - yEye)
    ELSEIF mode$ = "Z" THEN
        proj = zEye + (yplane - yEye) * (z - zEye) / (y - yEye)
    ELSE
        proj = 0
    END IF

    Projection = proj

END FUNCTION



SUB InitializeAtoms (w, n) ' w is range of values num can be, from -w to w

    RANDOMIZE TIMER

    FOR i = 0 TO n - 1
        FOR j = 0 TO 2
            h(i, j) = INT(2 * w * (RND - 0.5))
        NEXT j
    NEXT i

    h(n, 0) = 0
    h(n, 1) = 0
    h(n, 2) = 0

    FOR i = 0 TO n
        q(i) = 1
        m(i) = 1
    NEXT i

END SUB



