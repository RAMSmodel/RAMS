!##############################################################################
Subroutine ll_xy (qlat,qlon,polelat,polelon,x,y)

implicit none

real :: qlat,qlon,polelat,polelon,x,y
real, parameter :: erad=6.367e6,erad2=1.2734e7,pi180=3.14159265/180.
real :: sinplat,cosplat,sinplon,cosplon,sinqlat,cosqlat,sinqlon,cosqlon &
       ,x3p,y3p,z3p, z3q,x3q,y3q, xq,yq,zq, t
       
! Evaluate sine and cosine of latitude and longitude of pole point p and
! input point q.

sinplat = sin(polelat * pi180)
cosplat = cos(polelat * pi180)
sinplon = sin(polelon * pi180)
cosplon = cos(polelon * pi180)

sinqlat = sin(qlat * pi180)
cosqlat = cos(qlat * pi180)
sinqlon = sin(qlon * pi180)
cosqlon = cos(qlon * pi180)

! Compute (x3,y3,z3) coordinates where the origin is the center of the earth,
! the z axis is the north pole, the x axis is the equator and prime
! meridian, and the y axis is the equator and 90 E.

! For the pole point, these are:

x3p = erad * cosplat * cosplon
y3p = erad * cosplat * sinplon
z3p = erad * sinplat

! For the given lat,lon point, these are:

z3q = erad * sinqlat
x3q = erad * cosqlat * cosqlon
y3q = erad * cosqlat * sinqlon

! Transform q point from (x3,y3,z3) coordinates in the above grid to
! polar stereographic coordinates (x,y,z):

xq = - sinplon * (x3q-x3p) + cosplon * (y3q-y3p)
yq =   cosplat * (z3q-z3p)  &
     - sinplat * ( cosplon * (x3q-x3p) + sinplon * (y3q-y3p) )
zq =   sinplat * (z3q-z3p)  &
     + cosplat * ( cosplon * (x3q-x3p) + sinplon * (y3q-y3p) )

! Parametric equation for line from antipodal point at (0,0,-2 erad) to
! point q has the following parameter (t) value on the polar stereographic
! plane:

t = erad2 / (erad2 + zq)

! This gives the following x and y coordinates for the projection of point q
! onto the polar stereographic plane:

x = xq * t
y = yq * t

return
END SUBROUTINE ll_xy

!##############################################################################
Subroutine xy_ll (qlat,qlon,polelat,polelon,x,y)

implicit none

real :: qlat,qlon,polelat,polelon,x,y
real, parameter :: erad=6.367e6,erad2=1.2734e7,pi180=3.14159265/180.
real :: sinplat,cosplat,sinplon,cosplon &
       ,x3p,y3p,z3p, z3q,x3q,y3q, xq,yq,zq, t,d,alpha,r3q

! Evaluate sine and cosine of latitude and longitude of pole point p.

sinplat = sin(polelat * pi180)
cosplat = cos(polelat * pi180)
sinplon = sin(polelon * pi180)
cosplon = cos(polelon * pi180)

! Compute (x3,y3,z3) coordinates of the pole point where the origin is the
! center of the earth, the z axis is the north pole, the x axis is the
! equator and prime meridian, and the y axis is the equator and 90 E.

x3p = erad * cosplat * cosplon
y3p = erad * cosplat * sinplon
z3p = erad * sinplat

! Compute distance d from given point R on the polar stereographic plane
! to the pole point P:

d = sqrt (x ** 2 + y ** 2)

! Compute angle QCP where C is the center of the Earth.  This is twice
! angle QAP where A is the antipodal point.  Angle QAP is the same as
! angle RAP:

alpha = 2. * atan2(d,erad2)

! Compute zq, the height of Q relative to the polar stereographic plane:

zq = erad * (cos(alpha) - 1.)

! Compute the parameter t which is the the distance ratio AQ:AR

t = (erad2 + zq) / erad2

! Compute xq and yq, the x and y coordinates of Q in polar stereographic space:

xq = t * x
yq = t * y

! Transform location of Q from (x,y,z) coordinates to (x3,y3,z3):

x3q = x3p - xq * sinplon - yq * cosplon * sinplat  &
    + zq * cosplat * cosplon
y3q = y3p + xq * cosplon - yq * sinplon * sinplat  &
    + zq * cosplat * sinplon
z3q = z3p + yq * cosplat + zq * sinplat

! Compute the latitude and longitude of Q:

qlon = atan2(y3q,x3q) / pi180
r3q = sqrt(x3q ** 2 + y3q ** 2)
qlat = atan2(z3q,r3q) / pi180

return
END SUBROUTINE xy_ll

!##############################################################################
Subroutine uevetouv (u,v,ue,ve,qlat,qlon,polelat,polelon)

implicit none

real :: u,v,ue,ve,qlat,qlon,polelat,polelon
real :: angle,x0,y0,x1,y1

CALL ll_xy (qlat,qlon,polelat,polelon,x0,y0)
CALL ll_xy (qlat,qlon+.1,polelat,polelon,x1,y1)
angle = atan2(y1-y0,x1-x0)
u = ue * cos(angle) - ve * sin(angle)
v = ue * sin(angle) + ve * cos(angle)

return
END SUBROUTINE uevetouv

!##############################################################################
Subroutine uvtoueve (u,v,ue,ve,qlat,qlon,polelat,polelon)

implicit none

real :: u,v,ue,ve,qlat,qlon,polelat,polelon
real :: angle,x0,y0,x1,y1

CALL ll_xy (qlat,qlon,polelat,polelon,x0,y0)
CALL ll_xy (qlat,qlon+.1,polelat,polelon,x1,y1)
angle = -atan2(y1-y0,x1-x0)
ue = u * cos(angle) - v * sin(angle)
ve = u * sin(angle) + v * cos(angle)

return
END SUBROUTINE uvtoueve

!##############################################################################
Subroutine winduv (dd,ff,uu,vv)

implicit none

real :: dd,ff,uu,vv
real, parameter :: pi=3.14159,pi180=pi/180.,v180pi=180./pi

uu=-ff*sin(dd*pi180)
vv=-ff*cos(dd*pi180)
return
END SUBROUTINE winduv

!##############################################################################
Subroutine winddf (dd,ff,uu,vv)

implicit none

real :: dd,ff,uu,vv
real, parameter :: pi=3.14159,pi180=pi/180.,v180pi=180./pi
real :: u,v

u=uu
v=vv
ff=sqrt(u*u+v*v)
if(abs(u).lt.1.e-20)u=1.e-20
if(abs(v).lt.1.e-20)v=1.e-20
dd=atan2(-u,-v)*v180pi
if(dd < 0.)dd=dd+360.

return
END SUBROUTINE winddf
