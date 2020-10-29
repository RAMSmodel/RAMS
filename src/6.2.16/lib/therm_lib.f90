!##############################################################################
real Function rslf (p,t)

!     This func calculates the liquid saturation vapor mixing ratio as
!     a func of pressure and Kelvin temperature

implicit none

real esl,x,t,p,c0,c1,c2,c3,c4,c5,c6,c7,c8

parameter (c0= .6105851e+03,c1= .4440316e+02,c2= .1430341e+01)
parameter (c3= .2641412e-01,c4= .2995057e-03,c5= .2031998e-05)
parameter (c6= .6936113e-08,c7= .2564861e-11,c8=-.3704404e-13)

x    = max(-80.,t-273.15)
esl  = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
rslf = .622*esl/(p-esl)

return
END FUNCTION rslf

!##############################################################################
real Function rsif (p,t)

!     This func calculates the ice saturation vapor mixing ratio as a
!     func of pressure and Kelvin temperature

implicit none

real esi,x,t,p,c0,c1,c2,c3,c4,c5,c6,c7,c8

parameter (c0= .6114327e+03,c1= .5027041e+02,c2= .1875982e+01)
parameter (c3= .4158303e-01,c4= .5992408e-03,c5= .5743775e-05)
parameter (c6= .3566847e-07,c7= .1306802e-09,c8= .2152144e-12)

x    = max(-80.,t-273.15)
esi  = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
rsif = .622*esi/(p-esi)

return
END FUNCTION rsif

!##############################################################################
real Function rslif (p,t)

implicit none

real :: p,t
real, external :: rslf
real, external :: rsif

!     This func calculates the saturation vapor mixing ratio, over
!     liquid or ice depending on temperature, as a func of pressure 
!     and Kelvin temperature

if (t >= 273.15) then
   rslif = rslf(p,t)
else
   rslif = rsif(p,t)
endif

return
END FUNCTION rslif

!##############################################################################
real Function eslf (t)

!     This func calculates the liquid saturation vapor pressure as a
!     func of Celcius temperature

implicit none

real :: x,t

real, parameter ::c0= .6105851e+03,c1= .4440316e+02,c2= .1430341e+01
real, parameter ::c3= .2641412e-01,c4= .2995057e-03,c5= .2031998e-05
real, parameter ::c6= .6936113e-08,c7= .2564861e-11,c8=-.3704404e-13

x    = max(-80.,t)
eslf = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

return
END FUNCTION eslf

!##############################################################################
real Function esif (t)

!     This func calculates the ice saturation vapor pressure as a
!     func of Celsius temperature

implicit none

real :: x,t

real, parameter ::c0= .6114327e+03,c1= .5027041e+02,c2= .1875982e+01
real, parameter ::c3= .4158303e-01,c4= .5992408e-03,c5= .5743775e-05
real, parameter ::c6= .3566847e-07,c7= .1306802e-09,c8= .2152144e-12

x    = max(-80.,t)
esif = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

return
END FUNCTION esif

!##############################################################################
real Function eslpf (t)

!     This func calculates the partial derivative of liquid saturation vapor
!     pressure with respect to temperature as a func of Celsius temperature

implicit none

real :: x,t

real, parameter ::d0= .4443216e+02,d1= .2861503e+01,d2= .7943347e-01
real, parameter ::d3= .1209650e-02,d4= .1036937e-04,d5= .4058663e-07
real, parameter ::d6=-.5805342e-10,d7=-.1159088e-11,d8=-.3189651e-14

x     = max(-80.,t)
eslpf = d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))

return
END FUNCTION eslpf

!##############################################################################
real Function esipf (t)

!     This func calculates the partial derivative of ice saturation vapor
!     pressure with respect to temperature as a func of Celsius temperature

implicit none

real  :: x,t

real, parameter ::d0= .5036342e+02,d1= .3775758e+01,d2= .1269736e+00
real, parameter ::d3= .2503052e-02,d4= .3163761e-04,d5= .2623881e-06
real, parameter ::d6= .1392546e-08,d7= .4315126e-11,d8= .5961476e-14

x     = max(-80.,t)
esipf = d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))

return
END FUNCTION esipf

!##############################################################################
Subroutine mrsl (n1,p,t,rsl)

implicit none

integer :: n,n1
real :: rsl(n1),p(n1),t(n1)
real, external :: rslf

do n=1,n1
   rsl(n) = rslf(p(n),t(n))
enddo

return
END SUBROUTINE mrsl

!##############################################################################
Subroutine mrsi (n1,p,t,rsi)

implicit none

integer :: n,n1
real :: rsi(n1),p(n1),t(n1)
real, external :: rsif

do n=1,n1
   rsi(n) = rsif(p(n),t(n))
enddo

return
END SUBROUTINE mrsi

!##############################################################################
real Function tdewpt (p,rs)

implicit none

real :: rr,rs,es,esln,p

rr=rs+1e-8
es=p*rr/(.622+rr)
esln=log(es)
tdewpt = (35.86*esln-4947.2325)/(esln-23.6837)

return
END FUNCTION tdewpt

!##############################################################################
real Function rsatmix (p,t)

implicit none

real :: p,t,es

es = 610.78*exp(17.269*(t-273.15)/(t-35.86))
rsatmix = .622*es/(p-es)

return
END FUNCTION rsatmix

!##############################################################################
Subroutine thetae (p,t,rv,the)

implicit none

real :: p,t,rv,the
real, parameter :: cp=1004.,g=9.8,r=287.,alvl=2.35e6,cpg=cp/g
real :: pit,tupo,ttd,dz,tupn,tmn
integer :: itter
real, external :: tdewpt

pit=p
tupo=t
ttd=tdewpt(p,rv)
dz=cpg*(t-ttd)
if(dz.le.0.) goto 20
do itter=1,50
   tupn=t-g/cp*dz
   tmn=(tupn+t)*.5*(1.+.61*rv)
   pit=p*exp(-g*dz/(r*tmn))
   if(abs(tupn-tupo).lt.0.001) goto 20
   ttd=tdewpt(pit,rv)
   tupo=tupn
   dz=dz+cpg*(tupn-ttd)
enddo
stop 10
20 continue
the=tupo*(1e5/pit)**.286*exp(alvl*rv/(cp*tupo))

return
END SUBROUTINE thetae

!##############################################################################
Subroutine the2t (the,p,th,t,r)

implicit none

real :: the,p,th,t,r
real, parameter :: cp=1004.,alvl=2.350e6
real :: pi,to,tn
integer :: itter
real, external :: rsatmix

pi=(p*1e-5)**.286
to=the/exp(alvl*.012/(cp*295.))*pi
do itter=1,50
   r=rsatmix(p,to)
   th=the/exp(alvl*r/(cp*to))
   tn=th*pi
   if(abs(to-tn).lt.0.005) goto 12
   to=to+(tn-to)*.3
enddo
write(6,1) the,p,to,tn,th,r
1 format(' stop in routine the2t '/' the,p,to,tn,th,r',6e15.6)
stop 10
12 continue
t=tn

return
END SUBROUTINE the2t

!##############################################################################
Subroutine qtk (q,tempk,fracliq)

implicit none

real :: q,tempk,fracliq
real,parameter :: r4186=1./4186.,r2093=1./2093.,r334000=1./334000.

!     Input:
!        q        internal energy [J/kg]
!     Outputs:
!       tempk    temperature [K]
!       fracliq  liquid fraction [dimensionless]
!     Local Constants:
!       4186     specific heat of liquid [J/(kg K)]
!       2093     specific heat of ice [J/(kg K)]
!       334000   latent heat of fusion [J/kg]
!       273.15   conversion from temp [C] to temp [K]

if (q .le. 0.) then
   fracliq = 0.
   tempk = q * r2093 + 273.15
elseif (q .ge. 334000.) then
   fracliq = 1.
   tempk = q * r4186 + 193.36
else
   fracliq = q * r334000
   tempk = 273.15
endif

return
END SUBROUTINE qtk

!##############################################################################
Subroutine qtc (q,tempc,fracliq)

implicit none

real :: q,tempc,fracliq
real,parameter :: r4186=1./4186.,r2093=1./2093.,r334000=1./334000.

!     Input:
!        q        internal energy [J/kg]
!     Outputs:
!        tempc    temperature [C]
!        fracliq  liquid fraction [dimensionless]
!     Local Constants:
!        4186     specific heat of liquid [J/(kg K)]
!        2093     specific heat of ice [J/(kg K)]
!        334000   latent heat of fusion [J/kg]
!        273.15   conversion from temp [C] to temp [K]

if (q .le. 0.) then
   fracliq = 0.
   tempc = q * r2093
elseif (q .ge. 334000.) then
   fracliq = 1.
   tempc = q * r4186 - 80.
else
   fracliq = q * r334000
   tempc = 0.
endif

return
END SUBROUTINE qtc

!##############################################################################
Subroutine qwtk (qw,w,dryhcap,tempk,fracliq)

implicit none

real :: qw,w,dryhcap,tempk,fracliq
real, parameter :: r4186=1./4186.,r2093=1./2093.,r334000=1./334000.
real :: qwliq0

!     Inputs:
!        qw       internal energy [J/m^2] or [J/m^3]
!        w        mass [kg/m^2] or [kg/m^3]
!        dryhcap  heat capacity of nonwater part [J/(m^2 K)] or [J/(m^3 K)]
!     Outputs:
!        tempk    temperature [K]
!        fracliq  liquid fraction [dimensionless]
!     Local Constants:
!        4186     specific heat of liquid [J/(kg K)]
!        2093     specific heat of ice [J/(kg K)]
!        334000   latent heat of fusion [J/kg]
!        273.15   conversion from temp [C] to temp [K]

qwliq0 = w * 334000.
if (qw .le. 0.) then
   fracliq = 0.
   tempk = qw / (2093. * w + dryhcap) + 273.15
elseif (qw .ge. qwliq0) then
   fracliq = 1.
   tempk = (qw - qwliq0) / (4186. * w + dryhcap) + 273.15
else
   fracliq = qw / qwliq0
   tempk = 273.15
endif

return
END SUBROUTINE qwtk
