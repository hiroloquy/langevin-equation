reset

#=================== Parameter ====================
m = 1
gamma = 1
r = 1
dt = 1.0e-2		# Time step
dh = dt/6.0		# coefficient for RK4
dt2 = dt/2.0	# coefficient for RK4

t_fin = 100
limit = ceil(t_fin / dt)	# Number of calculation

skip = 40					# skip frame

R = 5e1

#=================== Function ====================
# Round off to the i decimal place.
round(x, i) = 1 / (10.**(i+1)) * floor(x * (10.**(i+1)) + 0.5)

# Langevin equation (2D)
f1(x, y, vx, vy, t) = vx	# dx/dt
f2(x, y, vx, vy, t) = vy	# dy/dt
f3(x, y, vx, vy, t) = - gamma/m * vx + R * invnorm(rand(0))		# dvx/dt (white Gaussian noise)
f4(x, y, vx, vy, t) = - gamma/m * vy + R * invnorm(rand(0))		# dvy/dt

# Runge-Kutta 4th (Define rk_i(x, y, vx, vy, t))
do for[i=1:4]{
	rki = "rk"
	fi  = "f".sprintf("%d", i)
	rki = rki.sprintf("%d(x, y, vx, vy, t) = (\
	    k1 = %s(x, y, vx, vy, t),\
	    k2 = %s(x + dt2*k1, y + dt2*k1, vx + dt2*k1, vy + dt2*k1, t + dt2),\
	    k3 = %s(x + dt2*k2, y + dt2*k2, vx + dt2*k2, vy + dt2*k2, t + dt2),\
	    k4 = %s(x + dt*k3, y + dt*k3, vx + dt*k3, vy + dt*k3, t + dt),\
	    dh * (k1 + 2*k2 + 2*k3 + k4))", i, fi, fi, fi, fi)
	eval rki
}

# Time t
Time(t) = sprintf("{/TimesNewRoman:Italic t} = %.1f s", t) 

#========== Calculation ==========
# Prepare DAT file
outputfile = "langevin_data_gif.dat"
set print outputfile

print sprintf("# dt=%.2f, m=%.2f, gamma=%.2f, R=%.2f", dt, m, gamma, R)	# parameters
print "# t\t x\t y\t vx\t vy"	# items


# t = 0 (Initial value)
x	=	0		# position
y	=	0.0
vx	=	10		# velocity
vy	=	10
t	=	0.0		# time

print round(t, 2), round(x, 2), round(y, 2), round(vx, 2), round(vy, 2)	# in outputfile

# t > 0
do for [i = 1:limit] {
    t  = t  + dt
    x  = x  + rk1(x, y, vx, vy, t)
    y  = y  + rk2(x, y, vx, vy, t)
    vx = vx + rk3(x, y, vx, vy, t)
    vy = vy + rk4(x, y, vx, vy, t)

    print round(t, 2), round(x, 2), round(y, 2), round(vx, 2), round(vy, 2)	# in outputfile
}

unset print
print "Finish calculation!"

#=================== Animation ====================
set term gif animate delay 8 size 960, 720
set output 'Langevin_eq.gif'

set grid
set nokey
L = 1.5e2
set xr[-L/2:L/2]
set yr[-L/2:L/2]
set xlabel '{/TimesNewRoman:Italic=22 x}'
set ylabel '{/TimesNewRoman:Italic=22 y}'
set size ratio -1

do for [i = 0:limit:skip] { 	#[start:end:increment]
	t = i * dt
	set title Time(t) font 'TimesNewRoman:Normal, 22'

	plot outputfile u 2:3 every ::i::i w p pt 7 ps 2 lc rgb 'royalblue' title "", \
		 "" u 2:3 every ::::i w l ls 4 lc rgb 'royalblue' title ""
}

set out
print "Finish animation!"