import math
ener = 0.997;
vdwradius = 3.40;
r = 10.0;
i = r
with open('L-J.out', 'w') as f:
	while i > 1 :
		term1 = (math.pow(vdwradius/i, 12))
		term2 = (math.pow(vdwradius/i, 6))
		val = (term1-term2) * 4.0 * ener 
		f.write(i, val)
		f.close()
		i-=0.1
