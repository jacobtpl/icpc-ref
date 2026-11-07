/**
 * Author: Ulf Lundstrom
 * Date: 2009-04-07
 * License: CC0
 * Source: My geometric reasoning
 * Description: Returns the shortest distance on the sphere with radius radius between the points with azimuthal angles (longitude) f1 ($\phi_1$) and f2 ($\phi_2$) from x axis and zenith angles (latitude) t1 ($\theta_1$) and t2 ($\theta_2$) from z axis (0 = north pole). All angles measured in radians. The algorithm starts by converting the spherical coordinates to cartesian coordinates so if that is what you have you can use only the two last rows. d and s are the lengths of the difference and the sum of the two unit vectors; the angle between them is 2 atan2(d, s), which stays accurate (and never NaN) for nearby and antipodal points alike.
 * Status: tested on kattis:airlinehub
 */
#pragma once

double sphericalDistance(double f1, double t1,
		double f2, double t2, double radius) {
	double x1 = sin(t1)*cos(f1), y1 = sin(t1)*sin(f1);
	double x2 = sin(t2)*cos(f2), y2 = sin(t2)*sin(f2);
	double z1 = cos(t1), z2 = cos(t2);
	typedef long double ld; // no FMA contraction: exact 0 if equal
	ld dx = (ld)x2 - x1, dy = (ld)y2 - y1, dz = (ld)z2 - z1;
	ld ax = (ld)x2 + x1, ay = (ld)y2 + y1, az = (ld)z2 + z1;
	ld d = sqrtl(dx*dx + dy*dy + dz*dz);
	ld s = sqrtl(ax*ax + ay*ay + az*az);
	return (double)(radius*2*atan2l(d, s));
}
